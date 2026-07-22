#!/usr/bin/env python3
"""Run the primary/secondary adaptive-bias proof on one remoll ROOT file.

The current proof fully exercises the secondary branch: recover stored event
histories, record real parent states entering the responsible physical volume,
and compare a physical-bank replay with a parent-neighborhood replay.  No
result is merged into the original ROOT file.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import random
import shutil
import subprocess
import sys
from array import array
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = next(path for path in (HERE, *HERE.parents)
            if (path / "geometry/mollerMother.gdml").is_file())
ANALYZE = HERE / "analyze_surface_scores.C"
DIAGNOSE = HERE / "diagnose_tile.C"
DEFAULT_INPUT = REPO / (
    "_local_additions_archive/calibration/neff_phase_space/"
    "neff_carbon_US_2p2_100k_sieveout_sieve_out.root"
)
DEFAULT_OUTPUT = REPO / (
    "_local_additions_archive/background_campaign/runs/us_2p2_twofold_poc"
)


def run(command: list[str], *, cwd: Path = REPO, log: Path | None = None) -> str:
    if log is None:
        result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
        output = result.stdout + result.stderr
    else:
        with log.open("w") as stream:
            result = subprocess.run(command, cwd=cwd, text=True,
                                    stdout=stream, stderr=subprocess.STDOUT)
        output = log.read_text(errors="replace")
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n"
                           f"{output[-5000:]}")
    return output


def root_call(macro: Path, arguments: list[str]) -> str:
    expression = f"{macro}({','.join(arguments)})"
    return run([str(REPO / "build/reroot"), "-l", "-b", "-q", expression])


def quote(value: str | Path) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def analyze(path: Path, mode: str, scale: int, output: Path) -> dict[str, dict]:
    root_call(ANALYZE, [quote(path), quote(mode), str(scale), quote(output)])
    result = {}
    with output.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            n = int(row["events"])
            total = float(row["sum"])
            total2 = float(row["sum2"])
            variance = max(0.0, (total2 - total * total / n) / (n - 1))
            mean = total / n
            result[row["key"]] = {
                "events": n,
                "signal": int(row["signal"]),
                "estimate": mean,
                "rse": math.sqrt(variance / n) / abs(mean) if mean else math.inf,
                "neff": total * total / total2 if total2 else 0.0,
                "maximum_fraction": float(row["maximum"]) / abs(total) if total else math.inf,
            }
    return result


def choose_worst(stats: dict[str, dict]) -> tuple[str, dict]:
    candidates = [(value["rse"], key, value) for key, value in stats.items()
                  if key.startswith("tile_")]
    _, key, value = max(candidates)
    return key, value


def diagnose(input_path: Path, output_dir: Path, key: str) -> tuple[dict, dict]:
    tile, observable = key.split("/")
    detector = int(tile.removeprefix("tile_"))
    prefix = output_dir / f"{tile}_{observable}"
    root_call(DIAGNOSE, [quote(input_path), str(detector), quote(observable),
                         quote(prefix), "1000", "100"])
    histories = list(csv.DictReader(
        (prefix.with_suffix(".histories.tsv")).open(), delimiter="\t"))
    hits = list(csv.DictReader(
        (prefix.with_suffix(".hits.tsv")).open(), delimiter="\t"))
    if not histories or not hits:
        raise RuntimeError("worst tile has no ranked history/hit details")
    return histories[0], next(hit for hit in hits
                              if hit["entry"] == histories[0]["entry"])


def export_states(input_path: Path, original_entries: list[int], directory: Path) -> None:
    import ROOT
    directory.mkdir(parents=True, exist_ok=True)
    source = ROOT.TFile.Open(str(input_path), "READ")
    tree = source.Get("T")
    run_leaf = tree.GetLeaf("seed.fRunNo")
    event_leaf = tree.GetLeaf("seed.fEvtNo")
    if run_leaf is None or event_leaf is None:
        raise RuntimeError("input ROOT file does not expose stored seed identities")
    seed_identities = []
    for entry in original_entries:
        tree.GetEntry(entry)
        seed_identities.append((int(run_leaf.GetValue()), int(event_leaf.GetValue())))
    source.Close()
    selection = "||".join(f"Entry$=={entry}" for entry in original_entries)
    command = (
        f'TFile f({quote(input_path)}); auto* T=(TTree*)f.Get("T"); '
        f'T->Draw("seed.Save()","{selection}","goff");'
    )
    run([str(REPO / "build/reroot"), "-l", "-b", "-e", command], cwd=directory)
    for replay_event, (run_number, event_number) in enumerate(seed_identities):
        state = directory / f"run{run_number}evt{event_number}.state"
        if not state.is_file():
            raise RuntimeError(f"stored RNG state was not exported: {state}")
        values = [line.strip() for line in state.read_text().splitlines()
                  if line.strip() and "MixMaxRng" not in line]
        if len(values) < 3:
            raise RuntimeError(f"invalid MixMax state: {state}")
        vector, counter, sum_total = values[:-2], values[-2], values[-1]
        converted = (
            "mixmax state, file version 1.0\n"
            f"N={len(vector)}; V[N]={{{', '.join(vector)}}}; "
            f"counter={counter}; sumtot={sum_total};\n"
        )
        (directory / f"run0evt{replay_event}.rndm").write_text(converted)


def base_setup() -> list[str]:
    return [
        f"/remoll/setgeofile {REPO / 'geometry/mollerMother.gdml'}",
        "/remoll/physlist/register FTFP_BERT",
        "/remoll/physlist/parallel/enable",
        f"/remoll/parallel/setfile {REPO / 'geometry/mollerParallel.gdml'}",
        "/run/initialize",
        f"/remoll/addfield {REPO / 'map_directory/V2U.1a.50cm.parallel.txt'}",
        f"/remoll/addfield {REPO / 'map_directory/DS_TM1-4_CoilA-G_ll_TM2-4_out3mm.txt'}",
    ]


def detector_setup(detector: int) -> list[str]:
    return [
        "/remoll/SD/disable_all",
        f"/remoll/SD/enable {detector}",
        f"/remoll/SD/detect secondaries {detector}",
        f"/remoll/SD/detect boundaryhits {detector}",
        f"/remoll/SD/detect lowenergyneutral {detector}",
    ]


def write_exact_macro(path: Path, result: Path, interactions: Path,
                      detector: int, volume: str, events: int) -> None:
    lines = base_setup() + [
        "/remoll/geometry/absolute_position targetLadder (0,500,0)",
        "/remoll/target/mother Optics1",
        "/remoll/target/volume DSC",
        "/remoll/geometry/relative_position USsieve (200,0,0)",
        "/remoll/kryptonite/volume donutSystem_logical",
        "/remoll/kryptonite/volume logic_DSpipe_vacuumtube",
        "/remoll/evgen/set elasticC12",
        "/remoll/oldras false",
        "/remoll/beam_corrph 0.02134",
        "/remoll/beam_corrth 0.02134",
        "/remoll/evgen/thmin 0.1 deg",
        "/remoll/evgen/thmax 1.25 deg",
        "/remoll/evgen/phmin 0 deg",
        "/remoll/evgen/phmax 360 deg",
        "/remoll/beamene 2.2 GeV",
        "/remoll/beamcurr 1 microampere",
    ] + detector_setup(detector) + [
        f"/remoll/interaction/output {interactions}",
        f"/remoll/interaction/physicalVolume {volume}",
        "/remoll/interaction/includeTransportation true",
        "/remoll/interaction/enable true",
        "/random/resetEngineFromEachEvent 1",
        f"/remoll/filename {result}",
        f"/run/beamOn {events}",
    ]
    path.write_text("\n".join(lines) + "\n")


def read_event_rates(root_path: Path) -> list[float]:
    import ROOT
    source = ROOT.TFile.Open(str(root_path), "READ")
    tree = source.Get("T")
    rates = []
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        rates.append(float(tree.rate))
    source.Close()
    return rates


def entry_states(table: Path, volume: str, rates: list[float]) -> list[dict]:
    states = []
    with table.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            event = int(row["event"])
            if (row["role"] != "continuation" or row["post_volume"] != volume
                    or row["pre_volume"] == volume or int(row["pdg"]) != 11
                    or float(row["post_pz_mev"]) <= 0.0 or event >= len(rates)):
                continue
            states.append({
                "event": event, "interaction": int(row["interaction"]),
                "track": int(row["track"]), "pid": 11,
                "rate": rates[event] * float(row["weight"]),
                "x": float(row["post_x_mm"]), "y": float(row["post_y_mm"]),
                "z": float(row["post_z_mm"]),
                "px": float(row["post_px_mev"]), "py": float(row["post_py_mev"]),
                "pz": float(row["post_pz_mev"]), "k": float(row["post_ke_mev"]),
            })
    return states


def resolve_parent_entry(table: Path, volume: str, immediate_parent: int,
                         target_hit: dict, states: list[dict]) -> dict:
    """Find the last lineage entry before the selected secondary is created."""
    with table.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream, delimiter="\t")
                if int(row["event"]) == 0]
    target_pid = int(target_hit["pid"])
    target_k = float(target_hit["k_mev"])
    target_xyz = tuple(float(target_hit[name]) for name in ("vx_mm", "vy_mm", "vz_mm"))
    target_process = target_hit["creator_process"]
    candidates = []
    for row in rows:
        if (int(row["track"]) != immediate_parent or row["role"] != "secondary"
                or int(row["pdg"]) != target_pid or row["process"] != target_process):
            continue
        xyz = tuple(float(row[name]) for name in ("post_x_mm", "post_y_mm", "post_z_mm"))
        position_error = math.sqrt(sum((a - b) ** 2 for a, b in zip(xyz, target_xyz)))
        if position_error < 1e-3 and abs(float(row["post_ke_mev"]) - target_k) < 1e-6:
            candidates.append(int(row["interaction"]))
    if not candidates:
        raise RuntimeError("exact replay did not recover the selected secondary creation step")
    creation_interaction = min(candidates)
    parents = {}
    for row in rows:
        parents.setdefault(int(row["track"]), int(row["parent"]))

    def find_creation_interaction(child: int, parent: int) -> int:
        child_rows = [row for row in rows
                      if int(row["track"]) == child and row["role"] == "continuation"]
        if not child_rows:
            raise RuntimeError(f"no recorded first state for ancestry track {child}")
        first = min(child_rows, key=lambda row: int(row["interaction"]))
        child_pid = int(first["pdg"])
        child_state = tuple(float(first[name]) for name in
                            ("pre_x_mm", "pre_y_mm", "pre_z_mm", "pre_ke_mev"))
        matches = []
        for row in rows:
            if (int(row["track"]) != parent or row["role"] != "secondary"
                    or int(row["pdg"]) != child_pid):
                continue
            produced = tuple(float(row[name]) for name in
                             ("post_x_mm", "post_y_mm", "post_z_mm", "post_ke_mev"))
            position_error = math.sqrt(sum((produced[i] - child_state[i]) ** 2
                                           for i in range(3)))
            if position_error < 1e-3 and abs(produced[3] - child_state[3]) < 1e-6:
                matches.append(int(row["interaction"]))
        if not matches:
            raise RuntimeError(f"could not resolve creation of track {child} by {parent}")
        return min(matches)

    track = immediate_parent
    seen = set()
    for _ in range(1000):
        if track <= 0 or track in seen:
            break
        seen.add(track)
        eligible = [state for state in states
                    if state["event"] == 0 and state["track"] == track
                    and state["interaction"] < creation_interaction]
        if eligible:
            center = max(eligible, key=lambda state: state["interaction"])
            center["creation_interaction"] = creation_interaction
            center["immediate_parent"] = immediate_parent
            return center
        parent = parents.get(track, 0)
        if parent <= 0:
            break
        creation_interaction = find_creation_interaction(track, parent)
        track = parent
    raise RuntimeError(
        f"no parent/ancestor of track {immediate_parent} entered {volume} before "
        f"the selected secondary was created at interaction {creation_interaction}")


def feature_scales(states: list[dict]) -> tuple[float, float, float, float, float]:
    def spread(values: list[float], floor: float) -> float:
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / len(values)
        return max(floor, 0.35 * math.sqrt(variance))
    return (
        spread([s["x"] for s in states], 2.0),
        spread([s["y"] for s in states], 2.0),
        spread([math.log(max(s["k"], 1e-9)) for s in states], 0.05),
        spread([s["px"] / s["pz"] for s in states], 5e-4),
        spread([s["py"] / s["pz"] for s in states], 5e-4),
    )


def write_source(path: Path, states: list[dict], center: dict, draws: int,
                 targeted: bool, physical_fraction: float, seed: int) -> None:
    import ROOT
    ROOT.gSystem.Load(str(REPO / "build/libremoll.dylib"))
    masses = [abs(state["rate"]) for state in states]
    normalization = sum(masses)
    physical = [mass / normalization for mass in masses]
    sx, sy, sk, stx, sty = feature_scales(states)
    affinity = []
    for state in states:
        distance2 = (
            ((state["x"] - center["x"]) / sx) ** 2
            + ((state["y"] - center["y"]) / sy) ** 2
            + ((math.log(max(state["k"], 1e-9))
                - math.log(max(center["k"], 1e-9))) / sk) ** 2
            + ((state["px"] / state["pz"] - center["px"] / center["pz"]) / stx) ** 2
            + ((state["py"] / state["pz"] - center["py"] / center["pz"]) / sty) ** 2
        )
        affinity.append(math.exp(-0.5 * distance2))
    target_norm = sum(mass * value for mass, value in zip(masses, affinity))
    target_pdf = [mass * value / target_norm for mass, value in zip(masses, affinity)]
    if targeted:
        proposal = [physical_fraction * p + (1.0 - physical_fraction) * t
                    for p, t in zip(physical, target_pdf)]
    else:
        proposal = physical
    cumulative, running = [], 0.0
    for value in proposal:
        running += value
        cumulative.append(running)

    output = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("T", "Observed parent-volume-entry replay source")
    event = ROOT.remollEvent_t()
    hits = ROOT.std.vector("remollGenericDetectorHit_t")()
    output_rate = array("d", [0.0])
    source_event = array("i", [0])
    source_track = array("i", [0])
    proposal_pdf = array("d", [0.0])
    physical_pdf = array("d", [0.0])
    tree.Branch("rate", output_rate, "rate/D")
    tree.Branch("ev", event)
    tree.Branch("hit", hits)
    tree.Branch("source_event", source_event, "source_event/I")
    tree.Branch("source_track", source_track, "source_track/I")
    tree.Branch("proposal_pdf", proposal_pdf, "proposal_pdf/D")
    tree.Branch("physical_pdf", physical_pdf, "physical_pdf/D")
    generator = random.Random(seed)
    for _ in range(draws):
        index = bisect.bisect_left(cumulative, generator.random() * running)
        index = min(index, len(states) - 1)
        state = states[index]
        hit = ROOT.remollGenericDetectorHit_t()
        hit.det = 990035
        hit.pid = 11
        hit.trid = 1
        momentum = math.sqrt(state["px"] ** 2 + state["py"] ** 2 + state["pz"] ** 2)
        hit.x = state["x"] + 0.01 * state["px"] / momentum
        hit.y = state["y"] + 0.01 * state["py"] / momentum
        hit.z = state["z"] + 0.01 * state["pz"] / momentum
        hit.px, hit.py, hit.pz = state["px"], state["py"], state["pz"]
        hit.p, hit.k = momentum, state["k"]
        hits.clear()
        hits.push_back(hit)
        output_rate[0] = state["rate"] / (draws * proposal[index])
        source_event[0], source_track[0] = state["event"], state["track"]
        proposal_pdf[0], physical_pdf[0] = proposal[index], physical[index]
        tree.Fill()
    tree.Write()
    output.Close()


def write_replay_macro(path: Path, source: Path, result: Path,
                       detector: int, events: int, seed: int) -> None:
    lines = base_setup() + [
        "/remoll/evgen/set external",
        "/remoll/evgen/copyRate 1",
        f"/remoll/evgen/external/file {source}",
        "/remoll/evgen/external/detid 990035",
        "/remoll/evgen/external/startEvent 0",
    ] + detector_setup(detector) + [
        f"/remoll/seed {seed}",
        f"/remoll/filename {result}",
        f"/run/beamOn {events}",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bank-events", type=int, default=200)
    parser.add_argument("--replay-events", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--physical-fraction", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=7163022)
    args = parser.parse_args()
    args.input = args.input.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    original = analyze(args.input, "original", 0, args.output_dir / "original_stats.tsv")
    key, worst = choose_worst(original)
    top_history, top_hit = diagnose(args.input, args.output_dir, key)
    detector = int(key.split("/")[0].removeprefix("tile_"))
    observable = key.split("/")[1]
    print(f"Worst: {key}, RSE {worst['rse']:.2%}, Neff {worst['neff']:.2f}; "
          f"top entry {top_history['entry']} contributes "
          f"{float(top_history['fraction_of_total']):.2%}.")
    if int(top_hit["mtrid"]) == 0:
        raise RuntimeError(
            "the selected history is primary; hand it to campaign.py's full-support "
            "generator-neighborhood proposal (the secondary volume branch was not selected)")

    volume = top_hit["creator_volume"]
    parent_track = int(top_hit["mtrid"])
    target_entry = int(top_history["entry"])
    generator = random.Random(args.seed)
    population = list(range(original[key]["events"]))
    population.remove(target_entry)
    original_entries = [target_entry] + generator.sample(population, args.bank_events - 1)
    seed_dir = args.output_dir / "event_states"
    if seed_dir.exists():
        shutil.rmtree(seed_dir)
    export_states(args.input, original_entries, seed_dir)

    exact_root = args.output_dir / "exact_bank.root"
    interaction_table = args.output_dir / "volume_entries.tsv"
    exact_macro = args.output_dir / "exact_bank.mac"
    write_exact_macro(exact_macro, exact_root, interaction_table,
                      detector, volume, len(original_entries))
    run([str(REPO / "build/remoll"), "-t", str(args.threads), "-m", str(exact_macro)],
        cwd=seed_dir, log=args.output_dir / "exact_bank.log")
    rates = read_event_rates(exact_root)
    states = entry_states(interaction_table, volume, rates)
    center = resolve_parent_entry(interaction_table, volume, parent_track,
                                  top_hit, states)
    if len(states) < 10:
        raise RuntimeError(f"only {len(states)} forward electron entry states were recovered")
    print(f"Recovered lineage ancestor track {center['track']} at {volume} entry "
          f"before immediate parent {parent_track} created the selected secondary, and "
          f"{len(states)} real forward electron entry states in the physical bank.")

    results = {}
    for targeted, label in ((False, "physical"), (True, "targeted")):
        source = args.output_dir / f"{label}_source.root"
        result = args.output_dir / f"{label}_result.root"
        macro = args.output_dir / f"{label}.mac"
        write_source(source, states, center, args.replay_events, targeted,
                     args.physical_fraction, args.seed + (1 if targeted else 0))
        write_replay_macro(macro, source, result, detector,
                           args.replay_events, args.seed + (11 if targeted else 10))
        run([str(REPO / "build/remoll"), "-t", str(args.threads), "-m", str(macro)],
            log=args.output_dir / f"{label}.log")
        stats = analyze(result, "replay", args.replay_events,
                        args.output_dir / f"{label}_stats.tsv")
        results[label] = stats[key]

    print("\nEqual-cost parent-volume replay comparison")
    print("proposal   signal histories    Neff       RSE       estimate")
    for label in ("physical", "targeted"):
        value = results[label]
        print(f"{label:9s} {value['signal']:16d} {value['neff']:8.2f} "
              f"{value['rse']:9.2%} {value['estimate']:14.6g}")
    if results["physical"]["signal"]:
        print(f"Useful-history gain: "
              f"{results['targeted']['signal'] / results['physical']['signal']:.3g}x")
    print("This is a proposal-efficiency pilot over an observed physical entry bank; "
          "it is not yet a production normalization or a merged final estimate.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
