#!/usr/bin/env python3
"""Restartable 1k-event ShowerMax response campaign and live dashboard."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import queue
import re
import signal
import subprocess
import sys
import threading
import time
import tomllib
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ANALYZER = HERE / "analyze_response_batch.C"
DASHBOARD = HERE / "dashboard.html"
RESPONSE_DATA = REPO / "analysis/showermax/response_lookup/data"
TAIL = 201


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def load_config(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def resolve(path: str) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else REPO / candidate


def load_window(table: Path, cell: str, septant: int) -> dict[str, str]:
    with table.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream, delimiter="\t")
                if row["cell"] == cell and int(row["septant"]) == septant]
    if len(rows) != 1:
        raise RuntimeError(
            f"expected exactly one kinematics row for cell={cell!r}, septant={septant}; "
            f"found {len(rows)} in {table}")
    return rows[0]


def target_macro(target: str) -> str:
    return {
        "lh2": "macros/target/LH2.mac",
        "c12_us": "macros/target/Optics1.mac",
        "c12_ms": "macros/target/Optics3.mac",
        "c12_ds": "macros/target/Optics2.mac",
    }[target]


def generator(channel: str) -> str:
    return {
        "moller": "moller",
        "ep_elastic": "elastic",
        "ep_inelastic": "inelastic",
        "c12_elastic": "elasticC12",
        "c12_inelastic": "inelasticC12",
    }[channel]


def validate_proposal(config: dict[str, Any], row: dict[str, str]) -> None:
    """Reject malformed compact-table proposals before starting Geant4."""
    fraction = float(config["physical_fraction"])
    if not math.isfinite(fraction) or not 0.0 < fraction <= 1.0:
        raise RuntimeError("physical_fraction must be finite and in (0,1]")
    ranges = {
        "theta": ("theta_q01_deg", "theta_q99_deg"),
        "beamp": ("beamp_q01_gev", "beamp_q99_gev"),
        "vertexz": ("vertexz_q01_fraction", "vertexz_q99_fraction"),
    }
    if row["channel"] == "ep_inelastic":
        ranges["outgoinge"] = (
            "outgoinge_q01_fraction", "outgoinge_q99_fraction")
    for axis, (low_key, high_key) in ranges.items():
        low, high = float(row[low_key]), float(row[high_key])
        if not math.isfinite(low) or not math.isfinite(high):
            raise RuntimeError(
                f"{row['cell']} has non-finite {axis} proposal bounds")
        if high < low:
            raise RuntimeError(
                f"{row['cell']} has reversed {axis} proposal bounds: "
                f"{low} to {high}")
        if axis in {"vertexz", "outgoinge"} and not (
                0.0 <= low <= high <= 1.0):
            raise RuntimeError(
                f"{row['cell']} has invalid fractional {axis} bounds: "
                f"{low} to {high}")
        if axis == "beamp" and (low < 0.0 or high > float(row["energy_mev"])/1000):
            raise RuntimeError(
                f"{row['cell']} has beam-momentum bounds outside [0,Ebeam]: "
                f"{low} to {high} GeV")


def macro_text(config: dict[str, Any], row: dict[str, str], sieve: str,
               root_file: Path, seed: int) -> str:
    validate_proposal(config, row)
    base = resolve(config["base_macro"]).read_text().rstrip()
    gen = generator(row["channel"])
    events = int(config["batch_events"])
    axes = ["theta", "beamp", "vertexz"]
    if row["channel"] == "ep_inelastic":
        axes.append("outgoinge")
    # Every axis independently retains the configured physical component.
    # This keeps full joint support without diluting the targeted component.
    fraction = float(config["physical_fraction"])
    lines = [
        base, "", "# ShowerMax live campaign (analysis controls only)",
        f"/control/execute {target_macro(row['target'])}",
    ]
    if sieve != "none":
        lines.append(f"/control/execute macros/sieve/sieve_{sieve}.mac")
    lines += [
        f"/remoll/evgen/set {gen}",
        f"/remoll/beamene {float(row['energy_mev']) / 1000.0:.12g} GeV",
        f"/remoll/beamcurr {float(config['current_uA']):.12g} microampere",
    ]
    if row["channel"] == "moller":
        lines += ["/remoll/evgen/thcommin 30 deg",
                  "/remoll/evgen/thcommax 150 deg"]
    else:
        theta_max = 3.0 if row["channel"] == "ep_elastic" else 5.0
        lines += ["/remoll/evgen/thmin 0.1 deg",
                  f"/remoll/evgen/thmax {theta_max:.12g} deg"]
    lines += ["/remoll/evgen/phmin 0 deg", "/remoll/evgen/phmax 360 deg"]
    if row["channel"] == "ep_elastic":
        lines.append("/remoll/evgen/elastic/applyScreening true")
    commands = {
        # SampleThetaWithBias is shared by all supported generators. Its
        # historical public messenger name is thcom even when the sampled
        # variable is a laboratory scattering angle.
        "theta": ("/remoll/bias/thcom", "min", "max",
                  row["theta_q01_deg"], row["theta_q99_deg"], "deg"),
        "beamp": ("/remoll/bias/beamp", "min", "max",
                  row["beamp_q01_gev"], row["beamp_q99_gev"], "GeV"),
        "vertexz": ("/remoll/bias/vertexz", "minFraction", "maxFraction",
                    row["vertexz_q01_fraction"], row["vertexz_q99_fraction"], ""),
        "outgoinge": ("/remoll/bias/outgoinge", "minFraction", "maxFraction",
                      row["outgoinge_q01_fraction"], row["outgoinge_q99_fraction"], ""),
    }
    for axis in axes:
        prefix, low_name, high_name, low, high, unit = commands[axis]
        # A compact-table quantile window can collapse for a monochromatic or
        # otherwise fixed variable (notably beam momentum in carbon Møller
        # cells). A uniform proposal over zero width is undefined. Leaving that
        # axis in its physical mode is exact and loses no targeted support.
        if float(high) <= float(low):
            lines.append(
                f"# {prefix} left physical: compact window is degenerate "
                f"({low} to {high} {unit})".rstrip())
            continue
        lines += [
            f"{prefix}/mode mixture",
            f"{prefix}/{low_name} {low} {unit}".rstrip(),
            f"{prefix}/{high_name} {high} {unit}".rstrip(),
            f"{prefix}/physicalFraction {fraction:.12g}",
        ]
    lines += [
        f"/remoll/interval {max(1, events // 20)}",
        "/remoll/SD/disable_all",
        "/remoll/SD/enable 30",
        "/remoll/SD/detect surfacehits 30",
        "/remoll/SD/detect lowenergyneutral 30",
    ]
    # G4GenericMessenger's G4TwoVector parser is unreliable for the range
    # commands on this build ("istream ended before second value"). Configure
    # each response plane explicitly so a silently empty run is impossible.
    for detector in range(73001, 73029):
        lines += [
            f"/remoll/SD/enable {detector}",
            f"/remoll/SD/detect surfacehits {detector}",
            f"/remoll/SD/detect lowenergyneutral {detector}",
        ]
    ring_suffixes = (10, 20, 30, 40, 51, 52, 53, 60)
    for base in (110000, 150000):
        for segment in range(1, 15):
            for suffix in ring_suffixes:
                detector = base + 100 * segment + suffix
                lines += [
                    f"/remoll/SD/enable {detector}",
                    f"/remoll/SD/detect surfacehits {detector}",
                    f"/remoll/SD/detect lowenergyneutral {detector}",
                ]
    lines += [
        f"/remoll/seed {seed}",
        f"/remoll/filename {root_file}",
        f"/run/beamOn {events}", "",
    ]
    return "\n".join(lines)


def parse_analysis(output: str, elapsed: float) -> dict[str, Any]:
    result: dict[str, Any] = {
        "elapsed_seconds": elapsed, "statistics": {}, "tiles": {}, "heatmap": {},
        "species": {}, "energy_spectra": {}, "global_maps": {},
        "global_hits": {}}
    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "LIVE_ERROR":
            raise RuntimeError("response analyzer: " + " ".join(fields[1:]))
        if fields[0] == "LIVE_BATCH":
            result["events"] = int(fields[1])
            result["unsupported_hits"] = int(fields[2])
            result["file_events"] = int(fields[3])
            result["start_entry"] = int(fields[4])
        elif fields[0] == "LIVE_STAT":
            key = f"{fields[1]}/{fields[2]}"
            result["statistics"][key] = {
                "detector": int(fields[3]), "signal_events": int(fields[4]),
                "sum": float(fields[5]), "sum2": float(fields[6]),
                "sum3": float(fields[7]), "sum4": float(fields[8]),
                "maximum": float(fields[9]),
                "largest": [] if fields[10] == "-" else
                    [float(item) for item in fields[10].split(",")],
            }
        elif fields[0] == "LIVE_TILE":
            result["tiles"][fields[1]] = {
                "detector": int(fields[2]), "hits": int(fields[3]),
                "response_rate": float(fields[4]),
            }
        elif fields[0] == "LIVE_BIN":
            result["heatmap"].setdefault(fields[1], []).append(
                [int(fields[2]), int(fields[3]), float(fields[4])])
        elif fields[0] == "LIVE_SPECIES":
            result["species"][fields[1]] = {
                "hits": int(fields[2]), "ignored": int(fields[3]),
                "rate": float(fields[4]), "ignored_rate": float(fields[5]),
            }
        elif fields[0] == "LIVE_ENERGY":
            result["energy_spectra"].setdefault(fields[1], []).append(
                [int(fields[2]), int(fields[3]), int(fields[4]),
                 float(fields[5]), float(fields[6])])
        elif fields[0] == "LIVE_GLOBAL":
            result["global_hits"][fields[1]] = int(fields[2])
        elif fields[0] == "LIVE_GLOBAL_BIN":
            result["global_maps"].setdefault(fields[1], []).append(
                [int(fields[2]), int(fields[3]), float(fields[4])])
    if "events" not in result:
        raise RuntimeError("ROOT helper produced no LIVE_BATCH record")
    return result


def analyze(root_file: Path, septant: int, elapsed: float, start: int = 0,
            events: int = -1) -> dict[str, Any]:
    reroot = next((path for path in
                   (REPO/"build-surface-replay/reroot", REPO/"build-develop/reroot",
                    REPO/"bin/reroot") if path.is_file()), None)
    if not reroot:
        raise RuntimeError("reroot executable was not found")
    invocation = (
        f'{ANALYZER}("{root_file}",{septant},"{RESPONSE_DATA}",{start},{events})')
    environment = os.environ.copy()
    if reroot.parent.name.startswith("build-"):
        environment["DYLD_LIBRARY_PATH"] = str(reroot.parent) + (
            ":" + environment["DYLD_LIBRARY_PATH"]
            if environment.get("DYLD_LIBRARY_PATH") else "")
    run = subprocess.run([str(reroot), "-l", "-b", "-q", invocation],
                         cwd=REPO, text=True, capture_output=True, env=environment)
    if run.returncode:
        raise RuntimeError("ROOT helper failed:\n" + run.stdout[-2000:] + run.stderr[-2000:])
    return parse_analysis(run.stdout + run.stderr, elapsed)


def slope(points: list[tuple[float, float]]) -> float | None:
    if len(points) < 3:
        return None
    xm = sum(x for x, _ in points) / len(points)
    ym = sum(y for _, y in points) / len(points)
    den = sum((x-xm)**2 for x, _ in points)
    return None if den == 0 else sum((x-xm)*(y-ym) for x, y in points) / den


def metrics(history: list[dict[str, Any]], key: str, cfg: dict[str, Any]) -> dict[str, Any]:
    settings = cfg["reliability"]
    total = {name: 0.0 for name in ("sum", "sum2", "sum3", "sum4")}
    events = nonzero = 0
    elapsed = 0.0
    maximum = 0.0
    largest: list[float] = []
    checkpoints = []
    batch_estimates = []
    for batch in history:
        value = batch["statistics"][key]
        n = batch["events"]
        events += n
        elapsed += batch["elapsed_seconds"]
        nonzero += value["signal_events"]
        for name in total:
            total[name] += value[name]
        maximum = max(maximum, value["maximum"])
        largest = sorted(largest + value["largest"], reverse=True)[:TAIL]
        mean = total["sum"] / events
        second = max(0.0, total["sum2"] - total["sum"]**2/events)
        variance = second/(events-1) if events > 1 else 0
        rse = math.sqrt(variance/events)/abs(mean) if mean and events > 1 else None
        fourth = max(0.0, total["sum4"] - 4*mean*total["sum3"]
                     + 6*mean*mean*total["sum2"] - 3*events*mean**4)
        vov = max(0.0, fourth/second**2 - 1/events) if second else None
        fom = 1/(rse*rse*elapsed) if rse and elapsed else None
        checkpoints.append({"events": events, "estimate": mean, "rse": rse,
                            "vov": vov, "fom": fom})
        bmean = value["sum"]/n
        bsecond = max(0.0, value["sum2"] - value["sum"]**2/n)
        batch_estimates.append((bmean, bsecond/(n*(n-1)) if n > 1 else 0))
    final = checkpoints[-1] if checkpoints else {}
    rse_slope = slope([(math.log(p["events"]), math.log(p["rse"]))
                       for p in checkpoints if p["rse"] not in (None, 0)])
    vov_slope = slope([(math.log(p["events"]), math.log(p["vov"]))
                       for p in checkpoints if p["vov"] not in (None, 0)])
    foms = [p["fom"] for p in checkpoints if p["fom"] is not None]
    recent = foms[len(foms)//2:]
    fom_range = ((max(recent)-min(recent))/(sum(recent)/len(recent))
                 if len(recent) >= 2 and sum(recent) > 0 else None)
    tail_index = None
    if len(largest) >= TAIL and largest[TAIL-1] > 0:
        logs = sum(math.log(x/largest[TAIL-1]) for x in largest[:TAIL-1])
        tail_index = (TAIL-1)/logs if logs > 0 else None
    next_fraction = None
    if events and total["sum"]:
        next_fraction = abs((total["sum"]+maximum)/(events+1)-total["sum"]/events) / \
                        abs(total["sum"]/events)
    event_fraction = maximum/abs(total["sum"]) if total["sum"] else None
    effective_sample_size = (
        total["sum"]**2/total["sum2"] if total["sum2"] > 0 else 0.0)
    stability = None
    if len(batch_estimates) >= 6:
        old_n = len(batch_estimates)-3
        old_mean = sum(x for x, _ in batch_estimates[:old_n])/old_n
        new_mean = sum(x for x, _ in batch_estimates[-3:])/3
        old_var = sum(v for _, v in batch_estimates[:old_n])/(old_n*old_n)
        new_var = sum(v for _, v in batch_estimates[-3:])/9
        stability = abs(new_mean-old_mean)/math.sqrt(old_var+new_var) \
                    if old_var+new_var else (0.0 if new_mean == old_mean else None)
    checks = {
        "relative_error": final.get("rse") is not None
            and final["rse"] <= float(cfg["target_rse"]),
        "maximum_event_fraction": event_fraction is not None
            and event_fraction <= float(settings["maximum_single_event_fraction"]),
        "stability": stability is not None
            and stability <= float(settings["stability_sigma"]),
        "minimum_batches": len(history) >= int(settings["minimum_batches"]),
        "nonzero_histories": nonzero >= int(settings["minimum_nonzero_histories"]),
        "variance_of_variance": final.get("vov") is not None
            and final["vov"] <= float(settings["maximum_variance_of_variance"]),
        "rse_scaling": rse_slope is not None
            and float(settings["rse_slope_minimum"]) <= rse_slope
            <= float(settings["rse_slope_maximum"]),
        "vov_scaling": vov_slope is not None
            and float(settings["vov_slope_minimum"]) <= vov_slope
            <= float(settings["vov_slope_maximum"]),
        "fom_stability": fom_range is not None
            and fom_range <= float(settings["maximum_fom_relative_range"]),
        "tail_index": tail_index is not None
            and tail_index >= float(settings["minimum_tail_index"]),
        "next_history": next_fraction is not None
            and next_fraction <= float(settings["maximum_next_history_fraction"]),
    }
    return {**final, "relative_error": final.get("rse"),
            "variance_of_variance": final.get("vov"),
            "batches": len(history), "nonzero_histories": nonzero,
            "effective_sample_size": effective_sample_size,
            "maximum_event_fraction": event_fraction, "stability_z": stability,
            "rse_slope": rse_slope, "vov_slope": vov_slope,
            "fom_relative_range": fom_range, "tail_index": tail_index,
            "next_history_fraction": next_fraction, "checks": checks,
            "passed": all(checks.values()), "checkpoints": checkpoints}


def batch_records(run_dir: Path) -> list[dict[str, Any]]:
    def records_from(directory: Path) -> list[dict[str, Any]]:
        found = []
        for path in sorted(directory.glob("batch_*.json")):
            try:
                value = json.loads(path.read_text())
                if value.get("complete"):
                    found.extend(value["blocks"] if "blocks" in value else [value])
            except (OSError, json.JSONDecodeError):
                pass
        return found

    records = records_from(run_dir)
    campaign_file = run_dir/"campaign.json"
    if campaign_file.is_file():
        try:
            metadata = json.loads(campaign_file.read_text())
            for legacy in metadata.get("legacy_runs", []):
                legacy_dir = Path(legacy)
                if not legacy_dir.is_absolute():
                    legacy_dir = REPO/legacy_dir
                records.extend(records_from(legacy_dir))
        except (OSError, json.JSONDecodeError):
            pass
    return records


def production_batch_count(run_dir: Path) -> int:
    count = 0
    for path in run_dir.glob("batch_*.json"):
        try:
            if json.loads(path.read_text()).get("complete"):
                count += 1
        except (OSError, json.JSONDecodeError):
            pass
    return count


def produced_event_count(run_dir: Path) -> int:
    """Count accepted and analysis-pending histories, never partial runs."""
    by_batch: dict[int, int] = {}
    for path in run_dir.glob("batch_*.json"):
        match = re.fullmatch(r"batch_(\d+)(?:\.produced)?\.json", path.name)
        if not match:
            continue
        try:
            record = json.loads(path.read_text())
            if record.get("complete") or path.name.endswith(".produced.json"):
                by_batch[int(match.group(1))] = int(record["events"])
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass
    return sum(by_batch.values())


def next_batch_number(run_dir: Path) -> int:
    # A failed remoll invocation can leave a macro, log, and partial ROOT file.
    # Those artifacts do not constitute a batch. Reuse the first number which
    # has neither an accepted analysis record nor a produced marker.
    number = 1
    while ((run_dir/f"batch_{number:05d}.json").is_file() or
           (run_dir/f"batch_{number:05d}.produced.json").is_file()):
        number += 1
    return number


def upgrade_region_statistics(run_dir: Path, septant: int) -> None:
    """Reanalyze accepted ROOT files made before ring-region tallies existed."""
    required = "group/main_detector_open"
    for path in sorted(run_dir.glob("batch_*.json")):
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not record.get("complete"):
            continue
        old_blocks = record.get("blocks", [record])
        if old_blocks and all(
                required in block.get("statistics", {}) for block in old_blocks):
            continue
        root_file = Path(record["root_file"])
        file_events = int(record.get("events", 1000))
        if not root_file.is_file() or file_events % 1000:
            raise RuntimeError(
                f"cannot upgrade regional statistics for {path.name}: "
                "the accepted ROOT file is missing or not divisible into "
                "1,000-history blocks")
        print(f"  upgrading {path.name} with exact ring-region statistics "
              "(no simulation rerun)", flush=True)
        blocks = []
        for start_entry in range(0, file_events, 1000):
            old = old_blocks[min(start_entry // 1000, len(old_blocks)-1)]
            block = analyze(root_file, septant,
                            float(old.get("elapsed_seconds", 0.0)),
                            start_entry, 1000)
            block.update({
                "batch": int(record["batch"]),
                "block": start_entry // 1000 + 1,
                "seed": record["seed"],
                "root_file": record["root_file"],
                "macro_file": record["macro_file"],
            })
            blocks.append(block)
        record["blocks"] = blocks
        atomic_json(path, record)


def snapshot(run_dir: Path, cfg: dict[str, Any], metadata: dict[str, Any]) -> dict[str, Any]:
    history = batch_records(run_dir)
    groups = {}
    if history:
        for key in ("group/closed", "group/transition", "group/open",
                    "group/showermax_total", "group/ring1", "group/ring2",
                    "group/ring3", "group/ring4", "group/ring5",
                    "group/ring6", "group/main_detector",
                    "group/ring1_closed", "group/ring1_transition",
                    "group/ring1_open", "group/ring2_closed",
                    "group/ring2_transition", "group/ring2_open",
                    "group/ring3_closed", "group/ring3_transition",
                    "group/ring3_open", "group/ring4_closed",
                    "group/ring4_transition", "group/ring4_open",
                    "group/ring5_closed", "group/ring5_transition",
                    "group/ring5_open", "group/ring6_closed",
                    "group/ring6_transition", "group/ring6_open",
                    "group/main_detector_closed",
                    "group/main_detector_transition",
                    "group/main_detector_open"):
            applicable = [
                batch for batch in history if key in batch.get("statistics", {})]
            if applicable:
                groups[key.split("/")[1]] = metrics(applicable, key, cfg)
    tile_statistics = {}
    tile_keys = sorted({
        key for batch in history for key in batch.get("statistics", {})
        if key.startswith("tile/")})
    for key in tile_keys:
        applicable = [
            batch for batch in history if key in batch.get("statistics", {})]
        tile_statistics[key.split("/", 1)[1]] = metrics(
            applicable, key, cfg)
    tiles = {}
    heatmap: dict[str, dict[str, float]] = {}
    species: dict[str, dict[str, float | int]] = {}
    energy_spectra: dict[str, dict[int, list[float | int]]] = {}
    global_maps: dict[str, dict[str, float]] = {}
    global_hits: dict[str, int] = {}
    tile_blocks = 0
    heatmap_blocks: dict[str, int] = {}
    species_blocks: dict[str, int] = {}
    spectrum_blocks: dict[str, int] = {}
    global_map_blocks: dict[str, int] = {}
    for batch in history:
        if batch.get("tiles"):
            tile_blocks += 1
        for name, value in batch["tiles"].items():
            target = tiles.setdefault(name, {"detector": value["detector"],
                                              "hits": 0, "response_rate": 0.0})
            target["hits"] += value["hits"]
            target["response_rate"] += value["response_rate"]
        for name, bins in batch.get("heatmap", {}).items():
            heatmap_blocks[name] = heatmap_blocks.get(name, 0) + 1
            target_bins = heatmap.setdefault(name, {})
            for x, y, weight in bins:
                key = f"{x},{y}"
                target_bins[key] = target_bins.get(key, 0.0) + weight
        for pid, value in batch.get("species", {}).items():
            species_blocks[pid] = species_blocks.get(pid, 0) + 1
            target = species.setdefault(pid, {
                "hits": 0, "ignored": 0, "rate": 0.0, "ignored_rate": 0.0})
            target["hits"] += value["hits"]
            target["ignored"] += value["ignored"]
            target["rate"] += value["rate"]
            target["ignored_rate"] += value["ignored_rate"]
        for pid, bins in batch.get("energy_spectra", {}).items():
            spectrum_blocks[pid] = spectrum_blocks.get(pid, 0) + 1
            target_bins = energy_spectra.setdefault(pid, {})
            for bin_index, hits, ignored, rate, ignored_rate in bins:
                target = target_bins.setdefault(
                    bin_index, [0, 0, 0.0, 0.0])
                target[0] += hits
                target[1] += ignored
                target[2] += rate
                target[3] += ignored_rate
        for name, hits in batch.get("global_hits", {}).items():
            global_hits[name] = global_hits.get(name, 0) + hits
        for name, bins in batch.get("global_maps", {}).items():
            global_map_blocks[name] = global_map_blocks.get(name, 0) + 1
            target_bins = global_maps.setdefault(name, {})
            for x, y, weight in bins:
                key = f"{x},{y}"
                target_bins[key] = target_bins.get(key, 0.0) + weight
    if history:
        for value in tiles.values():
            value["response_rate"] /= tile_blocks
        for name, bins in heatmap.items():
            for key in bins:
                bins[key] /= heatmap_blocks[name]
        for pid, value in species.items():
            value["rate"] /= species_blocks[pid]
            value["ignored_rate"] /= species_blocks[pid]
        for pid, bins in energy_spectra.items():
            for value in bins.values():
                value[2] /= spectrum_blocks[pid]
                value[3] /= spectrum_blocks[pid]
        for name, bins in global_maps.items():
            for key in bins:
                bins[key] /= global_map_blocks[name]
    result = {"title": cfg["title"], "metadata": metadata,
              "batches": len(history),
              "production_batches": production_batch_count(run_dir),
              "events": sum(item["events"] for item in history),
              "ignored_response_hits": sum(
                  item.get("unsupported_hits", 0) for item in history),
              "elapsed_seconds": sum(item["elapsed_seconds"] for item in history),
              "groups": groups, "tiles": tiles,
              "tile_statistics": tile_statistics,
              "heatmap": {name: [[*map(int, key.split(",")), weight]
                                 for key, weight in bins.items()]
                          for name, bins in heatmap.items()},
              "species": species,
              "energy_spectra": {
                  pid: [[index, *value] for index, value in sorted(bins.items())]
                  for pid, bins in energy_spectra.items()},
              "global_hits": global_hits,
              "global_maps": {
                  name: [[*map(int, key.split(",")), weight]
                         for key, weight in bins.items()]
                  for name, bins in global_maps.items()},
              "targets": {"rse": cfg["target_rse"], **cfg["reliability"]},
              "state": (run_dir/"control").read_text().strip()
                       if (run_dir/"control").is_file() else "ready"}
    atomic_json(run_dir/"snapshot.json", result)
    return result


class DashboardHandler(BaseHTTPRequestHandler):
    run_dir: Path
    cfg: dict[str, Any]

    def selected_run_dir(self) -> Path:
        query = parse_qs(urlparse(self.path).query)
        name = query.get("run", [self.run_dir.name])[0]
        # Campaigns are direct children of one fixed output directory. Requiring
        # a bare directory name prevents the viewer endpoint from becoming an
        # arbitrary filesystem reader/writer.
        if Path(name).name != name:
            raise ValueError("invalid campaign directory")
        selected = self.run_dir.parent/name
        if not (selected/"campaign.json").is_file():
            raise ValueError("unknown campaign directory")
        return selected

    def runs_payload(self) -> bytes:
        runs = []
        for campaign_file in self.run_dir.parent.glob("*/campaign.json"):
            directory = campaign_file.parent
            try:
                metadata = json.loads(campaign_file.read_text())
                saved = json.loads((directory/"snapshot.json").read_text()) \
                    if (directory/"snapshot.json").is_file() else {}
            except (OSError, json.JSONDecodeError):
                continue
            active = False
            marker = directory/"active.json"
            if marker.is_file():
                try:
                    pid = int(json.loads(marker.read_text())["pid"])
                    os.kill(pid, 0)
                    active = True
                except (OSError, ValueError, KeyError, json.JSONDecodeError):
                    active = False
            control = (directory/"control").read_text().strip() \
                if (directory/"control").is_file() else "ready"
            state = "live" if active else (
                control if control in {"pause", "stop"} else "saved")
            runs.append({
                "name": directory.name,
                "cell": metadata.get("cell", directory.name),
                "sieve": metadata.get("sieve", "?"),
                "events": int(saved.get("events", 0)),
                "state": state,
                "diagnostic": directory.name.startswith("diagnostic"),
                "modified": campaign_file.stat().st_mtime,
            })
        runs.sort(key=lambda item: (
            item["diagnostic"], item["cell"], item["sieve"], item["name"]))
        return json.dumps(
            {"selected": self.run_dir.name, "runs": runs}).encode()

    def do_GET(self) -> None:
        route = urlparse(self.path).path
        if route in ("/", "/index.html"):
            data, kind = DASHBOARD.read_bytes(), "text/html; charset=utf-8"
        elif route == "/api/runs":
            data, kind = self.runs_payload(), "application/json"
        elif route == "/api/state":
            try:
                selected = self.selected_run_dir()
                snapshot_file = selected/"snapshot.json"
                if not snapshot_file.is_file():
                    metadata = json.loads((selected/"campaign.json").read_text())
                    snapshot(selected, self.cfg, metadata)
                data, kind = snapshot_file.read_bytes(), "application/json"
            except (OSError, ValueError, json.JSONDecodeError) as error:
                self.send_error(404, str(error)); return
        else:
            self.send_error(404); return
        self.send_response(200)
        self.send_header("Content-Type", kind)
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)
    def do_POST(self) -> None:
        action = urlparse(self.path).path.removeprefix("/api/")
        if action not in {"pause", "resume", "stop"}:
            self.send_error(404); return
        try:
            selected = self.selected_run_dir()
        except ValueError as error:
            self.send_error(404, str(error)); return
        (selected/"control").write_text(
            ("running" if action == "resume" else action) + "\n")
        self.send_response(204); self.end_headers()
    def log_message(self, *_: Any) -> None:
        pass


def serve(run_dir: Path, port: int, cfg: dict[str, Any]) -> ThreadingHTTPServer:
    handler = type("RunDashboardHandler", (DashboardHandler,), {
        "run_dir": run_dir, "cfg": cfg})
    server = ThreadingHTTPServer(("127.0.0.1", port), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server


def run_campaign(args: argparse.Namespace, cfg: dict[str, Any],
                 row: dict[str, str], run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"schema": "full_support_tiles_v5",
                "cell": args.cell, "septant": args.septant, "sieve": args.sieve,
                "energy_mev": row["energy_mev"], "target": row["target"],
                "channel": row["channel"], "batch_events": cfg["batch_events"]}
    existing_meta = run_dir/"campaign.json"
    if existing_meta.is_file():
        previous = json.loads(existing_meta.read_text())
        identity = ("schema", "cell", "septant", "sieve", "energy_mev",
                    "target", "channel")
        if any(previous.get(key) != metadata.get(key) for key in identity):
            raise RuntimeError(
                f"{run_dir} belongs to a different campaign selection")
        sizes = set(previous.get("production_batch_sizes",
                                 [previous.get("batch_events", 1000)]))
        sizes.add(int(metadata["batch_events"]))
        metadata["production_batch_sizes"] = sorted(sizes)
    else:
        metadata["production_batch_sizes"] = [int(metadata["batch_events"])]
    atomic_json(existing_meta, metadata)
    (run_dir/"control").write_text("running\n")
    atomic_json(run_dir/"active.json", {
        "pid": os.getpid(), "started": time.time()})
    upgrade_region_statistics(run_dir, args.septant)
    snapshot(run_dir, cfg, metadata)
    server = None
    if not args.no_dashboard:
        server = serve(run_dir, args.port, cfg)
        print(f"Dashboard: http://127.0.0.1:{args.port}")
        if not args.no_browser:
            webbrowser.open(f"http://127.0.0.1:{args.port}")
    # surface-replay is the build that contains the current kinematic-bias
    # messengers. build-develop may exist but predate those commands.
    remoll = next((path for path in (REPO/"build-surface-replay/remoll",
                   REPO/"build-develop/remoll", REPO/"bin/remoll")
                   if path.is_file()), None)
    if not remoll:
        raise RuntimeError("remoll executable was not found")
    analysis_queue: queue.Queue[dict[str, Any] | None] = queue.Queue()

    def analysis_worker() -> None:
        while True:
            pending = analysis_queue.get()
            if pending is None:
                return
            number = int(pending["batch"])
            try:
                root_file = Path(pending["root_file"])
                file_events = int(pending["events"])
                if file_events % 1000:
                    raise RuntimeError(
                        f"batch contains {file_events} events; statistical "
                        "analysis requires an exact multiple of 1000")
                blocks = []
                for start_entry in range(0, file_events, 1000):
                    block = analyze(
                        root_file, args.septant,
                        float(pending["elapsed_seconds"]) * 1000 / file_events,
                        start_entry, 1000)
                    block.update({
                        "batch": number, "block": start_entry // 1000 + 1,
                        "seed": pending["seed"], "root_file": pending["root_file"],
                        "macro_file": pending["macro_file"],
                    })
                    blocks.append(block)
                atomic_json(run_dir/f"batch_{number:05d}.json", {
                    "complete": True, "batch": number, "seed": pending["seed"],
                    "events": file_events, "root_file": pending["root_file"],
                    "macro_file": pending["macro_file"], "blocks": blocks,
                })
                Path(pending["log_file"]).unlink(missing_ok=True)
                (run_dir/f"batch_{number:05d}.produced.json").unlink(missing_ok=True)
                (run_dir/f"batch_{number:05d}.analysis_error.json").unlink(
                    missing_ok=True)
                view = snapshot(run_dir, cfg, metadata)
                closed = sum(group["passed"] for group in view["groups"].values())
                print(f"  analysis accepted batch {number}: {view['events']} "
                      f"analyzed events; {closed}/{len(view['groups'])} "
                      "detector groups converged", flush=True)
            except Exception as error:
                atomic_json(run_dir/f"batch_{number:05d}.analysis_error.json",
                            {"batch": number, "error": str(error)})
                print(f"ANALYSIS ERROR batch {number}: {error}", file=sys.stderr,
                      flush=True)
            finally:
                analysis_queue.task_done()

    threading.Thread(target=analysis_worker, daemon=True).start()
    # Resume any completely produced batch that was awaiting display analysis.
    for marker in sorted(run_dir.glob("batch_*.produced.json")):
        try:
            pending = json.loads(marker.read_text())
            if not (run_dir/f"batch_{int(pending['batch']):05d}.json").is_file():
                analysis_queue.put(pending)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass
    interrupted = False
    def handle_interrupt(*_: Any) -> None:
        nonlocal interrupted
        interrupted = True
        (run_dir/"control").write_text("pause\n")
    old = signal.signal(signal.SIGINT, handle_interrupt)
    try:
        while True:
            control = (run_dir/"control").read_text().strip()
            if interrupted or control == "stop":
                print(f"Campaign {'paused' if interrupted else 'stopped'} "
                      "at a safe batch boundary.")
                break
            if control == "pause":
                time.sleep(1.0)
                continue
            produced = production_batch_count(run_dir) + len(
                list(run_dir.glob("batch_*.produced.json")))
            if args.max_batches and produced >= args.max_batches:
                print("Requested batch limit reached.")
                break
            completed_events = produced_event_count(run_dir)
            if args.max_events and completed_events >= args.max_events:
                print(f"Requested event limit reached ({completed_events} "
                      "histories).")
                break
            batch_events = int(cfg["batch_events"])
            if args.max_events:
                batch_events = min(batch_events,
                                   args.max_events - completed_events)
            batch_cfg = {**cfg, "batch_events": batch_events}
            number = next_batch_number(run_dir)
            root_file = run_dir/f"batch_{number:05d}.root"
            macro_file = run_dir/f"batch_{number:05d}.mac"
            log_file = run_dir/f"batch_{number:05d}.log"
            seed = args.seed + number
            macro_file.write_text(
                macro_text(batch_cfg, row, args.sieve, root_file, seed))
            print(f"Batch {number}: {batch_events} events", flush=True)
            start = time.monotonic()
            environment = os.environ.copy()
            if remoll.parent.name.startswith("build-"):
                environment["DYLD_LIBRARY_PATH"] = str(remoll.parent) + (
                    ":" + environment["DYLD_LIBRARY_PATH"]
                    if environment.get("DYLD_LIBRARY_PATH") else "")
            maximum_attempts = 3
            for attempt in range(1, maximum_attempts + 1):
                # Only an accepted JSON record makes a batch durable. Remove
                # any partial ROOT artifact left by a failed invocation so it
                # cannot be mistaken for the current attempt's output.
                root_file.unlink(missing_ok=True)
                with log_file.open("w") as log:
                    process = subprocess.run(
                        [str(remoll), "-t", str(cfg["threads"]),
                         "-m", str(macro_file)], cwd=REPO, stdout=log,
                        stderr=subprocess.STDOUT, env=environment)
                log_text = log_file.read_text(errors="replace")
                command_error = (
                    "COMMAND NOT FOUND" in log_text or
                    "Batch is interrupted" in log_text or
                    "ERROR: invalid /remoll/bias/" in log_text or
                    "ERROR: unknown /remoll/bias/" in log_text)
                produced_root = root_file.is_file() and root_file.stat().st_size > 0
                if command_error:
                    raise RuntimeError(
                        f"remoll rejected a generated macro command; inspect "
                        f"{log_file}. This is a campaign configuration error "
                        "and will not be retried")
                if (process.returncode == 0 and produced_root) or interrupted:
                    break
                print(f"  remoll batch {number} attempt {attempt} failed "
                      f"(return code {process.returncode}); "
                      f"{'retrying the same batch and seed' if attempt < maximum_attempts else 'retry limit reached'}",
                      file=sys.stderr, flush=True)
            elapsed = time.monotonic() - start
            if process.returncode:
                if interrupted:
                    print(f"Batch {number} interrupted; its incomplete ROOT file "
                          "will not be analyzed or counted.", flush=True)
                    break
                raise RuntimeError(
                    f"remoll failed after {maximum_attempts} attempts; "
                    f"inspect {log_file}. Accepted earlier batches are intact; "
                    f"restart will retry batch {number}")
            if not root_file.is_file() or root_file.stat().st_size == 0:
                tail = log_file.read_text(errors="replace")[-2000:]
                raise RuntimeError(
                    f"remoll produced no ROOT output after {maximum_attempts} "
                    "attempts; "
                    f"inspect {log_file}\n{tail}")
            pending = {
                "batch": number, "seed": seed, "elapsed_seconds": elapsed,
                "events": batch_events,
                "root_file": str(root_file), "macro_file": str(macro_file),
                "log_file": str(log_file),
            }
            atomic_json(run_dir/f"batch_{number:05d}.produced.json", pending)
            analysis_queue.put(pending)
            print(f"  produced batch {number} in {elapsed:.1f}s; "
                  "starting the next batch without waiting for display analysis",
                  flush=True)
    finally:
        signal.signal(signal.SIGINT, old)
        # This wait occurs only while shutting down the producer. It never
        # delays launching the next remoll batch during normal running.
        analysis_queue.join()
        snapshot(run_dir, cfg, metadata)
        (run_dir/"active.json").unlink(missing_ok=True)
        if server is not None:
            server.shutdown()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=HERE/"campaign.toml")
    parser.add_argument("--cell", default="c12_us:2200:c12_elastic")
    parser.add_argument("--septant", type=int, default=0)
    parser.add_argument("--sieve", choices=("auto", "none", "in", "out"),
                        default="auto")
    parser.add_argument("--run-name", default="")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--seed", type=int, default=710000)
    parser.add_argument("--max-batches", type=int, default=0)
    parser.add_argument("--max-events", type=int, default=0,
                        help="maximum accepted or analysis-pending histories "
                             "in this campaign; must be a multiple of 1000")
    parser.add_argument("--batch-events", type=int, default=0,
                        help="remoll events per ROOT file; must be a multiple of 1000")
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--no-dashboard", action="store_true")
    parser.add_argument("command", choices=("check", "start", "status", "pause",
                                            "resume", "stop", "serve", "display"))
    args = parser.parse_args()
    if not 0 <= args.septant <= 6:
        parser.error("--septant must be 0 through 6")
    cfg = load_config(args.config)
    if args.batch_events:
        cfg["batch_events"] = args.batch_events
    if int(cfg["batch_events"]) <= 0 or int(cfg["batch_events"]) % 1000:
        parser.error("--batch-events must be a positive multiple of 1000")
    if args.max_events < 0 or args.max_events % 1000:
        parser.error("--max-events must be zero or a positive multiple of 1000")
    if args.command in {"serve", "display"}:
        output_root = resolve(cfg["output_directory"])
        if args.run_name:
            run_dir = output_root/args.run_name
            campaign_file = run_dir/"campaign.json"
        else:
            candidates = sorted(
                (path for path in output_root.glob("*/campaign.json")
                 if not path.parent.name.startswith("diagnostic")),
                key=lambda path: path.stat().st_mtime, reverse=True)
            if not candidates:
                raise RuntimeError(f"no saved campaigns exist under {output_root}")
            campaign_file = candidates[0]
            run_dir = campaign_file.parent
        if not campaign_file.is_file():
            raise RuntimeError(f"saved campaign metadata is missing: {campaign_file}")
        metadata = json.loads(campaign_file.read_text())
        snapshot(run_dir, cfg, metadata)
        server = serve(run_dir, args.port, cfg)
        print(f"Displaying saved campaign: {run_dir.name}")
        print(f"Dashboard: http://127.0.0.1:{args.port} (Ctrl-C to stop)")
        if not args.no_browser:
            webbrowser.open(f"http://127.0.0.1:{args.port}")
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        return 0
    row = load_window(resolve(cfg["kinematics_table"]), args.cell, args.septant)
    if args.sieve == "auto":
        args.sieve = "none" if row["target"] == "lh2" else "out"
    name = args.run_name or (
        f"{args.cell.replace(':','_')}_{args.sieve}_fullsupport_v5")
    run_dir = resolve(cfg["output_directory"])/name
    metadata = {"schema": "full_support_tiles_v5",
                "cell": args.cell, "septant": args.septant, "sieve": args.sieve,
                "energy_mev": row["energy_mev"], "target": row["target"],
                "channel": row["channel"], "batch_events": cfg["batch_events"]}
    if args.command == "check":
        # Generate the complete macro in memory. This exercises channel
        # selection and validates every compact-table proposal bound before a
        # production process is launched.
        checked_macro = macro_text(
            cfg, row, args.sieve, Path("/tmp/showermax_live_check.root"),
            args.seed)
        if "/remoll/bias/th/" in checked_macro:
            raise RuntimeError("internal validation found obsolete theta command")
        for path in (resolve(cfg["base_macro"]), RESPONSE_DATA, ANALYZER, DASHBOARD):
            if not path.exists(): raise RuntimeError(f"missing required path: {path}")
        print(f"OK: {args.cell}; all 28 ShowerMax planes, all six rings, "
              "and the main detector are active")
        print("Azimuth is sampled physically over 0-360 degrees; no phi or "
              "septant bias is applied.")
        print("The remaining target-kinematic proposal retains a 20% physical "
              "component on each biased non-azimuthal axis.")
        return 0
    if args.command == "start":
        run_campaign(args, cfg, row, run_dir); return 0
    if args.command in {"pause", "resume", "stop"}:
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir/"control").write_text(
            ("running" if args.command == "resume" else args.command) + "\n")
        print(f"{args.command}: {run_dir}"); return 0
    if args.command == "status":
        view = snapshot(run_dir, cfg, metadata)
        print(json.dumps(view, indent=2)); return 0
    raise RuntimeError(f"unhandled command: {args.command}")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
