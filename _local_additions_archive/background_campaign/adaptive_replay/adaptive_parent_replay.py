#!/usr/bin/env python3
"""Restartable adaptive ShowerMax parent-surface replay experiment.

This remains an experimental extension of ``campaign.py``.  For a secondary
detector response, replay candidates must be the secondary's parent (or an
older ancestor) at a real upstream surface; replaying the secondary's own
surface crossing cannot improve the probability that created it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
CAMPAIGN_DIR = HERE.parent
REPO = next(path for path in (HERE, *HERE.parents)
            if (path / "geometry/mollerMother.gdml").is_file())
sys.path.insert(0, str(CAMPAIGN_DIR))
from campaign import (  # noqa: E402
    batch_variance, merge_largest, new_merged_stat, reliability_result,
    statistic_result,
)

BUILD_MACRO = HERE / "build_surface_source.C"
PILOT_BUILD_MACRO = HERE / "build_surface_pilot.C"
PILOT_ANALYZE_MACRO = HERE / "analyze_surface_pilot.C"
ANALYZE_MACRO = HERE / "analyze_surface_scores.C"
DEFAULT_INPUT = REPO / (
    "_local_additions_archive/calibration/neff_phase_space/"
    "neff_carbon_US_2p2_100k_sieveout_sieve_out.root"
)
DEFAULT_OUTPUT = CAMPAIGN_DIR / "runs" / "adaptive_parent_replay"
STATE_VERSION = 2
TAIL_SIZE = 201


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def root_string(value: str | Path) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def run_macro(path: Path, arguments: list[str], require: str) -> str:
    invocation = f"{path}({','.join(arguments)})"
    result = subprocess.run(
        [str(REPO / "build/reroot"), "-l", "-b", "-q", invocation],
        cwd=REPO, text=True, capture_output=True,
    )
    output = result.stdout + result.stderr
    if result.returncode != 0 or require not in output or "_ERROR" in output:
        raise RuntimeError(f"ROOT helper failed:\n{output[-5000:]}")
    return output


def analyze_file(path: Path, mode: str, score_scale: int, output: Path) -> dict[str, Any]:
    run_macro(ANALYZE_MACRO, [root_string(path), root_string(mode),
                              str(score_scale), root_string(output)], "SCORE_ANALYSIS")
    rows: dict[str, Any] = {}
    with output.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            rows[row["key"]] = {
                "events": int(row["events"]),
                "signal_events": int(row["signal"]),
                "sum": float(row["sum"]), "sum2": float(row["sum2"]),
                "sum3": float(row["sum3"]), "sum4": float(row["sum4"]),
                "maximum": float(row["maximum"]),
                "largest": ([] if row["largest"] == "-" else
                            [float(item) for item in row["largest"].split(",")]),
            }
    if len(rows) != 62:
        raise RuntimeError(f"expected 62 ShowerMax statistics, found {len(rows)}")
    return rows


def direct_result(value: dict[str, Any]) -> dict[str, Any]:
    events = value["events"]
    mean = value["sum"] / events
    variance = batch_variance(events, value)
    standard_error = math.sqrt(variance / events) if events else 0.0
    return {
        "estimate": mean,
        "relative_error": standard_error / abs(mean) if mean else None,
        "effective_sample_size": (
            value["sum"] ** 2 / value["sum2"] if value["sum2"] > 0.0 else 0.0
        ),
        "signal_events": value["signal_events"],
        "largest_event_fraction": (
            value["maximum"] / abs(value["sum"]) if value["sum"] else None
        ),
    }


def settings(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "precision": {
            "relative_target": args.target_rse,
            "maximum_single_event_fraction": args.maximum_single_fraction,
            "stability_sigma": args.stability_sigma,
        },
        "reliability": {
            "minimum_ordinary_batches": args.minimum_batches,
            "minimum_nonzero_histories": args.minimum_nonzero,
            "maximum_variance_of_variance": args.maximum_vov,
            "rse_slope_minimum": args.rse_slope_minimum,
            "rse_slope_maximum": args.rse_slope_maximum,
            "vov_slope_minimum": args.vov_slope_minimum,
            "vov_slope_maximum": args.vov_slope_maximum,
            "maximum_fom_relative_range": args.maximum_fom_range,
            "minimum_tail_index": args.minimum_tail_index,
            "maximum_next_history_fraction": args.maximum_next_fraction,
        },
    }


def config_record(args: argparse.Namespace, input_hash: str) -> dict[str, Any]:
    return {
        "input": str(args.input.resolve()), "input_sha256": input_hash,
        "batch_events": args.batch_events, "threads": args.threads,
        "surface_detector": args.surface, "physical_fraction": args.physical_fraction,
        "control_interval": args.control_interval,
        "top_histories": args.top_histories, "physics_list": args.physics_list,
        "pilot_events": args.pilot_events, "pilot_unique_states": args.pilot_unique_states,
        "minimum_pilot_successes": args.minimum_pilot_successes,
        "maximum_pilots_per_target": args.maximum_pilots_per_target,
        "settings": settings(args),
    }


def merge_stat(target: dict[str, Any], value: dict[str, Any], batch_name: str,
               elapsed: float, proposal: dict[str, Any]) -> None:
    events = value["events"]
    target["signal_events"] += value["signal_events"]
    for field in ("sum", "sum2", "sum3", "sum4"):
        target[field] += value[field]
    target["maximum"] = max(target["maximum"], value["maximum"])
    target["largest"] = merge_largest(target["largest"], value["largest"])
    variance_numerator = events * batch_variance(events, value)
    target["variance_numerator"] += variance_numerator
    target["maximum_batch_sum"] = max(target["maximum_batch_sum"], abs(value["sum"]))
    target["history"].append({
        "batch": batch_name, "events": events, "sum": value["sum"],
        "sum2": value["sum2"], "sum3": value["sum3"], "sum4": value["sum4"],
        "variance_numerator": variance_numerator, "elapsed_seconds": elapsed,
        "proposal": proposal,
    })


def original_event_count(original: dict[str, Any]) -> int:
    counts = {int(value["events"]) for value in original.values()}
    if len(counts) != 1 or next(iter(counts)) <= 0:
        raise RuntimeError(f"original statistics have inconsistent event counts: {counts}")
    return next(iter(counts))


def original_aggregate(original: dict[str, Any]) -> dict[str, Any]:
    aggregate = {key: new_merged_stat() for key in original}
    proposal = {"kind": "original", "source": "unbiased target run"}
    for key, value in original.items():
        merge_stat(aggregate[key], value, "original", 0.0, proposal)
    return aggregate


def merge_existing_aggregate(target: dict[str, Any], source: dict[str, Any]) -> None:
    """Merge persisted raw moments while retaining their per-batch history."""
    target["signal_events"] += source["signal_events"]
    for field in ("sum", "sum2", "sum3", "sum4", "variance_numerator"):
        target[field] += source[field]
    target["maximum"] = max(target["maximum"], source["maximum"])
    target["maximum_batch_sum"] = max(
        target["maximum_batch_sum"], source["maximum_batch_sum"]
    )
    target["largest"] = merge_largest(target["largest"], source["largest"])
    target["history"].extend(source["history"])


def detector_group(detector: int) -> str | None:
    if detector < 70030 or detector > 72730 or (detector - 70030) % 100:
        return None
    stack = (detector - 70030) // 100
    if stack % 4 == 0:
        return "showermax_open"
    if stack % 4 == 2:
        return "showermax_closed"
    return "showermax_transition"


def history_targets_statistic(key: str, item: dict[str, Any]) -> bool:
    """Use only proposals deliberately optimized for this statistic."""
    if item.get("batch") == "original":
        return False
    proposal = item.get("proposal", {})
    if proposal.get("kind") != "targeted":
        return False
    observable = key.split("/", 1)[1]
    if proposal.get("observable") != observable:
        return False
    detector = int(proposal.get("target_detector", -1))
    if key.startswith("tile_"):
        return detector == int(key.split("/", 1)[0].removeprefix("tile_"))
    return detector_group(detector) == key.split("/", 1)[0]


def stratum_estimate(items: list[dict[str, Any]]) -> tuple[float, float] | None:
    """Return (mean, variance of mean) for one fixed proposal stratum."""
    events = sum(int(item["events"]) for item in items)
    if events <= 1:
        return None
    total = sum(float(item["sum"]) for item in items)
    total2 = sum(float(item["sum2"]) for item in items)
    variance = max(0.0, (total2 - total * total / events) / (events - 1))
    variance_of_mean = variance / events
    if variance_of_mean <= 0.0:
        return None
    return total / events, variance_of_mean


def precision_combination(key: str, aggregate: dict[str, Any]) -> dict[str, float] | None:
    """Combine the original estimator with useful targeted replay as strata.

    Incidental scores from proposals targeting other tiles remain available in
    the ROOT file, but they are not allowed to destabilize this tally.  The
    original and matched targeted estimators are each unbiased; inverse-variance
    stratification preserves the information from the much more efficient one.
    """
    original = [item for item in aggregate["history"] if item["batch"] == "original"]
    targeted = [item for item in aggregate["history"]
                if history_targets_statistic(key, item)]
    strata = [result for result in
              (stratum_estimate(original), stratum_estimate(targeted))
              if result is not None]
    if not strata:
        return None
    precision = sum(1.0 / variance for _, variance in strata)
    estimate = sum(mean / variance for mean, variance in strata) / precision
    standard_error = math.sqrt(1.0 / precision)
    relative_error = standard_error / abs(estimate) if estimate else None
    return {
        "estimate": estimate,
        "standard_error": standard_error,
        "relative_error": relative_error,
        "effective_sample_size": (
            1.0 / relative_error ** 2 if relative_error not in (None, 0.0) else 0.0
        ),
        "targeted_histories": sum(int(item["events"]) for item in targeted),
    }


def campaign_results(state: dict[str, Any]) -> dict[str, Any]:
    """Return moments for the original sample plus every accepted replay draw."""
    total_events = (state["original_events"]
                    + len(state["batches"]) * state["config"]["batch_events"])
    campaign_settings = state["config"]["settings"]
    results = {}
    for key, aggregate in state["aggregate"].items():
        value = dict(aggregate)
        value.update(statistic_result(total_events, value))
        value["effective_sample_size"] = (
            value["sum"] ** 2 / value["sum2"] if value["sum2"] > 0.0 else 0.0
        )
        combined = precision_combination(key, aggregate)
        if combined is not None:
            value.update(combined)
        control_value = dict(state["control_aggregate"][key])
        value["reliability"] = reliability_result(campaign_settings, control_value)
        results[key] = value
    return results


def gate(value: dict[str, Any], config: dict[str, Any]) -> tuple[bool, list[str]]:
    precision = config["settings"]["precision"]
    reliability = config["settings"]["reliability"]
    failed = []
    if value.get("estimate", 0.0) == 0.0:
        failed.append("zero_response")
    if value.get("signal_events", 0) < reliability["minimum_nonzero_histories"]:
        failed.append("nonzero_histories")
    if value.get("relative_error") is None or value["relative_error"] > precision["relative_target"]:
        failed.append("relative_error")
    fraction = value.get("largest_event_fraction")
    if fraction is None or fraction > precision["maximum_single_event_fraction"]:
        failed.append("single_history")
    stability = value.get("stability_z")
    if stability is None or stability > precision["stability_sigma"]:
        failed.append("stability")
    failed.extend(value.get("reliability", {}).get("failed_checks", []))
    failed = list(dict.fromkeys(failed))
    return not failed, failed


def required_keys() -> list[str]:
    keys = []
    for detector in range(70030, 72731, 100):
        keys.extend((f"tile_{detector}/rate", f"tile_{detector}/rate_energy"))
    for group in ("showermax_open", "showermax_closed", "showermax_transition"):
        keys.extend((f"{group}/rate", f"{group}/rate_energy"))
    return keys


def priority(value: dict[str, Any], config: dict[str, Any]) -> float:
    precision = config["settings"]["precision"]
    reliability = config["settings"]["reliability"]
    scores = [
        (value.get("relative_error") or 10.0) / precision["relative_target"],
        (value.get("largest_event_fraction") or 1.0)
        / precision["maximum_single_event_fraction"],
        reliability["minimum_nonzero_histories"] / max(1, value.get("signal_events", 0)),
    ]
    diagnostic = value.get("reliability", {})
    if diagnostic.get("variance_of_variance") is not None:
        scores.append(diagnostic["variance_of_variance"] /
                      reliability["maximum_variance_of_variance"])
    if diagnostic.get("next_largest_history_fraction") is not None:
        scores.append(diagnostic["next_largest_history_fraction"] /
                      reliability["maximum_next_history_fraction"])
    return max(scores)


def select_target(state: dict[str, Any], original: dict[str, Any]) -> tuple[int, str, str]:
    results = campaign_results(state)
    precision_candidates = []
    diagnostic_candidates = []
    for detector in range(70030, 72731, 100):
        for observable in ("rate", "rate_energy"):
            key = f"tile_{detector}/{observable}"
            passed, _ = gate(results[key], state["config"])
            if passed:
                continue
            rse = results[key].get("relative_error")
            rse_score = ((rse or 10.0) /
                         state["config"]["settings"]["precision"]["relative_target"])
            diagnostic_candidates.append((priority(results[key], state["config"]),
                                          detector, observable, key))
            if rse is None or rse > state["config"]["settings"]["precision"]["relative_target"]:
                precision_candidates.append((rse_score, detector, observable, key))
    candidates = precision_candidates or diagnostic_candidates
    if not candidates:
        return 70630, "rate_energy", "all tile gates passed"
    score, detector, observable, key = max(candidates)
    reason = (f"worst unresolved tile RSE {key}, {100*(results[key].get('relative_error') or 0):.2f}%"
              if precision_candidates else
              f"worst remaining diagnostic {key}, priority {score:.3g}")
    return detector, observable, reason


def target_key(detector: int, observable: str) -> str:
    return f"tile_{detector}/{observable}"


def initialize_target_lock(args: argparse.Namespace, state: dict[str, Any],
                           original: dict[str, Any]) -> None:
    """Migrate an existing run and establish one target that cannot be skipped."""
    if not state.get("target_lock_initialized", False):
        # Older controller versions moved on after precision alone.  Resume the
        # most recently targeted unfinished tally instead of silently selecting
        # yet another tile.
        results = campaign_results(state)
        for batch in reversed(state["batches"]):
            proposal = batch["proposal"]
            if proposal.get("kind") != "targeted":
                continue
            detector = int(proposal["target_detector"])
            observable = str(proposal["observable"])
            key = target_key(detector, observable)
            if not gate(results[key], state["config"])[0]:
                state["locked_target"] = {
                    "detector": detector, "observable": observable,
                    "since_batch": batch["name"],
                }
                break
        state["target_lock_initialized"] = True
        save_state(args, state)

    locked = state.get("locked_target")
    if locked:
        results = campaign_results(state)
        key = target_key(int(locked["detector"]), str(locked["observable"]))
        passed, _ = gate(results[key], state["config"])
        if passed:
            print(f"Full convergence passed for locked target {key}; releasing it.")
            state["locked_target"] = None
            save_state(args, state)


def locked_selection(args: argparse.Namespace, state: dict[str, Any],
                     original: dict[str, Any]) -> tuple[int, str, str]:
    initialize_target_lock(args, state, original)
    locked = state.get("locked_target")
    if not locked:
        detector, observable, reason = select_target(state, original)
        state["locked_target"] = {
            "detector": detector, "observable": observable,
            "since_batch": f"before_batch_{len(state['batches']) + 1:04d}",
        }
        save_state(args, state)
        return detector, observable, reason
    detector = int(locked["detector"])
    observable = str(locked["observable"])
    key = target_key(detector, observable)
    _, failed = gate(campaign_results(state)[key], state["config"])
    return detector, observable, (
        f"locked on {key} until full convergence; pending {','.join(failed)}"
    )


def generated_macro(config: dict[str, Any], source: Path, result: Path,
                    seed: int, events: int | None = None) -> str:
    event_count = config["batch_events"] if events is None else events
    lines = [
        "/remoll/setgeofile geometry/mollerMother.gdml",
        f"/remoll/physlist/register {config['physics_list']}",
        "/remoll/physlist/parallel/enable",
        "/remoll/parallel/setfile geometry/mollerParallel.gdml",
        "/run/initialize", "/control/execute macros/load_magnetic_fieldmaps.mac",
        "/remoll/evgen/set external", "/remoll/evgen/copyRate 1",
        f"/remoll/evgen/external/file {source}",
        "/remoll/evgen/external/detid 990035",
        "/remoll/evgen/external/startEvent 0", "/remoll/SD/disable_all",
    ]
    for detector in range(70030, 72731, 100):
        lines.extend((f"/remoll/SD/enable {detector}",
                      f"/remoll/SD/detect surfacehits {detector}",
                      f"/remoll/SD/detect boundaryhits {detector}",
                      f"/remoll/SD/detect secondaries {detector}",
                      f"/remoll/SD/detect lowenergyneutral {detector}"))
    lines.extend((f"/remoll/seed {seed}", f"/remoll/filename {result}",
                  f"/run/beamOn {event_count}", ""))
    return "\n".join(lines)


def learning_tree_name(detector: int, observable: str) -> str:
    return f"learned_{detector}_{observable}"


def root_entries(path: Path) -> int:
    result = subprocess.run(
        [str(REPO / "build/reroot"), "-l", "-b", "-q", str(path), "-e",
         'printf("SURFACE_ENTRIES %lld\\n",T ? T->GetEntries() : -1);'],
        cwd=REPO, text=True, capture_output=True,
    )
    match = re.search(r"SURFACE_ENTRIES (-?\d+)", result.stdout + result.stderr)
    return int(match.group(1)) if match else -1


def named_tree_entries(path: Path, name: str) -> int:
    expression = (
        f'auto* tree_named=(TTree*)gFile->Get("{name}");'
        'printf("SURFACE_NAMED_ENTRIES %lld\\n",tree_named ? tree_named->GetEntries() : -1);'
    )
    result = subprocess.run(
        [str(REPO / "build/reroot"), "-l", "-b", "-q", str(path), "-e", expression],
        cwd=REPO, text=True, capture_output=True,
    )
    match = re.search(r"SURFACE_NAMED_ENTRIES (-?\d+)", result.stdout + result.stderr)
    return int(match.group(1)) if match else -1


def merge_root(accumulated: Path, batch: Path, expected_previous: int) -> None:
    observed = root_entries(accumulated) if accumulated.is_file() else 0
    batch_entries = root_entries(batch)
    if batch_entries <= 0:
        raise RuntimeError(f"invalid replay batch ROOT file: {batch}")
    if observed != expected_previous:
        if observed == expected_previous + batch_entries:
            return
        raise RuntimeError(f"cumulative ROOT has {observed} entries; expected {expected_previous}")
    temporary = accumulated.with_suffix(".merge.root")
    if observed == 0:
        shutil.copy2(batch, temporary)
    else:
        result = subprocess.run(
            ["hadd", "-f", str(temporary), str(accumulated), str(batch)],
            cwd=REPO, text=True, capture_output=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"hadd failed:\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}")
    if root_entries(temporary) != expected_previous + batch_entries:
        raise RuntimeError("merged ROOT entry count is incorrect")
    temporary.replace(accumulated)


def write_comparison(path: Path, original: dict[str, Any], state: dict[str, Any]) -> None:
    combined = campaign_results(state)
    with path.open("w", newline="") as stream:
        fields = ["key", "before_estimate", "before_rse", "before_neff",
                  "before_signal", "after_estimate", "after_rse", "after_neff",
                  "after_signal", "largest_fraction", "stability_z", "vov",
                  "rse_slope", "vov_slope", "fom_range", "tail_index",
                  "next_history_fraction", "passed", "failed"]
        writer = csv.DictWriter(stream, fields, delimiter="\t")
        writer.writeheader()
        for key in required_keys():
            before = direct_result(original[key])
            after = combined.get(key, {})
            diagnostic = after.get("reliability", {})
            passed, failed = gate(after, state["config"]) if after else (False, ["not_run"])
            writer.writerow({
                "key": key, "before_estimate": before["estimate"],
                "before_rse": before["relative_error"],
                "before_neff": before["effective_sample_size"],
                "before_signal": before["signal_events"],
                "after_estimate": after.get("estimate"),
                "after_rse": after.get("relative_error"),
                "after_neff": after.get("effective_sample_size"),
                "after_signal": after.get("signal_events"),
                "largest_fraction": after.get("largest_event_fraction"),
                "stability_z": after.get("stability_z"),
                "vov": diagnostic.get("variance_of_variance"),
                "rse_slope": diagnostic.get("rse_log_slope"),
                "vov_slope": diagnostic.get("vov_log_slope"),
                "fom_range": diagnostic.get("fom_relative_range"),
                "tail_index": diagnostic.get("tail_index"),
                "next_history_fraction": diagnostic.get("next_largest_history_fraction"),
                "passed": passed, "failed": ",".join(failed),
            })


def display(original: dict[str, Any], state: dict[str, Any], limit: int = 10) -> None:
    combined = campaign_results(state)
    print(f"\nAfter {len(state['batches'])} batches / "
          f"{len(state['batches']) * state['config']['batch_events']} replay draws")
    total_histories = (state["original_events"]
                       + len(state["batches"]) * state["config"]["batch_events"])
    print(f"Combined convergence uses {state['original_events']} original + "
          f"{total_histories - state['original_events']} replay histories")
    locked = state.get("locked_target")
    if locked:
        locked_key = target_key(int(locked["detector"]), str(locked["observable"]))
        locked_passed, locked_failed = gate(combined[locked_key], state["config"])
        print(f"TARGET LOCK: {locked_key} -- "
              f"{'FULL CONVERGENCE PASSED' if locked_passed else 'not complete'}")
        if not locked_passed:
            print("  Remaining gates: " + ", ".join(locked_failed))
    if state["batches"]:
        latest = state["batches"][-1]["proposal"]
        if latest.get("kind") == "targeted":
            latest_key = f"tile_{latest['target_detector']}/{latest['observable']}"
            before_latest = direct_result(original[latest_key])
            after_latest = combined[latest_key]
            print(f"Latest targeted result: {latest_key}  "
                  f"{100*(before_latest['relative_error'] or 0):.2f}% original RSE -> "
                  f"{100*(after_latest['relative_error'] or 0):.2f}% combined RSE "
                  f"using {after_latest.get('targeted_histories', 0)} targeted histories")
        else:
            print("Latest batch was an unbiased physical reliability control; it does not "
                  "replace a targeted tile estimator.")
    targeted_keys = []
    for batch in state["batches"]:
        proposal = batch["proposal"]
        if proposal.get("kind") != "targeted":
            continue
        key = f"tile_{proposal['target_detector']}/{proposal['observable']}"
        if key not in targeted_keys:
            targeted_keys.append(key)
    if targeted_keys:
        target_rse = state["config"]["settings"]["precision"]["relative_target"]
        print("Completed targeted estimators:")
        for key in targeted_keys:
            before = direct_result(original[key])
            after = combined[key]
            precision_status = "RSE PASS" if (
                after.get("relative_error") is not None
                and after["relative_error"] <= target_rse
            ) else "NEEDS MORE"
            print(f"  {key:30s} {100*(before['relative_error'] or 0):6.2f}% -> "
                  f"{100*(after.get('relative_error') or 0):6.2f}% RSE, "
                  f"{after.get('targeted_histories', 0):5d} targeted histories  "
                  f"{precision_status}")
    print("Worst remaining precision rows:")
    print("statistic                         original RSE combined Neff combined RSE  estimate combined/original  status")
    rows = []
    for key in required_keys():
        before = direct_result(original[key])
        after = combined.get(key)
        if after:
            passed, failed = gate(after, state["config"])
            score = after.get("relative_error") or 10.0
            status = "PASS" if passed else (failed[0] if failed else "running")
            rows.append((score, key, before, after, status))
        else:
            rows.append((before["relative_error"] or 10.0, key, before, {}, "not_run"))
    for _, key, before, after, status in sorted(rows, reverse=True)[:limit]:
        def pct(value: float | None) -> str:
            return "-" if value is None else f"{100*value:8.2f}%"
        ratio = (after.get("estimate", 0.0) / before["estimate"]
                 if after and before["estimate"] else None)
        ratio_text = "-" if ratio is None else f"{ratio:12.4f}"
        print(f"{key:34s} {pct(before['relative_error'])} "
              f"{after.get('effective_sample_size', 0):13.1f} "
              f"{pct(after.get('relative_error'))} {ratio_text:>18s}  {status}")
    target_rse = state["config"]["settings"]["precision"]["relative_target"]
    precision_count = sum(
        value.get("estimate", 0.0) != 0.0
        and value.get("relative_error") is not None
        and value["relative_error"] <= target_rse
        for key, value in combined.items() if key in required_keys()
    )
    passed_count = sum(gate(value, state["config"])[0] for value in combined.values())
    print(f"RSE precision passed {precision_count}/{len(required_keys())}; full convergence "
          f"passed {passed_count}/{len(required_keys())} (also requires tail, stability, "
          "and physical-control reliability gates).\n")


def load_or_create(args: argparse.Namespace, original: dict[str, Any],
                   fingerprint: str) -> dict[str, Any]:
    state_path = args.output_dir / "state.json"
    config = config_record(args, fingerprint)
    if state_path.is_file():
        state = json.loads(state_path.read_text())
        saved_config = dict(state.get("config", {}))
        comparable_saved = dict(saved_config)
        comparable_current = dict(config)
        comparable_saved.pop("threads", None)
        comparable_current.pop("threads", None)
        if (state.get("version") != STATE_VERSION
                or comparable_saved != comparable_current):
            raise RuntimeError("existing run configuration differs; choose another --output-dir")
        if saved_config.get("threads") != config["threads"]:
            state["config"]["threads"] = config["threads"]
            atomic_json(state_path, state)
        # Version-2 runs originally accumulated replay moments alone.  Upgrade them
        # in place without rerunning any simulation by prepending the saved original
        # history moments exactly once.
        if not state.get("aggregate_includes_original", False):
            upgraded = original_aggregate(original)
            for key, replay in state["aggregate"].items():
                merge_existing_aggregate(upgraded[key], replay)
            state["aggregate"] = upgraded
            state["original_events"] = original_event_count(original)
            state["aggregate_includes_original"] = True
            atomic_json(state_path, state)
        return state
    state = {
        "version": STATE_VERSION, "created": utc_now(), "updated": utc_now(),
        "config": config, "batches": [],
        "original_events": original_event_count(original),
        "aggregate_includes_original": True,
        "aggregate": original_aggregate(original), "control_aggregate": {
            key: new_merged_stat() for key in original
        }, "training": {}, "training_batches": [], "active": None,
        "target_lock_initialized": True, "locked_target": None,
    }
    atomic_json(state_path, state)
    return state


def save_state(args: argparse.Namespace, state: dict[str, Any]) -> None:
    state["updated"] = utc_now()
    atomic_json(args.output_dir / "state.json", state)


def run_training(args: argparse.Namespace, state: dict[str, Any], detector: int,
                 observable: str, reason: str) -> None:
    key = f"tile_{detector}/{observable}"
    record = state["training"].get(key, {"pilots": 0, "total_successes": 0})
    if record["pilots"] >= args.maximum_pilots_per_target:
        raise RuntimeError(
            f"surface {args.surface} failed to reproduce {key} after "
            f"{record['pilots']} exact-state pilots and {record['total_successes']} successes; "
            "the controller stopped rather than continuing an ineffective density"
        )
    pilot_number = len(state["training_batches"]) + 1
    seed = args.seed + 50_000_000 + pilot_number * 1013
    pilot_source = args.output_dir / "current_pilot_source.root"
    pilot_result = args.output_dir / "current_pilot_result.root"
    pilot_macro = args.output_dir / "current_pilot.mac"
    pilot_log = args.output_dir / "current_pilot.log"
    learning_file = args.output_dir / "learning.root"
    tree_name = learning_tree_name(detector, observable)
    print(f"Training pilot {pilot_number}: exact states for detector {detector} "
          f"{observable}; {reason}")
    output = run_macro(PILOT_BUILD_MACRO, [
        root_string(args.input.resolve()), str(detector), root_string(observable),
        str(args.surface), root_string(pilot_source), str(args.pilot_events),
        str(args.pilot_unique_states), str(seed), str(args.top_histories),
    ], "PILOT_BUILD")
    build_line = next(line for line in output.splitlines() if line.startswith("PILOT_BUILD"))
    print("  " + build_line)
    pilot_macro.write_text(generated_macro(
        state["config"], pilot_source, pilot_result, seed, args.pilot_events
    ))
    start = time.monotonic()
    with pilot_log.open("w") as stream:
        completed = subprocess.run(
            [str(REPO / "build/remoll"), "-t", str(args.threads), "-m", str(pilot_macro)],
            cwd=REPO, stdout=stream, stderr=subprocess.STDOUT,
        )
    elapsed = time.monotonic() - start
    if completed.returncode != 0 or root_entries(pilot_result) != args.pilot_events:
        raise RuntimeError(f"exact-state pilot failed; inspect {pilot_log}")
    analysis_output = run_macro(PILOT_ANALYZE_MACRO, [
        root_string(pilot_source), root_string(pilot_result), str(detector),
        root_string(observable), root_string(learning_file), root_string(tree_name),
    ], "PILOT_ANALYSIS")
    analysis_line = next(
        line for line in analysis_output.splitlines() if line.startswith("PILOT_ANALYSIS")
    )
    fields = dict(item.split("=", 1) for item in analysis_line.split()[1:] if "=" in item)
    total_successes = int(fields["total_successes"])
    successful_states = int(fields["successful_states"])
    record = {
        "pilots": record["pilots"] + 1,
        "total_trials": int(fields["total_trials"]),
        "total_successes": total_successes,
        "successful_states": successful_states,
        "tree": tree_name,
    }
    state["training"][key] = record
    state["training_batches"].append({
        "pilot": pilot_number, "target": key, "events": args.pilot_events,
        "elapsed_seconds": elapsed, "build": build_line, "analysis": analysis_line,
    })
    save_state(args, state)
    print("  " + analysis_line)
    pilot_source.unlink(missing_ok=True); pilot_result.unlink(missing_ok=True)
    pilot_log.unlink(missing_ok=True)


def training_ready(args: argparse.Namespace, state: dict[str, Any], detector: int,
                   observable: str) -> bool:
    record = state["training"].get(f"tile_{detector}/{observable}", {})
    return (record.get("total_successes", 0) >= args.minimum_pilot_successes
            and record.get("successful_states", 0) > 0)


def execute_batch(args: argparse.Namespace, state: dict[str, Any],
                  original: dict[str, Any],
                  selected: tuple[int, str, str] | None = None) -> None:
    batch_number = len(state["batches"]) + 1
    detector, observable, reason = selected or select_target(state, original)
    control = batch_number % args.control_interval == 0
    batch_physical_fraction = 1.0 if control else args.physical_fraction
    seed = args.seed + batch_number * 1009
    source = args.output_dir / "current_source.root"
    macro = args.output_dir / "current_replay.mac"
    result_root = args.output_dir / "current_batch.root"
    stats_file = args.output_dir / "current_stats.tsv"
    log = args.output_dir / "current_remoll.log"
    accumulated = args.output_dir / "replay_accumulated.root"
    proposal = {"kind": "physical_control" if control else "targeted",
                "target_detector": detector, "observable": observable,
                "surface_detector": args.surface, "reason": reason,
                "physical_fraction": batch_physical_fraction, "seed": seed}
    previous_entries = len(state["batches"]) * args.batch_events
    state["active"] = {"batch": batch_number, "proposal": proposal,
                       "previous_entries": previous_entries, "phase": "building"}
    save_state(args, state)

    label = "physical reliability control" if control else "targeted replay"
    print(f"Batch {batch_number}: {label}; detector {detector} {observable}; {reason}")
    output = run_macro(BUILD_MACRO, [
        root_string(args.input.resolve()), str(detector), root_string(observable),
        str(args.surface), root_string(source), str(args.batch_events),
        f"{batch_physical_fraction:.17g}", str(seed), str(args.top_histories),
        root_string(args.output_dir / "learning.root") if not control else root_string(""),
        root_string(learning_tree_name(detector, observable)) if not control else root_string(""),
    ], "SURFACE_BUILD")
    build_line = next(line for line in output.splitlines() if line.startswith("SURFACE_BUILD"))
    print("  " + build_line)
    macro.write_text(generated_macro(state["config"], source, result_root, seed))
    state["active"].update({"phase": "running", "build": build_line,
                            "started_epoch": time.time()})
    save_state(args, state)

    start = time.monotonic()
    with log.open("w") as stream:
        completed = subprocess.run(
            [str(REPO / "build/remoll"), "-t", str(args.threads), "-m", str(macro)],
            cwd=REPO, stdout=stream, stderr=subprocess.STDOUT,
        )
    elapsed = time.monotonic() - start
    if completed.returncode != 0 or root_entries(result_root) != args.batch_events:
        raise RuntimeError(f"remoll batch failed; inspect {log}")
    alerts = []
    for line in log.read_text(errors="replace").splitlines():
        if re.search(r"illegal parameter|fatal exception|command not found|segmentation", line, re.I):
            alerts.append(line[:500])
    if alerts:
        raise RuntimeError(f"remoll reported fatal alerts; inspect {log}: {alerts[:3]}")
    batch_stats = analyze_file(result_root, "replay", args.batch_events, stats_file)
    state["active"].update({"phase": "analyzed", "elapsed_seconds": elapsed})
    save_state(args, state)

    merge_root(accumulated, result_root, previous_entries)
    state["active"]["phase"] = "merged"
    save_state(args, state)
    batch_name = f"batch_{batch_number:04d}"
    for key, value in batch_stats.items():
        merge_stat(state["aggregate"][key], value, batch_name, elapsed, proposal)
        if control:
            merge_stat(state["control_aggregate"][key], value, batch_name,
                       elapsed, proposal)
    state["batches"].append({"name": batch_name, "events": args.batch_events,
                             "elapsed_seconds": elapsed, "proposal": proposal,
                             "build": build_line})
    state["active"] = None
    save_state(args, state)
    result_root.unlink(missing_ok=True)
    stats_file.unlink(missing_ok=True)
    log.unlink(missing_ok=True)
    write_comparison(args.output_dir / "comparison.tsv", original, state)
    display(original, state)


def recover_active(args: argparse.Namespace, state: dict[str, Any],
                   original: dict[str, Any]) -> None:
    """Finish a complete interrupted batch without ever adding it twice."""
    active = state.get("active")
    if not active:
        return
    batch_number = int(active["batch"])
    previous_entries = int(active["previous_entries"])
    expected_after = previous_entries + args.batch_events
    accumulated = args.output_dir / "replay_accumulated.root"
    result_root = args.output_dir / "current_batch.root"
    stats_file = args.output_dir / "current_stats.tsv"
    log = args.output_dir / "current_remoll.log"
    accumulated_entries = root_entries(accumulated) if accumulated.is_file() else 0
    result_entries = root_entries(result_root) if result_root.is_file() else -1

    if accumulated_entries not in (previous_entries, expected_after):
        raise RuntimeError(
            f"cannot recover: cumulative ROOT has {accumulated_entries} entries; "
            f"expected {previous_entries} or {expected_after}"
        )
    if accumulated_entries == previous_entries and result_entries != args.batch_events:
        # The interrupted result is incomplete. It has never been accumulated,
        # Discard only controller-owned temporaries.  Do not call execute_batch
        # here: target selection may have changed after incorporating completed
        # histories, and production must never bypass the pilot-training gate in
        # main().  Returning there preserves the normal select/train/run order.
        result_root.unlink(missing_ok=True)
        stats_file.unlink(missing_ok=True)
        state["active"] = None
        save_state(args, state)
        print(f"Recovering interrupted batch {batch_number}: no complete result was present; "
              "discarded its temporary output and returned it to the training queue.")
        return
    if result_entries != args.batch_events:
        raise RuntimeError(
            "the cumulative file contains the active batch but its temporary result is "
            "missing; retain the directory and inspect state.json before proceeding"
        )

    print(f"Recovering completed batch {batch_number} without rerunning remoll.")
    batch_stats = analyze_file(result_root, "replay", args.batch_events, stats_file)
    if accumulated_entries == previous_entries:
        merge_root(accumulated, result_root, previous_entries)
    elapsed = float(active.get("elapsed_seconds", 0.0))
    if elapsed <= 0.0:
        started = float(active.get("started_epoch", result_root.stat().st_mtime))
        elapsed = max(1.0, result_root.stat().st_mtime - started)
    proposal = active["proposal"]
    batch_name = f"batch_{batch_number:04d}"
    for key, value in batch_stats.items():
        merge_stat(state["aggregate"][key], value, batch_name, elapsed, proposal)
        if proposal.get("kind") == "physical_control":
            merge_stat(state["control_aggregate"][key], value, batch_name,
                       elapsed, proposal)
    state["batches"].append({"name": batch_name, "events": args.batch_events,
                             "elapsed_seconds": elapsed, "proposal": proposal,
                             "build": active.get("build", "recovered")})
    state["active"] = None
    save_state(args, state)
    result_root.unlink(missing_ok=True)
    stats_file.unlink(missing_ok=True)
    log.unlink(missing_ok=True)
    write_comparison(args.output_dir / "comparison.tsv", original, state)
    display(original, state)


def validate(args: argparse.Namespace) -> None:
    missing = [path for path in (args.input, BUILD_MACRO, PILOT_BUILD_MACRO,
                                  PILOT_ANALYZE_MACRO, ANALYZE_MACRO,
                                  REPO / "build/remoll", REPO / "build/reroot")
               if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing required files: {missing}")
    if shutil.which("hadd") is None:
        raise RuntimeError("ROOT hadd was not found")
    with tempfile.TemporaryDirectory(prefix="surface_replay_check_") as directory:
        temporary = Path(directory)
        stats = analyze_file(args.input, "original", 0, temporary / "original.tsv")
        worst = max(
            (direct_result(value)["relative_error"] or 0.0, key)
            for key, value in stats.items() if key.startswith("tile_")
        )
        source = temporary / "source.root"
        run_macro(BUILD_MACRO, [
            root_string(args.input.resolve()), worst[1].split("/")[0].removeprefix("tile_"),
            root_string(worst[1].split("/")[1]), str(args.surface), root_string(source),
            "10", f"{args.physical_fraction:.17g}", str(args.seed), "50",
        ], "SURFACE_BUILD")
        if root_entries(source) != 10:
            raise RuntimeError("generated source did not contain 10 validation entries")
        if named_tree_entries(source, "density") <= 10:
            raise RuntimeError("generated source did not retain the complete density bank")
    print("CHECK PASSED: old-schema analysis, all 62 tallies, multi-species density, "
          "weighted source schema, remoll binaries, and hadd are available.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batch-events", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1,
                        help="Geant4 worker threads for pilot and production replay")
    parser.add_argument("--surface", type=int, default=35)
    parser.add_argument("--physical-fraction", type=float, default=0.20)
    parser.add_argument("--control-interval", type=int, default=5,
                        help="every Nth batch uses q=p for MCNP-style reliability checks")
    parser.add_argument("--top-histories", type=int, default=400)
    parser.add_argument("--pilot-events", type=int, default=1000)
    parser.add_argument("--pilot-unique-states", type=int, default=200)
    parser.add_argument("--minimum-pilot-successes", type=int, default=20)
    parser.add_argument("--maximum-pilots-per-target", type=int, default=3)
    parser.add_argument("--physics-list", default="QGSP_BERT_HP")
    parser.add_argument("--seed", type=int, default=2202701)
    parser.add_argument("--target-rse", type=float, default=0.05)
    parser.add_argument("--maximum-single-fraction", type=float, default=0.02)
    parser.add_argument("--stability-sigma", type=float, default=3.0)
    parser.add_argument("--minimum-batches", type=int, default=10)
    parser.add_argument("--minimum-nonzero", type=int, default=201)
    parser.add_argument("--maximum-vov", type=float, default=0.10)
    parser.add_argument("--rse-slope-minimum", type=float, default=-0.75)
    parser.add_argument("--rse-slope-maximum", type=float, default=-0.25)
    parser.add_argument("--vov-slope-minimum", type=float, default=-1.50)
    parser.add_argument("--vov-slope-maximum", type=float, default=-0.50)
    parser.add_argument("--maximum-fom-range", type=float, default=0.30)
    parser.add_argument("--minimum-tail-index", type=float, default=3.0)
    parser.add_argument("--maximum-next-fraction", type=float, default=0.02)
    parser.add_argument("--max-batches", type=int, default=100)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    args.input = args.input.resolve()
    args.output_dir = args.output_dir.resolve()
    if (args.batch_events <= 0 or args.threads <= 0 or args.max_batches <= 0
            or args.control_interval < 2 or args.pilot_events <= 0
            or args.pilot_unique_states <= 0 or args.minimum_pilot_successes <= 0
            or args.maximum_pilots_per_target <= 0):
        parser.error("batch events, threads, and max batches must be positive")
    if not 0.0 < args.physical_fraction <= 1.0:
        parser.error("physical fraction must be in (0,1]")
    return args


def main() -> int:
    args = parse_args()
    if args.check:
        validate(args)
        return 0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fingerprint = sha256(args.input)
    original_path = args.output_dir / "original_stats.tsv"
    if not original_path.is_file():
        original = analyze_file(args.input, "original", 0, original_path)
    else:
        original = analyze_file(args.input, "original", 0, original_path)
    state = load_or_create(args, original, fingerprint)
    recover_active(args, state, original)
    initialize_target_lock(args, state, original)
    write_comparison(args.output_dir / "comparison.tsv", original, state)
    display(original, state)
    while len(state["batches"]) < args.max_batches:
        results = campaign_results(state)
        if results and all(gate(results[key], state["config"])[0]
                           for key in required_keys()):
            print("All ShowerMax tile and group convergence gates passed.")
            return 0
        selected = locked_selection(args, state, original)
        next_batch = len(state["batches"]) + 1
        control = next_batch % args.control_interval == 0
        if not control and not training_ready(args, state, selected[0], selected[1]):
            run_training(args, state, *selected)
            continue
        execute_batch(args, state, original, selected)
    print(f"Stopped at --max-batches={args.max_batches}; the same output directory is "
          "resumable. Rerun with a larger --max-batches value to continue.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted; the original ROOT file was not modified.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
