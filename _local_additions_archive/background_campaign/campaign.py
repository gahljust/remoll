#!/usr/bin/env python3
"""Run reproducible remoll batches without streaming remoll output to a terminal."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import tomllib
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_CAMPAIGN = HERE / "campaign.toml"
DETECTORS_FILE = HERE / "detectors.toml"
STATE_VERSION = 1
SUMMARY_VERSION = 3
ANALYSIS_VERSION = 2
TAIL_HISTORY_COUNT = 201
FAILED_LOG_LIMIT = 5_000_000
ACTIVE_LOG_LIMIT = 5_000_000


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def find_repo_root() -> Path:
    for path in (HERE, *HERE.parents):
        if (path / "geometry/mollerMother.gdml").is_file():
            return path
    raise RuntimeError("cannot locate remoll repository root")


REPO = find_repo_root()


def load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def resolve_repo_path(value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO / path


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def detector_groups() -> tuple[dict[str, list[int]], set[int], set[int]]:
    definitions = load_toml(DETECTORS_FILE)
    main = definitions["main_detector"]
    groups: dict[str, list[int]] = {}
    for ring, suffixes in main["rings"].items():
        ids = []
        for base in (main["beam_facing_base"], main["far_facing_base"]):
            for segment in range(1, main["segments"] + 1):
                for suffix in suffixes:
                    ids.append(base + segment * main["segment_multiplier"] + suffix)
        groups[ring] = sorted(ids)

    showermax = definitions["showermax"]
    all_showermax = list(
        range(showermax["first_id"], showermax["last_id"] + 1, showermax["step"])
    )
    groups["showermax_open"] = showermax["open"]
    groups["showermax_closed"] = showermax["closed"]
    groups["showermax_transition"] = showermax["transition"]
    if sorted(sum((groups[name] for name in (
        "showermax_open", "showermax_closed", "showermax_transition"
    )), [])) != all_showermax:
        raise ValueError("ShowerMax region groups do not partition all 28 stack planes")

    diagnostic = set(definitions["diagnostic"]["gem"])
    diagnostic.update(definitions["diagnostic"]["ring_virtual"])
    precision = set().union(*groups.values())
    return groups, precision, diagnostic


def geometry_detector_ids() -> set[int]:
    ids: set[int] = set()
    paths = list((REPO / "geometry/detector/ThinQuartz/DetectorArray").glob("*.gdml"))
    paths.append(REPO / "geometry/mollerParallel.gdml")
    for path in paths:
        text = re.sub(r"<!--.*?-->", "", path.read_text(), flags=re.DOTALL)
        for value in re.findall(
            r'<auxiliary\s+auxtype="DetNo"\s+auxvalue="(\d+)"', text
        ):
            ids.add(int(value))
    return ids


def active_lines(path: Path) -> list[str]:
    lines = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            lines.append(line)
    return lines


def energy_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def all_configurations(campaign: dict[str, Any]) -> list[dict[str, Any]]:
    configs = [dict(item) for item in campaign.get("configuration", [])]
    for matrix in campaign.get("configuration_matrix", []):
        for target, sieve, interaction, energy in itertools.product(
            matrix["targets"], matrix["sieves"], matrix["interactions"],
            matrix["energies_gev"],
        ):
            values = {
                "target": target,
                "sieve": sieve,
                "interaction": interaction,
                "energy_gev": float(energy),
                "energy_label": energy_label(float(energy)),
            }
            configs.append({
                "name": matrix["name_format"].format(**values),
                "enabled": bool(matrix.get("enabled", True)),
                "priority": int(matrix.get("priority", 999)),
                "current_uA": float(matrix["current_uA"]),
                **values,
            })

    resolved = []
    for item in configs:
        config = dict(item)
        target = campaign.get("target", {}).get(config.get("target"))
        sieve = campaign.get("sieve", {}).get(config.get("sieve"))
        interaction = campaign.get("interaction", {}).get(config.get("interaction"))
        if target is not None:
            config["target_settings"] = target
        if sieve is not None:
            config["sieve_settings"] = sieve
        if interaction is not None:
            config["interaction_settings"] = interaction
        resolved.append(config)
    return resolved


def config_by_name(campaign: dict[str, Any], name: str) -> dict[str, Any]:
    matches = [item for item in all_configurations(campaign) if item.get("name") == name]
    if len(matches) != 1:
        raise ValueError(f"configuration {name!r} was not found exactly once")
    return matches[0]


def check_campaign(path: Path, require_enabled: bool = False) -> list[str]:
    errors: list[str] = []
    try:
        campaign = load_toml(path)
    except Exception as error:
        return [f"cannot read campaign: {error}"]

    for key in (
        "title", "output_directory", "base_macro", "batch_events", "threads",
        "precision", "reliability",
    ):
        if key not in campaign:
            errors.append(f"missing campaign setting: {key}")
    for key in ("batch_events", "threads", "max_batches_per_configuration"):
        if int(campaign.get(key, 0)) <= 0:
            errors.append(f"{key} must be positive")
    adaptive = campaign.get("adaptive", {})
    if int(adaptive.get("ordinary_control_interval", 0)) < 2:
        errors.append("adaptive ordinary_control_interval must be at least 2")
    if int(adaptive.get("minimum_trials_before_reject", 0)) <= 0:
        errors.append("adaptive minimum_trials_before_reject must be positive")
    if float(adaptive.get("rejection_cost_ratio", 0.0)) <= 1.0:
        errors.append("adaptive rejection_cost_ratio must be greater than 1")
    reliability = campaign.get("reliability", {})
    required_reliability = (
        "minimum_ordinary_batches", "minimum_nonzero_histories",
        "maximum_variance_of_variance", "rse_slope_minimum", "rse_slope_maximum",
        "vov_slope_minimum", "vov_slope_maximum", "maximum_fom_relative_range",
        "minimum_tail_index", "maximum_next_history_fraction",
    )
    for key in required_reliability:
        if key not in reliability:
            errors.append(f"missing reliability setting: {key}")

    base_macro = resolve_repo_path(str(campaign.get("base_macro", "")))
    if not base_macro.is_file():
        errors.append(f"base macro not found: {base_macro}")
    else:
        lines = active_lines(base_macro)
        if any(line.startswith("/run/beamOn ") for line in lines):
            errors.append("base macro must not contain an active /run/beamOn")
        if not any(line.startswith("/run/initialize") for line in lines):
            errors.append("base macro has no active /run/initialize")

    names: set[str] = set()
    enabled = 0
    try:
        configs = all_configurations(campaign)
    except Exception as error:
        errors.append(f"cannot expand configuration inventory: {error}")
        configs = []
    for config in configs:
        name = config.get("name", "")
        if not re.fullmatch(r"[a-z0-9][a-z0-9_.-]*", name):
            errors.append(f"invalid configuration name: {name!r}")
        if name in names:
            errors.append(f"duplicate configuration name: {name}")
        names.add(name)
        enabled += bool(config.get("enabled", False))
        for reference in ("target_settings", "sieve_settings", "interaction_settings"):
            if reference not in config:
                errors.append(f"{name}: unknown {reference.removesuffix('_settings')}")
        target_settings = config.get("target_settings", {})
        if not all(key in target_settings for key in ("z_min_mm", "z_max_mm")):
            errors.append(f"{name}: target requires z_min_mm and z_max_mm")
        elif float(target_settings["z_max_mm"]) <= float(target_settings["z_min_mm"]):
            errors.append(f"{name}: target z bounds are reversed or empty")
        for reference in ("target_settings", "sieve_settings"):
            macro_value = config.get(reference, {}).get("macro", "")
            if macro_value and not resolve_repo_path(str(macro_value)).is_file():
                errors.append(f"{name}: macro not found: {resolve_repo_path(str(macro_value))}")
        if float(config.get("energy_gev", 0.0)) <= 0.0:
            errors.append(f"{name}: energy_gev must be positive")
        if float(config.get("current_uA", 0.0)) <= 0.0:
            errors.append(f"{name}: current_uA must be positive")

    if require_enabled and enabled == 0:
        errors.append("no configurations are enabled")
    expected = int(campaign.get("expected_configurations", len(configs)))
    if len(configs) != expected:
        errors.append(f"configuration inventory has {len(configs)} entries; expected {expected}")

    remoll = REPO / "build/remoll"
    if not os.access(remoll, os.X_OK):
        errors.append(f"remoll is not executable: {remoll}")

    try:
        _, precision, diagnostic = detector_groups()
        missing = (precision | diagnostic) - geometry_detector_ids()
        if missing:
            errors.append(f"detector IDs absent from current GDML: {sorted(missing)}")
    except Exception as error:
        errors.append(f"detector definitions are invalid: {error}")
    return errors


def output_directory(campaign: dict[str, Any]) -> Path:
    return resolve_repo_path(campaign["output_directory"])


def state_path(campaign: dict[str, Any]) -> Path:
    return output_directory(campaign) / "state.json"


def read_state(campaign: dict[str, Any]) -> dict[str, Any]:
    path = state_path(campaign)
    if path.is_file():
        return json.loads(path.read_text())
    return {
        "version": STATE_VERSION,
        "campaign": campaign["title"],
        "control": "ready",
        "active": None,
        "updated": utc_now(),
    }


def write_state(campaign: dict[str, Any], state: dict[str, Any]) -> None:
    state["updated"] = utc_now()
    atomic_json(state_path(campaign), state)


STRIP_COMMANDS = (
    "/run/beamOn",
    "/remoll/filename",
    "/remoll/seed",
    "/remoll/SD/",
)


def proposal_commands(proposal: dict[str, Any]) -> list[str]:
    if proposal.get("kind") != "mixture":
        return []
    axis = proposal["axis"]
    fraction = proposal["physical_fraction"]
    if axis == "theta":
        prefix, minimum, maximum, unit = "/remoll/bias/thcom", "min", "max", "deg"
    elif axis == "phi":
        prefix, minimum, maximum, unit = "/remoll/bias/phi", "min", "max", "deg"
    elif axis == "beamp":
        prefix, minimum, maximum, unit = "/remoll/bias/beamp", "min", "max", "GeV"
    elif axis == "vertexz":
        return [
            "/remoll/bias/vertexz/mode mixture",
            f"/remoll/bias/vertexz/minFraction {proposal['minimum']:.12g}",
            f"/remoll/bias/vertexz/maxFraction {proposal['maximum']:.12g}",
            f"/remoll/bias/vertexz/physicalFraction {fraction:.12g}",
        ]
    elif axis == "outgoinge":
        return [
            "/remoll/bias/outgoinge/mode mixture",
            f"/remoll/bias/outgoinge/minFraction {proposal['minimum']:.12g}",
            f"/remoll/bias/outgoinge/maxFraction {proposal['maximum']:.12g}",
            f"/remoll/bias/outgoinge/physicalFraction {fraction:.12g}",
        ]
    else:
        raise ValueError(f"unknown proposal axis: {axis}")
    return [
        f"{prefix}/mode mixture",
        f"{prefix}/{minimum} {proposal['minimum']:.12g} {unit}",
        f"{prefix}/{maximum} {proposal['maximum']:.12g} {unit}",
        f"{prefix}/physicalFraction {fraction:.12g}",
    ]


def generated_macro_text(source: Path, config: dict[str, Any], proposal: dict[str, Any],
                         root_file: Path, seed: int, events: int) -> str:
    _, precision, diagnostic = detector_groups()
    retained = []
    for raw in source.read_text().splitlines():
        stripped = raw.strip()
        if stripped and not stripped.startswith("#"):
            if any(stripped.startswith(command) for command in STRIP_COMMANDS):
                continue
            if stripped.startswith("/control/execute") and "enable_sm_vplanes" in stripped:
                continue
        retained.append(raw)

    interaction = config["interaction_settings"]
    theta_command = interaction["theta_command"]
    theta_prefix = "thcom" if theta_command == "thcom" else "th"
    interval = max(1, events // 100)
    lines = retained + [
        "",
        "# Physics configuration",
        f"/control/execute {config['target_settings']['macro']}",
    ]
    if config["sieve_settings"].get("macro"):
        lines.append(f"/control/execute {config['sieve_settings']['macro']}")
    lines.extend([
        f"/remoll/evgen/set {interaction['generator']}",
        f"/remoll/evgen/{theta_prefix}min {interaction['theta_min_deg']:.12g} deg",
        f"/remoll/evgen/{theta_prefix}max {interaction['theta_max_deg']:.12g} deg",
        "/remoll/evgen/phmin 0 deg",
        "/remoll/evgen/phmax 360 deg",
        f"/remoll/beamene {config['energy_gev']:.12g} GeV",
        f"/remoll/beamcurr {config['current_uA']:.12g} microampere",
    ])
    if interaction.get("apply_screening"):
        lines.append("/remoll/evgen/elastic/applyScreening true")
    lines.extend(proposal_commands(proposal))
    lines.extend([
        "",
        "# Detector-plane campaign controls",
        f"/remoll/interval {interval}",
        "/remoll/SD/disable_all",
    ])
    for detector in sorted(precision | diagnostic):
        lines.append(f"/remoll/SD/enable {detector}")
        lines.append(f"/remoll/SD/detect surfacehits {detector}")
        lines.append(f"/remoll/SD/detect lowenergyneutral {detector}")
    lines.extend((
        f"/remoll/seed {seed}",
        f"/remoll/filename {root_file}",
        f"/run/beamOn {events}",
        "",
    ))
    return "\n".join(lines)


def next_batch_number(directory: Path) -> int:
    numbers = []
    for path in directory.glob("batch_*.json"):
        match = re.fullmatch(r"batch_(\d+)\.json", path.name)
        if match:
            numbers.append(int(match.group(1)))
    return max(numbers, default=0) + 1


def random_seed() -> int:
    return int.from_bytes(os.urandom(4), "little") & 0x7FFFFFFF


def disk_used_gb(path: Path) -> float:
    if not path.exists():
        return 0.0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file()) / 1e9


def on_ac_power() -> bool:
    if sys.platform != "darwin" or not shutil.which("pmset"):
        return True
    result = subprocess.run(["pmset", "-g", "batt"], text=True, capture_output=True)
    return "AC Power" in result.stdout


def root_integrity(path: Path) -> tuple[int | None, bool]:
    reroot = REPO / "build/reroot"
    expression = (
        'printf("CAMPAIGN_INTEGRITY %lld %d\\n",T ? T->GetEntries() : -1,'
        'T && T->GetBranch("ev") && T->GetBranch("part") && '
        'T->GetBranch("hit") && T->GetBranch("rate") && T->GetBranch("bm"));'
    )
    result = subprocess.run(
        [str(reroot), "-l", "-b", str(path), "-q", "-e", expression],
        cwd=REPO, text=True, capture_output=True,
    )
    match = re.search(r"CAMPAIGN_INTEGRITY (-?\d+) ([01])", result.stdout + result.stderr)
    return (int(match.group(1)), bool(int(match.group(2)))) if match else (None, False)


def extract_alerts(log: Path) -> list[str]:
    alerts = []
    pattern = re.compile(r"warning|error|fatal|exception|abort|segmentation|cannot|failed", re.I)
    ignored = (
        "high precision neutron", "Some can be safely ignore", "GDML file validation",
        "warnings that cannot be avoided", "G4Material WARNING: duplicate name",
    )
    for line in log.read_text(errors="replace").splitlines():
        if pattern.search(line) and not any(text in line for text in ignored):
            if line not in alerts:
                alerts.append(line[:1000])
    return alerts[:100]


def cap_file(path: Path, limit: int) -> None:
    if not path.is_file() or path.stat().st_size <= limit:
        return
    with path.open("rb") as stream:
        stream.seek(-limit, os.SEEK_END)
        tail = stream.read()
    path.write_bytes(b"[earlier output removed]\n" + tail)


def analyze_root(path: Path, config: dict[str, Any]) -> tuple[
        int, dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    macro = HERE / "analyze_batch.C"
    interaction = config["interaction_settings"]
    target = config["target_settings"]
    invocation = (
        f'{macro}("{path}",{config["energy_gev"]:.12g},'
        f'{interaction["theta_min_deg"]:.12g},{interaction["theta_max_deg"]:.12g},'
        f'{target["z_min_mm"]:.12g},{target["z_max_mm"]:.12g})'
    )
    result = subprocess.run(
        [str(REPO / "build/reroot"), "-l", "-b", "-q", invocation],
        cwd=REPO, text=True, capture_output=True,
    )
    if result.returncode != 0 or "CAMPAIGN_ERROR" in result.stdout:
        raise RuntimeError(f"ROOT analysis failed for {path}")
    events = 0
    stats: dict[str, dict[str, float]] = {}
    cells: dict[str, dict[str, float]] = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if fields[:1] == ["CAMPAIGN_BATCH"]:
            events = int(fields[1])
        elif fields[:1] == ["CAMPAIGN_STAT"] and len(fields) == 11:
            key = "/".join(fields[1:4])
            stats[key] = {
                "signal_events": int(fields[4]),
                "sum": float(fields[5]),
                "sum2": float(fields[6]),
                "sum3": float(fields[7]),
                "sum4": float(fields[8]),
                "maximum": float(fields[9]),
                "largest": (
                    [] if fields[10] == "-"
                    else [float(item) for item in fields[10].split(",")]
                ),
            }
        elif fields[:1] == ["CAMPAIGN_CELL"] and len(fields) == 12:
            key = "/".join((fields[1], fields[2], fields[5], fields[6]))
            cells[key] = {
                "minimum": float(fields[3]),
                "maximum_bound": float(fields[4]),
                "sample_events": int(fields[7]),
                "signal_events": int(fields[8]),
                "sum": float(fields[9]),
                "sum2": float(fields[10]),
                "maximum": float(fields[11]),
            }
    if events <= 0 or not stats:
        raise RuntimeError(f"ROOT analysis returned no statistics for {path}")
    if not cells:
        raise RuntimeError(f"ROOT analysis returned no adaptive cells for {path}")
    return events, stats, cells


def batch_variance(events: int, value: dict[str, float]) -> float:
    if events <= 1:
        return 0.0
    return max(0.0, (value["sum2"] - value["sum"] ** 2 / events) / (events - 1))


def statistic_result(events: int, value: dict[str, Any]) -> dict[str, float | bool | None]:
    mean = value["sum"] / events if events else 0.0
    standard_error = math.sqrt(value["variance_numerator"]) / events if events else 0.0
    relative_error = standard_error / abs(mean) if mean else None
    largest_fraction = value["maximum"] / abs(value["sum"]) if value["sum"] else None
    largest_batch_fraction = (
        value["maximum_batch_sum"] / abs(value["sum"]) if value["sum"] else None
    )
    history = value["history"]
    recent_count = 3
    stable = False
    stability_z = None
    if len(history) > recent_count:
        older = history[:-recent_count]
        recent = history[-recent_count:]
        old_events = sum(item["events"] for item in older)
        recent_events = sum(item["events"] for item in recent)
        old_mean = sum(item["sum"] for item in older) / old_events
        recent_mean = sum(item["sum"] for item in recent) / recent_events
        old_se = math.sqrt(sum(item["variance_numerator"] for item in older)) / old_events
        recent_se = math.sqrt(sum(item["variance_numerator"] for item in recent)) / recent_events
        scale = math.hypot(old_se, recent_se)
        stability_z = abs(recent_mean - old_mean) / scale if scale else (
            0.0 if recent_mean == old_mean else None
        )
        stable = stability_z is not None
    return {
        "estimate": mean,
        "standard_error": standard_error,
        "relative_error": relative_error,
        "largest_event_fraction": largest_fraction,
        "largest_batch_fraction": largest_batch_fraction,
        "stability_z": stability_z,
        "stable": stable,
    }


def apply_completion_states(campaign: dict[str, Any], result: dict[str, Any]) -> None:
    if not result["batches"]:
        result["state"] = "pending"
        return
    precision = campaign["precision"]
    target = float(precision["relative_target"])
    max_event = float(precision["maximum_single_event_fraction"])
    stability_sigma = float(precision["stability_sigma"])
    min_ordinary = int(precision["minimum_ordinary_batches"])
    min_total = int(precision["minimum_total_batches"])
    negligible_fraction = float(precision["negligible_fraction"])
    minimum_bound_events = int(precision["minimum_signal_events_for_bound"])
    at_budget = (
        len(result["batches"]) >= int(campaign["max_batches_per_configuration"])
        or result["events"] >= int(campaign["max_events_per_configuration"])
    )
    enough_batches = (
        result["ordinary_batches"] >= min_ordinary
        and len(result["batches"]) >= min_total
    )
    confidence = float(precision["confidence"])
    z_value = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    required = []
    groups = detector_groups()[0]
    detector_to_group = {
        detector: group for group, detectors in groups.items() for detector in detectors
    }
    for group in groups:
        if group.startswith("ring"):
            value = result["statistics"].get(f"{group}/all/rate_energy")
            if value is not None:
                value["state"] = "not_relevant"
                value["reason"] = "rings use crossing rate, not rate times energy"
    required_tallies: list[tuple[str, float]] = []
    for observable in ("rate", "rate_energy"):
        observable_groups = [group for group in groups
                             if observable == "rate" or group.startswith("showermax_")]
        group_reference = max(
            (result["statistics"].get(f"{group}/all/{observable}", {}).get("estimate", 0.0)
             for group in observable_groups), default=0.0,
        )
        required_tallies.extend(
            (f"{group}/all/{observable}", group_reference) for group in observable_groups
        )
        for detector, parent in detector_to_group.items():
            if observable == "rate_energy" and not parent.startswith("showermax_"):
                continue
            parent_estimate = result["statistics"].get(
                f"{parent}/all/{observable}", {}
            ).get("estimate", 0.0)
            required_tallies.append((f"tile_{detector}/all/{observable}", parent_estimate))

    for key, reference in required_tallies:
        value = result["statistics"].get(key, {})
        estimate = value.get("estimate", 0.0)
        error = value.get("standard_error", 0.0)
        relative = value.get("relative_error")
        largest = value.get("largest_event_fraction")
        stable = (
            value.get("stable", False)
            and value.get("stability_z") is not None
            and value["stability_z"] <= stability_sigma
        )
        upper = estimate + z_value * error
        reliability = value.get("reliability", {})
        if not enough_batches:
            state, reason = "running", "ordinary pilot or minimum batch count incomplete"
        elif estimate == 0.0:
            state, reason = "rare_tail", "no signal-producing event observed"
        elif value.get("signal_events", 0) >= minimum_bound_events \
                and reference > 0.0 and upper <= negligible_fraction * reference:
            state, reason = "negligible", "confidence upper bound is below importance threshold"
        elif relative is not None and relative <= target and largest is not None \
                and largest <= max_event and stable and reliability.get("passed", False):
            state, reason = "measured", "precision and tally reliability checks passed"
        elif relative is not None and relative > target:
            state, reason = "running", "relative error above target"
        elif largest is not None and largest > max_event:
            state, reason = "rare_tail", "single-event contribution is too large"
        elif not reliability.get("passed", False):
            failed = ", ".join(reliability.get("failed_checks", []))
            state, reason = "running", f"tally reliability checks incomplete: {failed}"
        else:
            state, reason = "running", "recent batches are not yet stable"
        if at_budget and state not in ("measured", "negligible"):
            state, reason = "statistics_limited", "configuration budget reached"
        value["state"] = state
        value["reason"] = reason
        value["confidence_upper"] = upper
        required.append(state)
    result["state"] = (
        "measured" if required and all(state in ("measured", "negligible") for state in required)
        else "statistics_limited" if at_budget
        else "running"
    )


def new_merged_stat() -> dict[str, Any]:
    return {
        "signal_events": 0, "sum": 0.0, "sum2": 0.0, "sum3": 0.0,
        "sum4": 0.0, "maximum": 0.0, "largest": [],
        "maximum_batch_sum": 0.0, "variance_numerator": 0.0, "history": [],
    }


def merge_largest(first: list[float], second: list[float]) -> list[float]:
    return sorted((*first, *second), reverse=True)[:TAIL_HISTORY_COUNT]


def linear_slope(points: list[tuple[float, float]]) -> float | None:
    if len(points) < 3:
        return None
    x_mean = sum(point[0] for point in points) / len(points)
    y_mean = sum(point[1] for point in points) / len(points)
    denominator = sum((point[0] - x_mean) ** 2 for point in points)
    if denominator == 0.0:
        return None
    return sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in points
    ) / denominator


def central_moments(events: int, value: dict[str, Any]) -> tuple[float, float]:
    if events <= 0:
        return 0.0, 0.0
    mean = value["sum"] / events
    second = max(0.0, value["sum2"] - value["sum"] ** 2 / events)
    fourth = (
        value["sum4"] - 4.0 * mean * value["sum3"]
        + 6.0 * mean ** 2 * value["sum2"]
        - 4.0 * mean ** 3 * value["sum"] + events * mean ** 4
    )
    return second, max(0.0, fourth)


def reliability_result(campaign: dict[str, Any], value: dict[str, Any]) -> dict[str, Any]:
    settings = campaign["reliability"]
    history = value["history"]
    checkpoints = []
    events = 0
    elapsed = 0.0
    cumulative = {"sum": 0.0, "sum2": 0.0, "sum3": 0.0, "sum4": 0.0}
    for item in history:
        events += item["events"]
        elapsed += item["elapsed_seconds"]
        for field in cumulative:
            cumulative[field] += item[field]
        second, fourth = central_moments(events, cumulative)
        variance = second / (events - 1) if events > 1 else 0.0
        mean = cumulative["sum"] / events if events else 0.0
        rse = math.sqrt(variance / events) / abs(mean) if mean and events > 1 else None
        vov = fourth / second ** 2 - 1.0 / events if second > 0.0 else None
        fom = 1.0 / (rse ** 2 * elapsed) if rse and elapsed > 0.0 else None
        checkpoints.append({
            "events": events, "estimate": mean, "relative_error": rse,
            "variance_of_variance": max(0.0, vov) if vov is not None else None,
            "figure_of_merit": fom, "elapsed_seconds": elapsed,
        })

    positive_rse = [
        (math.log(point["events"]), math.log(point["relative_error"]))
        for point in checkpoints if point["relative_error"] not in (None, 0.0)
    ]
    positive_vov = [
        (math.log(point["events"]), math.log(point["variance_of_variance"]))
        for point in checkpoints if point["variance_of_variance"] not in (None, 0.0)
    ]
    rse_slope = linear_slope(positive_rse)
    vov_slope = linear_slope(positive_vov)
    foms = [point["figure_of_merit"] for point in checkpoints
            if point["figure_of_merit"] is not None]
    recent_foms = foms[len(foms) // 2:]
    fom_relative_range = (
        (max(recent_foms) - min(recent_foms)) / (sum(recent_foms) / len(recent_foms))
        if len(recent_foms) >= 2 and sum(recent_foms) > 0.0 else None
    )

    tail = value["largest"]
    tail_index = None
    if len(tail) >= TAIL_HISTORY_COUNT and tail[TAIL_HISTORY_COUNT - 1] > 0.0:
        threshold = tail[TAIL_HISTORY_COUNT - 1]
        logarithms = [math.log(item / threshold) for item in tail[:TAIL_HISTORY_COUNT - 1]]
        if sum(logarithms) > 0.0:
            tail_index = (TAIL_HISTORY_COUNT - 1) / sum(logarithms)

    next_history_fraction = None
    next_history_rse = None
    if events > 1 and value["sum"] != 0.0:
        maximum = value["maximum"]
        original_mean = value["sum"] / events
        next_mean = (value["sum"] + maximum) / (events + 1)
        next_history_fraction = abs(next_mean - original_mean) / abs(original_mean)
        augmented = {
            "sum": value["sum"] + maximum,
            "sum2": value["sum2"] + maximum ** 2,
            "sum3": value["sum3"] + maximum ** 3,
            "sum4": value["sum4"] + maximum ** 4,
        }
        second, _ = central_moments(events + 1, augmented)
        variance = second / events
        next_history_rse = math.sqrt(variance / (events + 1)) / abs(next_mean)

    ordinary_batches = len(history)
    final = checkpoints[-1] if checkpoints else {}
    checks = {
        "minimum_ordinary_batches": ordinary_batches >= int(settings["minimum_ordinary_batches"]),
        "minimum_nonzero_histories": value["signal_events"] >= int(settings["minimum_nonzero_histories"]),
        "variance_of_variance": (
            final.get("variance_of_variance") is not None
            and final["variance_of_variance"] <= float(settings["maximum_variance_of_variance"])
        ),
        "relative_error_scaling": (
            rse_slope is not None
            and float(settings["rse_slope_minimum"]) <= rse_slope
            <= float(settings["rse_slope_maximum"])
        ),
        "variance_of_variance_scaling": (
            vov_slope is not None
            and float(settings["vov_slope_minimum"]) <= vov_slope
            <= float(settings["vov_slope_maximum"])
        ),
        "figure_of_merit_stability": (
            fom_relative_range is not None
            and fom_relative_range <= float(settings["maximum_fom_relative_range"])
        ),
        "high_score_tail": (
            tail_index is not None and tail_index >= float(settings["minimum_tail_index"])
        ),
        "next_largest_history": (
            next_history_fraction is not None
            and next_history_fraction <= float(settings["maximum_next_history_fraction"])
        ),
    }
    return {
        "ordinary_events": events,
        "ordinary_batches": ordinary_batches,
        "ordinary_signal_events": value["signal_events"],
        "ordinary_estimate": final.get("estimate"),
        "ordinary_relative_error": final.get("relative_error"),
        "variance_of_variance": final.get("variance_of_variance"),
        "rse_log_slope": rse_slope,
        "vov_log_slope": vov_slope,
        "figure_of_merit": final.get("figure_of_merit"),
        "fom_relative_range": fom_relative_range,
        "tail_index": tail_index,
        "tail_scores": len(tail),
        "next_largest_history_fraction": next_history_fraction,
        "next_largest_history_rse": next_history_rse,
        "checks": checks,
        "passed": all(checks.values()),
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "checkpoints": checkpoints,
    }


def proposal_signature(proposal: dict[str, Any]) -> str:
    return "|".join((
        str(proposal["axis"]), str(proposal["cell"]),
        str(proposal["target_statistic"]),
    ))


def add_performance_sample(target: dict[str, Any], events: int,
                           variance_numerator: float, elapsed: float,
                           signal_events: int) -> None:
    target["batches"] = target.get("batches", 0) + 1
    target["events"] = target.get("events", 0) + events
    target["elapsed_seconds"] = target.get("elapsed_seconds", 0.0) + elapsed
    target["variance_numerator"] = (
        target.get("variance_numerator", 0.0) + variance_numerator
    )
    target["signal_events"] = target.get("signal_events", 0) + signal_events


def variance_cost(value: dict[str, Any]) -> float | None:
    if value.get("events", 0) <= 0 or value.get("elapsed_seconds", 0.0) <= 0.0:
        return None
    return (
        value["variance_numerator"] * value["elapsed_seconds"]
        / (value["events"] ** 2)
    )


def rebuild_summary(campaign: dict[str, Any], only: str | None = None) -> dict[str, Any]:
    summary_path = output_directory(campaign) / "summary.json"
    if only and summary_path.is_file():
        summary = json.loads(summary_path.read_text())
        summary["updated"] = utc_now()
    else:
        summary = {
            "version": SUMMARY_VERSION,
            "campaign": campaign["title"],
            "updated": utc_now(),
            "configurations": {},
        }
    for config in all_configurations(campaign):
        name = config["name"]
        if only and name != only:
            continue
        directory = output_directory(campaign) / "configurations" / name
        merged: dict[str, dict[str, float]] = {}
        ordinary_merged: dict[str, dict[str, float]] = {}
        merged_cells: dict[str, dict[str, float]] = {}
        total_events = 0
        ordinary_batches = 0
        batches = []
        proposal_history = []
        ordinary_performance: dict[str, dict[str, Any]] = {}
        biased_performance: dict[str, dict[str, Any]] = {}
        for record_path in sorted(directory.glob("batch_*.json")) if directory.exists() else []:
            record = json.loads(record_path.read_text())
            if not record.get("valid"):
                continue
            root_file = Path(record["root_file"])
            if record.get("analysis", {}).get("version") == ANALYSIS_VERSION:
                events = int(record["analysis"]["events"])
                stats = record["analysis"]["statistics"]
                cells = record["analysis"].get("cells")
            else:
                events, stats, cells = analyze_root(root_file, config)
                record["analysis"] = {
                    "version": ANALYSIS_VERSION, "events": events,
                    "statistics": stats, "cells": cells,
                }
                atomic_json(record_path, record)
            if cells is None:
                events, stats, cells = analyze_root(root_file, config)
                record["analysis"] = {
                    "version": ANALYSIS_VERSION, "events": events,
                    "statistics": stats, "cells": cells,
                }
                atomic_json(record_path, record)
            if events != record["events_requested"]:
                raise RuntimeError(f"event mismatch in {root_file}")
            total_events += events
            batches.append(record_path.name)
            proposal = record.get("proposal", {"kind": "ordinary"})
            elapsed = float(record.get("elapsed_seconds", 0.0))
            proposal_history.append({
                "batch": record_path.name,
                "proposal": proposal,
                "events": events,
                "elapsed_seconds": elapsed,
            })
            if proposal.get("kind", "ordinary") == "ordinary":
                ordinary_batches += 1
            for key, value in stats.items():
                target = merged.setdefault(key, new_merged_stat())
                target["signal_events"] += value["signal_events"]
                target["sum"] += value["sum"]
                target["sum2"] += value["sum2"]
                target["sum3"] += value["sum3"]
                target["sum4"] += value["sum4"]
                target["maximum"] = max(target["maximum"], value["maximum"])
                target["largest"] = merge_largest(target["largest"], value["largest"])
                variance_numerator = events * batch_variance(events, value)
                target["variance_numerator"] += variance_numerator
                target["maximum_batch_sum"] = max(
                    target["maximum_batch_sum"], abs(value["sum"])
                )
                target["history"].append({
                    "batch": record_path.name,
                    "events": events,
                    "sum": value["sum"],
                    "sum2": value["sum2"],
                    "sum3": value["sum3"],
                    "sum4": value["sum4"],
                    "variance_numerator": variance_numerator,
                    "elapsed_seconds": elapsed,
                    "proposal": proposal,
                })
                if proposal.get("kind", "ordinary") == "ordinary":
                    ordinary = ordinary_merged.setdefault(key, new_merged_stat())
                    ordinary["signal_events"] += value["signal_events"]
                    ordinary["sum"] += value["sum"]
                    ordinary["sum2"] += value["sum2"]
                    ordinary["sum3"] += value["sum3"]
                    ordinary["sum4"] += value["sum4"]
                    ordinary["maximum"] = max(ordinary["maximum"], value["maximum"])
                    ordinary["largest"] = merge_largest(
                        ordinary["largest"], value["largest"]
                    )
                    ordinary["history"].append({
                        "batch": record_path.name, "events": events,
                        "sum": value["sum"], "sum2": value["sum2"],
                        "sum3": value["sum3"], "sum4": value["sum4"],
                        "elapsed_seconds": elapsed,
                    })
                    add_performance_sample(
                        ordinary_performance.setdefault(key, {}), events,
                        variance_numerator, elapsed, value["signal_events"],
                    )
                elif key == proposal.get("target_statistic"):
                    signature = proposal_signature(proposal)
                    performance = biased_performance.setdefault(signature, {
                        "axis": proposal["axis"],
                        "cell": proposal["cell"],
                        "target_statistic": key,
                        "minimum": proposal["minimum"],
                        "maximum": proposal["maximum"],
                    })
                    add_performance_sample(
                        performance, events, variance_numerator, elapsed,
                        value["signal_events"],
                    )
            for key, value in cells.items():
                target = merged_cells.setdefault(key, {
                    "minimum": value["minimum"],
                    "maximum_bound": value["maximum_bound"],
                    "sample_events": 0, "signal_events": 0, "sum": 0.0,
                    "sum2": 0.0, "maximum": 0.0,
                })
                target["sample_events"] += value["sample_events"]
                target["signal_events"] += value["signal_events"]
                target["sum"] += value["sum"]
                target["sum2"] += value["sum2"]
                target["maximum"] = max(target["maximum"], value["maximum"])
        for key, value in merged.items():
            value.update(statistic_result(total_events, value))
            value["reliability"] = reliability_result(
                campaign, ordinary_merged.get(key, new_merged_stat())
            )
        proposal_performance = {}
        minimum_trials = int(campaign["adaptive"]["minimum_trials_before_reject"])
        rejection_ratio = float(campaign["adaptive"]["rejection_cost_ratio"])
        for signature, value in biased_performance.items():
            target_key = value["target_statistic"]
            baseline = ordinary_performance.get(target_key, {})
            baseline_cost = variance_cost(baseline)
            proposal_cost = variance_cost(value)
            relative_efficiency = (
                baseline_cost / proposal_cost
                if baseline_cost is not None and proposal_cost not in (None, 0.0)
                else None
            )
            enough_trials = value["batches"] >= minimum_trials
            no_response = enough_trials and value["signal_events"] == 0
            too_expensive = (
                enough_trials and baseline_cost not in (None, 0.0)
                and proposal_cost is not None
                and proposal_cost >= rejection_ratio * baseline_cost
            )
            value.update({
                "ordinary_variance_cost": baseline_cost,
                "proposal_variance_cost": proposal_cost,
                "relative_efficiency": relative_efficiency,
                "rejected": no_response or too_expensive,
                "rejection_reason": (
                    "no target responses in repeated trials" if no_response
                    else "variance per runtime is worse than ordinary sampling"
                    if too_expensive else None
                ),
            })
            proposal_performance[signature] = value
        result = {
            "events": total_events,
            "batches": batches,
            "ordinary_batches": ordinary_batches,
            "proposal_history": proposal_history,
            "proposal_performance": proposal_performance,
            "statistics": merged,
            "cells": merged_cells,
        }
        apply_completion_states(campaign, result)
        summary["configurations"][name] = result
    atomic_json(summary_path, summary)
    return summary


def choose_proposal(campaign: dict[str, Any], name: str,
                    force_ordinary: bool = False) -> dict[str, Any]:
    ordinary = {"kind": "ordinary", "reason": "physical generator distribution"}
    if force_ordinary or not campaign.get("adaptive", {}).get("enabled", False):
        return ordinary
    summary_path = output_directory(campaign) / "summary.json"
    if not summary_path.is_file():
        return ordinary
    result = json.loads(summary_path.read_text()).get("configurations", {}).get(name)
    if not result or result.get("ordinary_batches", 0) < int(
            campaign["precision"]["minimum_ordinary_batches"]):
        return ordinary
    if result.get("state") in ("measured", "statistics_limited"):
        return ordinary

    control_interval = int(campaign["adaptive"]["ordinary_control_interval"])
    trailing_biased = 0
    for item in reversed(result.get("proposal_history", [])):
        if item["proposal"].get("kind", "ordinary") == "ordinary":
            break
        trailing_biased += 1
    if trailing_biased >= control_interval - 1:
        return {
            "kind": "ordinary",
            "reason": f"periodic ordinary control after {trailing_biased} biased batches",
        }

    candidates = []
    empty_candidates = []
    for key, value in result.get("statistics", {}).items():
        group, ancestry, observable = key.split("/")
        if ancestry != "all" or group not in detector_groups()[0]:
            continue
        if value.get("state") in ("measured", "negligible", "not_relevant"):
            continue
        if value.get("estimate", 0.0) == 0.0:
            empty_candidates.append((0.0, key, group, observable))
            continue
        else:
            score = max(
                value.get("relative_error") or 0.0,
                value.get("largest_event_fraction") or 0.0,
            )
        candidates.append((score, key, group, observable))
    candidate_statistics = candidates if candidates else empty_candidates
    if not candidate_statistics:
        return ordinary

    config = config_by_name(campaign, name)
    enabled_axes = set(campaign["adaptive"].get("axes", []))
    enabled_axes &= set(config["interaction_settings"].get("adaptive_axes", enabled_axes))
    rejected = {
        signature for signature, value in result.get("proposal_performance", {}).items()
        if value.get("rejected")
    }
    for _, statistic_key, group, observable in sorted(
            candidate_statistics, key=lambda item: item[0], reverse=True):
        cells = []
        for key, value in result.get("cells", {}).items():
            axis, index, cell_group, cell_observable = key.split("/")
            if axis not in enabled_axes or cell_group != group \
                    or cell_observable != observable:
                continue
            cells.append((axis, int(index), value))
        if result["statistics"][statistic_key].get("estimate", 0.0) == 0.0:
            cells.sort(key=lambda item: (item[2]["sample_events"], item[0], item[1]))
            reason = f"least-sampled cell for empty {statistic_key} tail"
        else:
            cells.sort(
                key=lambda item: item[2]["sum2"] / max(1, item[2]["sample_events"]),
                reverse=True,
            )
            reason = f"largest conditional second moment for {statistic_key}"
        for axis, index, cell in cells:
            proposal = {
                "kind": "mixture",
                "axis": axis,
                "cell": index,
                "minimum": cell["minimum"],
                "maximum": cell["maximum_bound"],
                "physical_fraction": float(campaign["adaptive"]["physical_fraction"]),
                "target_statistic": statistic_key,
                "reason": reason,
            }
            if proposal_signature(proposal) not in rejected:
                return proposal
    return {
        "kind": "ordinary",
        "reason": "all candidate biased proposals performed worse than ordinary sampling",
    }


def run_one(campaign_path: Path, name: str, dry_run: bool, events_override: int | None,
            force_ordinary: bool = False,
            proposal_override: dict[str, Any] | None = None) -> int:
    errors = check_campaign(campaign_path)
    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        return 2
    campaign = load_toml(campaign_path)
    config = config_by_name(campaign, name)
    proposal = proposal_override or choose_proposal(campaign, name, force_ordinary)
    events = events_override or int(campaign["batch_events"])
    threads = int(campaign["threads"])
    root_dir = output_directory(campaign)
    config_dir = root_dir / "configurations" / name
    number = next_batch_number(config_dir)
    stem = f"batch_{number:06d}"
    root_file = config_dir / f"{stem}.root"
    macro_file = config_dir / f"{stem}.mac"
    record_file = config_dir / f"{stem}.json"
    log_file = config_dir / f"{stem}.log"
    seed = random_seed()
    source_macro = resolve_repo_path(campaign["base_macro"])
    macro_text = generated_macro_text(source_macro, config, proposal, root_file, seed, events)
    command = [str(REPO / "build/remoll"), "-t", str(threads), "-m", str(macro_file)]

    if dry_run:
        print(f"configuration: {name}")
        print(f"events:        {events}")
        print(f"threads:       {threads}")
        print(f"source macro:  {source_macro}")
        print(f"proposal:      {proposal['kind']} ({proposal['reason']})")
        if proposal["kind"] == "mixture":
            print(f"target cell:   {proposal['axis']} [{proposal['minimum']:.6g}, "
                  f"{proposal['maximum']:.6g}]")
        print(f"output ROOT:   {root_file}")
        print("command:       " + " ".join(command))
        print(f"precision detector IDs: {len(detector_groups()[1])}")
        print(f"diagnostic detector IDs: {len(detector_groups()[2])}")
        return 0

    config_dir.mkdir(parents=True, exist_ok=True)
    if disk_used_gb(root_dir) >= float(campaign["max_disk_gb"]):
        print("error: campaign disk limit reached", file=sys.stderr)
        return 1
    if not campaign.get("run_on_battery", False) and not on_ac_power():
        print("error: campaign is configured to run only on AC power", file=sys.stderr)
        return 1

    macro_file.write_text(macro_text)
    state = read_state(campaign)
    state["active"] = {"configuration": name, "batch": number, "events": events,
                       "started": utc_now(), "pid": None}
    state["control"] = "running"
    write_state(campaign, state)
    started = time.monotonic()
    started_iso = utc_now()
    print(f"starting {name} batch {number}: {events} events, {threads} threads")
    with log_file.open("ab", buffering=0) as log:
        process = subprocess.Popen(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
        state["active"]["pid"] = process.pid
        write_state(campaign, state)
        run_alerts: list[str] = []
        while process.poll() is None:
            time.sleep(max(1, int(campaign["status_update_seconds"])))
            for alert in extract_alerts(log_file):
                if alert not in run_alerts and len(run_alerts) < 100:
                    run_alerts.append(alert)
            cap_file(log_file, ACTIVE_LOG_LIMIT)
            state = read_state(campaign)
            state["active"].update({
                "elapsed_seconds": round(time.monotonic() - started, 1),
                "root_bytes": root_file.stat().st_size if root_file.exists() else 0,
            })
            write_state(campaign, state)
            if state.get("control") == "stop_requested":
                process.send_signal(signal.SIGTERM)

    elapsed = time.monotonic() - started
    status = process.returncode
    entries, required_branches = root_integrity(root_file) if root_file.is_file() else (None, False)
    valid = status == 0 and entries == events and required_branches
    for alert in extract_alerts(log_file):
        if alert not in run_alerts and len(run_alerts) < 100:
            run_alerts.append(alert)
    alerts = run_alerts
    record = {
        "version": 1,
        "configuration": name,
        "batch": number,
        "source_macro": str(source_macro),
        "source_macro_sha256": sha256(source_macro),
        "generated_macro": str(macro_file),
        "generated_macro_sha256": sha256(macro_file),
        "root_file": str(root_file),
        "root_bytes": root_file.stat().st_size if root_file.exists() else 0,
        "root_entries": entries,
        "required_branches_present": required_branches,
        "seed": seed,
        "events_requested": events,
        "threads": threads,
        "proposal": proposal,
        "command": command,
        "started": started_iso,
        "finished": utc_now(),
        "elapsed_seconds": round(elapsed, 3),
        "exit_status": status,
        "valid": valid,
        "alerts": alerts,
    }
    atomic_json(record_file, record)
    state = read_state(campaign)
    state["active"] = None
    if state.get("control") == "stop_requested":
        state["control"] = "stopped"
    elif state.get("control") != "paused":
        state["control"] = "ready"
    state["last_batch"] = record
    write_state(campaign, state)
    if valid:
        log_file.unlink(missing_ok=True)
        rebuild_summary(campaign, name)
        print(f"completed {name} batch {number}: {entries} events in {elapsed:.1f}s")
        if alerts:
            print(f"alerts retained in batch record: {len(alerts)}")
        return 0
    cap_file(log_file, FAILED_LOG_LIMIT)
    print(f"failed {name} batch {number}: exit={status}, entries={entries}", file=sys.stderr)
    print(f"log: {log_file}", file=sys.stderr)
    return 1


def set_control(campaign_path: Path, value: str) -> int:
    campaign = load_toml(campaign_path)
    state = read_state(campaign)
    state["control"] = value
    write_state(campaign, state)
    print(value)
    return 0


def print_status(campaign_path: Path) -> int:
    campaign = load_toml(campaign_path)
    state = read_state(campaign)
    print(f"campaign: {campaign['title']}")
    print(f"state:    {state.get('control', 'unknown')}")
    active = state.get("active")
    if active:
        print(f"active:   {active['configuration']} batch {active['batch']}")
        print(f"events:   {active['events']}")
        print(f"elapsed:  {active.get('elapsed_seconds', 0)} s")
        print(f"ROOT:     {active.get('root_bytes', 0) / 1e6:.1f} MB")
    configs = all_configurations(campaign)
    summary_path = output_directory(campaign) / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {"configurations": {}}
    state_counts: dict[str, int] = {}
    for config in configs:
        value = summary.get("configurations", {}).get(config["name"], {})
        label = value.get("state", "pending")
        state_counts[label] = state_counts.get(label, 0) + 1
    print("inventory: " + ", ".join(f"{key}={value}" for key, value in sorted(state_counts.items())))
    print("\nconfiguration                         state               batches  events     group RSE  tile RSE")
    pending_shown = False
    for config in sorted(configs, key=lambda item: (item.get("priority", 999), item["name"])):
        directory = output_directory(campaign) / "configurations" / config["name"]
        count = sum(1 for path in directory.glob("batch_*.json")
                    if json.loads(path.read_text()).get("valid")) if directory.exists() else 0
        result = summary.get("configurations", {}).get(config["name"], {})
        config_state = result.get("state", "pending")
        if config_state == "pending":
            if pending_shown:
                continue
            pending_shown = True
            config_state = "next"
        events = result.get("events", 0)
        group_errors = [
            value["relative_error"]
            for key, value in result.get("statistics", {}).items()
            if key.endswith("/all/rate") and not key.startswith("tile_")
            and key.split("/", 1)[0] in detector_groups()[0]
            and value.get("estimate") and value.get("relative_error") is not None
        ]
        tile_errors = [
            value["relative_error"]
            for key, value in result.get("statistics", {}).items()
            if key.startswith("tile_") and key.endswith("/all/rate")
            and value.get("estimate") and value.get("relative_error") is not None
        ]
        group_worst = f"{max(group_errors) * 100:.2f}%" if group_errors else "--"
        tile_worst = f"{max(tile_errors) * 100:.2f}%" if tile_errors else "--"
        print(f"{config['name']:<37} {config_state:<19} {count:<7}  {events:<9}  "
              f"{group_worst:<9}  {tile_worst}")
    return 0


def inspect_configuration(campaign_path: Path, name: str, show_tiles: bool = False) -> int:
    campaign = load_toml(campaign_path)
    config_by_name(campaign, name)
    summary_path = output_directory(campaign) / "summary.json"
    if not summary_path.is_file():
        print("no analyzed batches")
        print("next proposal: ordinary physical-distribution pilot")
        return 0
    result = json.loads(summary_path.read_text()).get("configurations", {}).get(name)
    if not result or not result.get("events"):
        print("no analyzed batches")
        print("next proposal: ordinary physical-distribution pilot")
        return 0
    print(f"configuration: {name}")
    print(f"state:         {result['state']}")
    print(f"events:        {result['events']}")
    print(f"batches:       {len(result['batches'])}")
    proposal = choose_proposal(campaign, name)
    print(f"next proposal: {proposal['kind']} - {proposal['reason']}")
    if proposal["kind"] == "mixture":
        print(f"target:        {proposal['axis']} cell {proposal['cell']} "
              f"[{proposal['minimum']:.6g}, {proposal['maximum']:.6g}]")
    performance = result.get("proposal_performance", {})
    if performance:
        rejected = sum(bool(value.get("rejected")) for value in performance.values())
        print(f"bias trials:   {len(performance)} proposal choices evaluated, "
              f"{rejected} rejected")
    print("\ndetector group             rate estimate      rate RSE  E*rate RSE  secondary  state")
    groups = detector_groups()[0]
    for group in groups:
        value = result["statistics"].get(f"{group}/all/rate", {})
        energy_value = result["statistics"].get(f"{group}/all/rate_energy", {})
        secondary_value = result["statistics"].get(f"{group}/secondary/rate", {})
        estimate = value.get("estimate", 0.0)
        relative = value.get("relative_error")
        state = value.get("state", "running")
        rse = f"{relative * 100:6.2f}%" if relative is not None else "    -- "
        energy_relative = (
            energy_value.get("relative_error") if group.startswith("showermax_") else None
        )
        energy_rse = f"{energy_relative * 100:6.2f}%" if energy_relative is not None else "    -- "
        secondary = secondary_value.get("estimate", 0.0) / estimate if estimate else None
        secondary_text = f"{secondary * 100:6.1f}%" if secondary is not None else "    -- "
        if group.startswith("showermax_") and energy_value.get("state") \
                not in (None, "measured", "negligible"):
            state = energy_value["state"]
        print(f"{group:<26} {estimate:13.6e}  {rse}   {energy_rse}   "
              f"{secondary_text}   {state}")
    print("\nordinary-control reliability (rate tally)")
    print("detector group             histories   VOV      R slope  tail     next max  failed")
    for group in groups:
        reliability = result["statistics"].get(
            f"{group}/all/rate", {}
        ).get("reliability", {})
        vov = reliability.get("variance_of_variance")
        slope = reliability.get("rse_log_slope")
        tail = reliability.get("tail_index")
        next_fraction = reliability.get("next_largest_history_fraction")
        print(
            f"{group:<26} {reliability.get('ordinary_signal_events', 0):>9}  "
            f"{vov if vov is not None else float('nan'):7.3g}  "
            f"{slope if slope is not None else float('nan'):7.3f}  "
            f"{tail if tail is not None else float('nan'):7.3g}  "
            f"{next_fraction * 100 if next_fraction is not None else float('nan'):7.3f}%  "
            f"{len(reliability.get('failed_checks', []))}"
        )
    tile_statistics = [
        (key, value) for key, value in result["statistics"].items()
        if key.startswith("tile_")
    ]
    state_counts: dict[str, int] = {}
    for _, value in tile_statistics:
        state = value.get("state", "reported_only")
        state_counts[state] = state_counts.get(state, 0) + 1
    print("\nindividual tile tallies: " + ", ".join(
        f"{state}={count}" for state, count in sorted(state_counts.items())
    ))
    if show_tiles:
        print("tile tally                         estimate       RSE     VOV      tail    state")
        for key, value in sorted(tile_statistics):
            reliability = value.get("reliability", {})
            relative = value.get("relative_error")
            vov = reliability.get("variance_of_variance")
            tail = reliability.get("tail_index")
            print(
                f"{key:<34} {value.get('estimate', 0.0):11.4e}  "
                f"{relative * 100 if relative is not None else float('nan'):6.2f}%  "
                f"{vov if vov is not None else float('nan'):7.3g}  "
                f"{tail if tail is not None else float('nan'):7.3g}  "
                f"{value.get('state', 'reported_only')}"
            )
    return 0


def print_inventory(campaign_path: Path) -> int:
    campaign = load_toml(campaign_path)
    configs = sorted(all_configurations(campaign),
                     key=lambda item: (item.get("priority", 999), item["name"]))
    print(f"configurations: {len(configs)}")
    print("name                                      target  sieve  interaction   energy  current")
    for config in configs:
        print(f"{config['name']:<41} {config['target']:<6}  {config['sieve']:<5}  "
              f"{config['interaction']:<12}  {config['energy_gev']:>4g} GeV  "
              f"{config['current_uA']:>4g} uA")
    return 0


def start(campaign_path: Path) -> int:
    errors = check_campaign(campaign_path, require_enabled=True)
    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        return 2
    campaign = load_toml(campaign_path)
    enabled = sorted(
        (config for config in all_configurations(campaign) if config.get("enabled")),
        key=lambda item: (item.get("priority", 999), item["name"]),
    )
    rebuild_summary(campaign)
    for config in enabled:
        consecutive_failures = 0
        while True:
            state = read_state(campaign)
            if state.get("control") in ("paused", "stop_requested", "stopped"):
                return 0
            summary = rebuild_summary(campaign, config["name"])
            result = summary["configurations"][config["name"]]
            if result["state"] in ("measured", "statistics_limited"):
                break
            batch_result = run_one(campaign_path, config["name"], False, None)
            if batch_result == 0:
                consecutive_failures = 0
            else:
                consecutive_failures += 1
                if consecutive_failures >= 3:
                    print(f"configuration failed three consecutive batches: {config['name']}",
                          file=sys.stderr)
                    break
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="campaign")
    parser.add_argument("--campaign", type=Path, default=DEFAULT_CAMPAIGN)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("check")
    subparsers.add_parser("status")
    subparsers.add_parser("inventory")
    inspect_parser = subparsers.add_parser("inspect")
    inspect_parser.add_argument("configuration")
    inspect_parser.add_argument("--tiles", action="store_true")
    run_parser = subparsers.add_parser("run-one")
    run_parser.add_argument("configuration")
    run_parser.add_argument("--events", type=int)
    run_parser.add_argument("--dry-run", action="store_true")
    run_parser.add_argument("--ordinary", action="store_true")
    subparsers.add_parser("start")
    subparsers.add_parser("report")
    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("configuration", nargs="?")
    subparsers.add_parser("pause")
    subparsers.add_parser("resume")
    subparsers.add_parser("stop")
    args = parser.parse_args()
    campaign_path = args.campaign.resolve()

    if args.command == "check":
        errors = check_campaign(campaign_path)
        if errors:
            for error in errors:
                print(f"error: {error}")
            return 2
        groups, precision, diagnostic = detector_groups()
        print("campaign is valid")
        print(f"precision groups: {', '.join(groups)}")
        print(f"precision detector IDs: {len(precision)}")
        print(f"diagnostic detector IDs: {len(diagnostic)}")
        return 0
    if args.command == "status":
        return print_status(campaign_path)
    if args.command == "inventory":
        return print_inventory(campaign_path)
    if args.command == "inspect":
        return inspect_configuration(campaign_path, args.configuration, args.tiles)
    if args.command == "run-one":
        return run_one(campaign_path, args.configuration, args.dry_run, args.events, args.ordinary)
    if args.command == "start":
        return start(campaign_path)
    if args.command == "report":
        campaign = load_toml(campaign_path)
        print(output_directory(campaign) / "summary.json")
        return 0
    if args.command == "analyze":
        campaign = load_toml(campaign_path)
        rebuild_summary(campaign, args.configuration)
        print(output_directory(campaign) / "summary.json")
        return 0
    if args.command == "pause":
        return set_control(campaign_path, "paused")
    if args.command == "resume":
        return set_control(campaign_path, "ready")
    if args.command == "stop":
        return set_control(campaign_path, "stop_requested")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
