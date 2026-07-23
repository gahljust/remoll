#!/usr/bin/env bash
#
# Produce one compact table of target-level kinematic windows whose electrons
# survive COL1/COL2 into one selected septant. The GPU transport executable
# generates the rate nodes internally and reduces accepted nodes in memory:
# no rate ROOT, transport ROOT, hit tree, track tree, or per-node TSV is made.
#
# The q01--q99 intervals are the suggested targeted components for remoll's
# full-support mixture biasers:
#   theta       -> /remoll/bias/thcom/{min,max}
#   phi         -> /remoll/bias/phi/{min,max}
#   beamp       -> /remoll/bias/beamp/{min,max}
#   vertexz     -> /remoll/bias/vertexz/{minFraction,maxFraction}
#   outgoinge   -> /remoll/bias/outgoinge/{minFraction,maxFraction}
#
# Usage:
#   run_col2_kinematics.sh [--quick|--standard|--production]
#       [--septant 0..6] [--only REGEX] [--out DIR]
#       [--force] [--dry-run] [--list] [--build]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
GPU_BIN="${REPO_ROOT}/_local_additions_archive/analytic_solver/gpu/transport_gpu"
GPU_BUILD="${REPO_ROOT}/_local_additions_archive/analytic_solver/gpu/build_gpu.sh"
MASK_DIR="${REPO_ROOT}/_local_additions_archive/analytic_solver/transport"

MODE="production"
OUTBASE="${REPO_ROOT}/_local_additions_archive/analytic_solver/pipeline/col2_kinematics"
ONLY=""
SEPTANT=0
FORCE=0
DRY=0
DO_BUILD=0
LIST=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --quick) MODE="quick"; shift ;;
    --standard) MODE="standard"; shift ;;
    --production) MODE="production"; shift ;;
    --septant) SEPTANT="$2"; shift 2 ;;
    --only) ONLY="$2"; shift 2 ;;
    --out) OUTBASE="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --dry-run) DRY=1; shift ;;
    --build) DO_BUILD=1; shift ;;
    --list) LIST=1; shift ;;
    -h|--help)
      sed -n '2,22p' "$0" | sed 's/^# //'
      exit 0
      ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

[[ "${SEPTANT}" =~ ^[0-6]$ ]] || {
  echo "--septant must be an integer from 0 through 6" >&2
  exit 2
}

case "${MODE}" in
  quick)      VERTEX_NODES=8;  ENERGY_NODES=16; THETA_NODES=48; PHI_NODES=12 ;;
  standard)   VERTEX_NODES=12; ENERGY_NODES=24; THETA_NODES=64; PHI_NODES=24 ;;
  production) VERTEX_NODES=24; ENERGY_NODES=48; THETA_NODES=96; PHI_NODES=40 ;;
  *) echo "unknown mode: ${MODE}" >&2; exit 2 ;;
esac
PHI_FOLD=7

if [[ "${DO_BUILD}" == "1" ]]; then
  "${GPU_BUILD}"
fi
[[ -x "${GPU_BIN}" ]] || {
  echo "missing GPU transport: ${GPU_BIN} (run with --build)" >&2
  exit 1
}
[[ -f "${MASK_DIR}/COL1.mask" && -f "${MASK_DIR}/COL2.mask" ]] || {
  echo "COL1/COL2 masks are missing from ${MASK_DIR}" >&2
  exit 1
}

FIELD_ARGS=(
  --field "${REPO_ROOT}/map_directory/V2U.1a.50cm.parallel.txt"
  --field "${REPO_ROOT}/map_directory/DS_TM1-4_CoilA-G_ll_TM2-4_out3mm.txt"
)

target_args() {
  TARGET_ARGS=()
  case "$1" in
    lh2)
      TARGET_ARGS=(--target-Z 1 --target-A 1.008 --thickness-mm 1250
        --density-g-cm3 0.0715 --molar-mass 1.008
        --radiation-length-mm 8904 --target-z0-mm -5125
        --target-z-span-mm 1250)
      ;;
    c12_us)
      TARGET_ARGS=(--target-Z 6 --target-A 12 --thickness-mm 0.254
        --density-g-cm3 2.0 --molar-mass 12.0107
        --radiation-length-mm 213.5 --target-z0-mm -5124.627
        --target-z-span-mm 0.254)
      ;;
    c12_ms)
      TARGET_ARGS=(--target-Z 6 --target-A 12 --thickness-mm 0.254
        --density-g-cm3 2.0 --molar-mass 12.0107
        --radiation-length-mm 213.5 --target-z0-mm -4500.000
        --target-z-span-mm 0.254)
      ;;
    c12_ds)
      TARGET_ARGS=(--target-Z 6 --target-A 12 --thickness-mm 0.254
        --density-g-cm3 2.0 --molar-mass 12.0107
        --radiation-length-mm 213.5 --target-z0-mm -3875.627
        --target-z-span-mm 0.254)
      ;;
    *) echo "unknown target $1" >&2; return 1 ;;
  esac
}

channel_args() {
  CHANNEL_ARGS=()
  case "$1" in
    moller)
      CHANNEL_ARGS=(--theta-min-deg 30 --theta-max-deg 150)
      ;;
    ep_elastic)
      local floor
      floor="$(python3 -c "print(80.0*${2}/11000.0)")"
      CHANNEL_ARGS=(--theta-min-deg 0.1 --theta-max-deg 3.0
        --ep-internal-brems --eout-floor-mev "${floor}")
      ;;
    ep_inelastic|c12_elastic|c12_inelastic)
      CHANNEL_ARGS=(--theta-min-deg 0.1 --theta-max-deg 5.0)
      ;;
    *) echo "unknown channel $1" >&2; return 1 ;;
  esac
}

cells=("lh2:11000:moller" "lh2:11000:ep_elastic" "lh2:11000:ep_inelastic")
for target in c12_us c12_ms c12_ds; do
  for energy in 2200 4400 6600 8800 11000; do
    for channel in moller c12_elastic c12_inelastic; do
      cells+=("${target}:${energy}:${channel}")
    done
  done
done

if [[ -n "${ONLY}" ]]; then
  filtered=()
  for cell in "${cells[@]}"; do
    display_cell="${cell//:/_}"
    if [[ "${cell}" =~ ${ONLY} ]] || [[ "${display_cell}" =~ ${ONLY} ]]; then
      filtered+=("${cell}")
    fi
  done
  cells=("${filtered[@]}")
fi

if [[ "${LIST}" == "1" ]]; then
  printf '%s\n' "${cells[@]}"
  echo "total=${#cells[@]} mode=${MODE} septant=${SEPTANT}"
  exit 0
fi
[[ "${#cells[@]}" -gt 0 ]] || { echo "no configurations matched --only" >&2; exit 2; }

mkdir -p "${OUTBASE}"
TABLE="${OUTBASE}/col2_bias_kinematics.tsv"
if [[ "${FORCE}" == "1" ]]; then
  : > "${TABLE}"
fi

PLANES="$(mktemp "${TMPDIR:-/tmp}/remoll-col2-planes.XXXXXX")"
cleanup() { rm -f "${PLANES}" "${TEMP_SUMMARY:-}" "${TEMP_LOG:-}"; }
trap cleanup EXIT
cat > "${PLANES}" <<EOF
name,z_m,r_min_m,r_max_m,x_min_m,x_max_m,y_min_m,y_max_m,is_aperture,n_septant,sept_phi0_deg,sept_half_deg
COL1,0.600,0,0.20,-10,10,-10,10,1,1,0,180,${MASK_DIR}/COL1.mask
COL2,0.890,0,0.20,-10,10,-10,10,1,1,0,180,${MASK_DIR}/COL2.mask
COL2_exit,0.900,0,10,-10,10,-10,10,0,1,0,180
EOF

append_summary() {
  local cell="$1" target="$2" energy="$3" channel="$4" input="$5"
  python3 - "${TABLE}" "${cell}" "${target}" "${energy}" "${channel}" \
    "${SEPTANT}" "${MODE}" "${PHI_NODES}" "${PHI_FOLD}" "${input}" <<'PY'
import csv
import math
import sys

table, cell, target, energy, channel, septant, mode, phi_nodes, phi_fold, source = sys.argv[1:]
with open(source, newline="") as stream:
    row = next(csv.DictReader(stream, delimiter="\t"))

# The solve uses one symmetry sector. Rotate its generator-phi window to the
# requested physical septant; no per-node sevenfold expansion is necessary.
shift = int(septant) * 360.0 / int(phi_fold)
for key in ("phi_q01_deg", "phi_q50_deg", "phi_q99_deg"):
    row[key] = f"{(float(row[key]) + shift) % 360.0:.12g}"

metadata = {
    "cell": cell, "target": target, "energy_mev": energy, "channel": channel,
    "septant": septant, "resolution": mode, "phi_nodes": phi_nodes,
    "phi_fold": phi_fold,
}
fields = list(metadata) + [key for key in row if key != "plane"]
write_header = not __import__("pathlib").Path(table).is_file() or \
               __import__("pathlib").Path(table).stat().st_size == 0
with open(table, "a", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                            lineterminator="\n")
    if write_header:
        writer.writeheader()
    writer.writerow(metadata | {key: row[key] for key in row if key != "plane"})
PY
}

idx=0
for cell in "${cells[@]}"; do
  idx=$((idx + 1))
  IFS=: read -r target energy channel <<< "${cell}"
  tag="${target}_${energy}_${channel}"

  if [[ -s "${TABLE}" && "${FORCE}" == "0" ]] && \
     awk -F '\t' -v key="${cell}" 'NR>1 && $1==key {found=1} END{exit !found}' "${TABLE}"; then
    echo "[${idx}/${#cells[@]}] ${tag}: already in table, skipping"
    continue
  fi

  target_args "${target}"
  channel_args "${channel}" "${energy}"
  TEMP_SUMMARY="$(mktemp "${TMPDIR:-/tmp}/remoll-col2-summary.XXXXXX")"
  TEMP_LOG="$(mktemp "${TMPDIR:-/tmp}/remoll-col2-run.XXXXXX")"
  gpu_cmd=("${GPU_BIN}" --channel "${channel}" --beam-mev "${energy}" --current-uA 1
    "${TARGET_ARGS[@]}" "${CHANNEL_ARGS[@]}"
    --screening --vertex-nodes "${VERTEX_NODES}" --energy-nodes "${ENERGY_NODES}"
    --theta-nodes "${THETA_NODES}" --phi-nodes "${PHI_NODES}" --phi-fold "${PHI_FOLD}"
    "${FIELD_ARGS[@]}" --planes "${PLANES}" --z-stop-m 1.0 --r-stop-m 3
    --ds-field-m 0.01 --ds-drift-m 0.25
    --kinematics-summary COL2_exit --kinematics-output "${TEMP_SUMMARY}"
    --no-json --output "${OUTBASE}/.unused_${tag}")

  echo "[${idx}/${#cells[@]}] ${tag}"
  if [[ "${DRY}" == "1" ]]; then
    printf '  +'; printf ' %q' "${gpu_cmd[@]}"; printf '\n'
    rm -f "${TEMP_SUMMARY}" "${TEMP_LOG}"
    TEMP_SUMMARY=""; TEMP_LOG=""
    continue
  fi

  if "${gpu_cmd[@]}" > "${TEMP_LOG}" 2>&1; then
    append_summary "${cell}" "${target}" "${energy}" "${channel}" "${TEMP_SUMMARY}"
    echo "  appended ${TABLE}"
    rm -f "${TEMP_SUMMARY}" "${TEMP_LOG}"
    TEMP_SUMMARY=""; TEMP_LOG=""
  else
    echo "FAILED: ${tag}; last log lines:" >&2
    tail -40 "${TEMP_LOG}" >&2
    exit 1
  fi
done

if [[ "${DRY}" == "0" ]]; then
  echo "DONE: ${TABLE}"
  echo "Only the compact table is retained; temporary planes, logs, and summaries were removed."
fi
