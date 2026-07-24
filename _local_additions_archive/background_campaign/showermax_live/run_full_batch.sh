#!/usr/bin/env bash
set -euo pipefail

# Restartable full detector-support campaign. Each configuration owns its
# normal showermax_live run directory, so completed histories are reused. A new
# configuration is produced as one 150k ROOT file; a partial configuration gets
# one ROOT file containing only the exact remainder to 150k.

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
runner="${repo_dir}/_local_additions_archive/background_campaign/showermax_live/showermax_live.py"
event_cap="${EVENT_CAP:-150000}"
batch_events="${BATCH_EVENTS:-${event_cap}}"
failures=()

run_cell() {
  local cell="$1"
  echo
  echo "================================================================"
  echo "Campaign: ${cell}  sieve=out  cap=${event_cap}"
  echo "================================================================"
  if ! python3 "${runner}" \
    --cell "${cell}" \
    --sieve out \
    --batch-events "${batch_events}" \
    --max-events "${event_cap}" \
    --no-dashboard \
    --no-browser \
    start; then
    failures+=("${cell}")
    echo "FAILED: ${cell}; recorded for retry. Continuing to the next configuration." >&2
  fi
}

# Validate the complete 14-cell queue before launching the first long job.
queue=(
  "lh2:11000:ep_elastic"
  "lh2:11000:ep_inelastic"
)
for target in c12_us c12_ms c12_ds; do
  for energy in 2200 4400; do
    queue+=("${target}:${energy}:moller")
    queue+=("${target}:${energy}:c12_elastic")
  done
done

echo "Preflighting ${#queue[@]} campaign configurations..."
for cell in "${queue[@]}"; do
  python3 "${runner}" \
    --cell "${cell}" \
    --sieve out \
    --batch-events "${batch_events}" \
    --max-events "${event_cap}" \
    --no-dashboard \
    --no-browser \
    check
done
echo "All campaign configurations passed preflight."

# LH2 Møller is already complete in its existing campaign; do not create a
# second sieve-out Møller run. Finish the other two LH2 interactions.
run_cell "${queue[0]}"
run_cell "${queue[1]}"

# Carbon: target position order US, MS, DS. For now run only 2.2 and 4.4 GeV,
# and only Møller and C12 elastic. C12 inelastic and 6.6/8.8/11 GeV are
# intentionally excluded.
for target in c12_us c12_ms c12_ds; do
  for energy in 2200 4400; do
    run_cell "${target}:${energy}:moller"
    run_cell "${target}:${energy}:c12_elastic"
  done
done

echo
if ((${#failures[@]})); then
  echo "${#failures[@]} configuration(s) require another invocation:"
  printf '  %s\n' "${failures[@]}"
  echo "Completed configurations are intact; rerun this script to retry only the remainder."
  exit 1
fi
echo "All requested configurations reached ${event_cap} events."
echo "Use showermax_live.py --run-name RUN_NAME display to review any result."
