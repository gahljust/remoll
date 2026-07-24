# ShowerMax live response campaign

This is an analysis/controller layer. It does not modify remoll physics. Every
batch is a separate, complete ROOT file and is accepted only after ROOT
integrity and response analysis succeed. Restarting reuses all accepted batches.

All 28 ShowerMax stack planes are displayed and tallied individually. Aggregate
closed, open, transition, and total ShowerMax tallies are formed by summing the
applicable plane responses inside each generated event before calculating
moments. Response is the existing ShowerMax photoelectron lookup multiplied by
the event `rate`; `rate` contains remoll's non-azimuthal kinematic-bias
correction.

The runner reads the selected q01--q99 target-level window from the compact COL2
kinematics table and applies full-support mixture biasing. All 28 ShowerMax
response planes (73001--73028) are enabled and use entry-surface-only hit
storage; unrelated remoll sensitive detectors are disabled.

Azimuth is not biased. The generator covers 0--360 degrees physically, and no
`/remoll/bias/phi` command is emitted. The septant column in the compact
kinematics table only selects the target-level theta/energy/vertex window; it
does not restrict or bias generated azimuth.

Detector 30, the full ShowerMax virtual plane, is also recorded with
`surfacehits`. The dashboard draws its crossing-rate map and all 28 response
planes on a common global x/y scale using their actual positions and
orientations.

The beam-facing and far-facing main-detector surface planes are recorded for
all 14 segments. Dashboard tabs provide Ring 1 through Ring 6 and a combined
main-detector view. These systems use event-level rate-weighted scores only;
the ShowerMax tab alone applies the photoelectron-response lookup.

Every ring and the combined main detector also have closed, transition, open,
and total convergence tallies. Ring hits are assigned from their actual global
azimuth using the established remoll sevenfold partition, rather than inferred
from detector-number ordering.

Simulation and display analysis are separate producer/consumer paths. remoll
immediately starts the next production batch after atomically marking its ROOT
file complete. A background ROOT worker analyzes completed batches and updates
the dashboard; display work never blocks the simulation producer.

Run `python3 showermax_live.py check` before a campaign. Use `start` to run and
serve the dashboard, Ctrl-C or `pause` for a clean batch-boundary pause, and
the same `start` command to resume.

From the repository root, for the current US 2.2 GeV carbon-elastic,
sieve-out selection:

```bash
python3 _local_additions_archive/background_campaign/showermax_live/showermax_live.py \
  --cell c12_us:2200:c12_elastic --septant 0 --sieve out check

python3 _local_additions_archive/background_campaign/showermax_live/showermax_live.py \
  --cell c12_us:2200:c12_elastic --septant 0 --sieve out start
```

The second command opens `http://127.0.0.1:8765`. Dashboard pause/resume takes
effect at a batch boundary. After a process exit, rerun the identical `start`
command; accepted `batch_*.json`/ROOT pairs are reused.

Production ROOT files may contain any positive multiple of 1,000 events:

```bash
python3 _local_additions_archive/background_campaign/showermax_live/showermax_live.py \
  --cell lh2:11000:moller --septant 0 --batch-events 10000 start
```

Analysis and convergence checkpoints remain fixed at 1,000 histories. A
10,000-event ROOT file therefore contributes ten independent checkpoints, and
production batch size may change when resuming the same campaign. To separate
display serving from production completely, use `start --no-dashboard` in one
terminal and `serve` with the same selection in another.

To display saved data without running remoll:

```bash
python3 _local_additions_archive/background_campaign/showermax_live/showermax_live.py \
  display
```

This opens the most recently updated non-diagnostic campaign. Use
`--run-name NAME display` for a particular saved campaign. The snapshot and
all trend graphs are rebuilt from every accepted batch record whenever
`start` or `display` is launched.

The dashboard's **Setup** menu discovers every campaign directory under the
configured output path. Selecting another setup immediately switches all maps,
tile statistics, convergence trends, and particle spectra. The same menu also
shows event counts and whether a setup is saved, paused, stopped, diagnostic,
or actively running. A selected live campaign continues updating every two
seconds without restarting the viewer.

## Full batch campaign

`run_full_batch.sh` runs LH2 ep elastic and ep inelastic, followed by Møller and
C12 elastic at 2.2 and 4.4 GeV for the US, MS, and DS carbon targets. The
existing LH2 Møller campaign is left untouched. This is 14 queued
configurations, run sequentially without opening dashboards. Every new
configuration uses the sieve-out geometry and stops at 150,000 histories:

```bash
_local_additions_archive/background_campaign/showermax_live/run_full_batch.sh
```

The script is restartable and reuses accepted batches. A fresh configuration is
produced as one 150,000-event ROOT file. If accepted batches already exist, the
next ROOT file contains only the remainder needed to reach 150,000. Analysis
and convergence checkpoints remain fixed at 1,000 histories regardless of
ROOT-file size. `BATCH_EVENTS` and `EVENT_CAP` may override the production-file
size and cap; both must be multiples of 1,000.

Before transport begins, all 14 cells are preflighted for valid generators and
finite, correctly ordered bias ranges. Degenerate physical variables are left
unbiased. A runtime failure is retried by the campaign runner; if it still
fails, the batch script records that cell, continues through the remaining
queue, and prints a retry summary at the end. Rerunning the script preserves
completed work and retries only unfinished histories.
