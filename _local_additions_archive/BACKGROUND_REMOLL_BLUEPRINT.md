# Background remoll Campaign: Blueprint and Build Plan

## Purpose

Continuously run remoll configurations until detector-plane rates and
kinematic distributions are either:

- measured to the requested precision;
- bounded tightly enough to be negligible; or
- explicitly reported as unresolved within the allowed run budget.

The system must make it easy to see what is running, where the signal is
going, why a configuration has or has not converged, and what it will do next.

## Blueprint

### Minimal architecture

Use four components:

1. **Campaign file**: configurations, priorities, tolerances, and resource
   limits.
2. **Controller**: runs one independent remoll batch at a time, resumes after
   interruption, and remains on one configuration until it closes or exhausts
   its budget.
3. **Analyzer**: reads each untouched ROOT batch and updates mergeable
   event-level statistics and detector-plane distributions.
4. **Operator view**: concise `status`, expanded `inventory`, and detailed
   per-configuration `inspect` output generated from the same summary.

Keep implementation files together in one directory. Each campaign gets one
output directory containing its state, batch ROOT files, and report. Do not
create per-observable scripts or per-batch analysis files.

### Campaign file

Use one commented `campaign.toml` file. TOML is plain text, readable without
special tools, and easy to edit while still being strict enough to validate.
The controller must provide `campaign check` to explain mistakes before a run
starts.

```toml
title = "MOLLER detector-plane campaign"
batch_events = 10000
threads = 2
max_disk_gb = 500
run_on_battery = false

[precision]
relative_target = 0.02
confidence = 0.95
minimum_ordinary_batches = 5
minimum_total_batches = 10
maximum_single_event_fraction = 0.01
minimum_signal_events_for_bound = 30

[[detector_requirement]]
group = "all_rings"
relative_target = 0.02
negligible_fraction = 0.001

[[detector_requirement]]
group = "showermax_regions"
relative_target = 0.03
negligible_fraction = 0.001

[[configuration_matrix]]
name_format = "calibration_c12_{target}_{energy_label}gev_{interaction}_sieve_{sieve}"
targets = ["us", "ms", "ds"]
sieves = ["in", "out"]
interactions = ["moller", "c12_elastic"]
energies_gev = [2.2, 6.6, 8.8]
```

Frequently changed settings belong here. Detector definitions and statistical
methods belong in the implementation and report their version; they should not
be duplicated across campaign files.

### Complete data retention

Every successfully completed remoll ROOT batch is permanent campaign data.
Retain the complete original tree, including:

- `ev` and all generator-level kinematics;
- `part` and all generated particles;
- every recorded `hit`, including ancestry and source information;
- event rate/weight, beam/target data, seeds, and run metadata;
- any other branches written by remoll.

The analyzer writes only one mergeable campaign summary and report index. It
never strips, rewrites, or replaces a batch ROOT file. This preserves the
ability to define new kinematic observables after the campaign has run.

Disk use must be measured during the first tests. If full files are too large,
change what remoll records only after measuring branch and detector costs; do
not discard data speculatively.

### Detector recording

- Record Main Detector Rings 1--6 on the physical quartz entry surfaces.
- Record ShowerMax on the 28 tile-sized, beam-facing stack planes, grouped as
  open, closed, and transition.
- Retain GEM and ring virtual-plane hits as diagnostics only.
- Exclude aggregate detector IDs 28 and 30 from precision and acceptance
  results.
- Record crossings with remoll `surfacehits`: one entry record, no detector
  energy-deposition sum, and no interior-step records.
- Analyze primaries and secondaries separately while retaining their combined
  result.

Surface-hit mode changes only sensitive-detector recording. It does not stop
or simplify Geant4 transport. In particular, tracks are not killed when they
reach a detector: cascades produced in one detector region can reach another
region and are part of the required flux.

Before production, benchmark short runs including the intended final
`surfacehits` configuration.
Record runtime and ROOT bytes per event so the actual cost is known.

### Statistical unit

The independent unit is one generated remoll event, not one detector hit.
For event `i` and detector group `d`, retain at least:

\[
Y_{i,d}^{R}=w_i N_{i,d},
\qquad
Y_{i,d}^{E}=w_i\sum_{j\in(i,d)}K_j.
\]

Hits from one shower remain grouped so their correlations are included in the
uncertainty.

### What is monitored

For every configuration and detector plane/group:

- weighted crossing rate;
- `rate * kinetic energy` response;
- primary and secondary contributions;
- signal-producing event count;
- relative standard error and confidence interval;
- largest-event and largest-batch contribution;
- recent-batch stability;
- weighted energy, position, PID, ancestry, and source-location distributions.

`rate * kinetic energy` is a transport-response proxy, not deposited energy.
Ancestry analysis uses track and parent IDs plus the creator process, volume,
and material fields. The legacy `hit.gen` pointer-valued field is not used.

The report must also show the signal flow through the apparatus, so a weak
ShowerMax result can be interpreted alongside strong rates on upstream planes
or other detector rings.

### Continuous learning loop

Every completed batch must improve the understanding of the same accumulated
answer:

1. Merge its event-level statistics into every detector plane and group.
2. Update precision and stability for all tracked rates and distributions.
3. Locate generator-level kinematic cells contributing the most signal,
   variance, single-event dominance, or unresolved coverage.
4. Choose the next ordinary or biased remoll batch that most improves the
   campaign's worst important result.
5. Preserve the sampling probability and weight correction so biased and
   ordinary batches combine without bias.

The campaign is therefore not a list of unrelated 10k scans. It is one growing
measurement whose next sample is selected from everything learned so far.

### Completion states

Every required detector result receives one visible state:

- `measured`: significant result meeting relative-precision and stability
  requirements;
- `negligible`: confidence upper bound below its absolute importance limit;
- `not_relevant`: excluded by the campaign's declared measurement scope;
- `running`: more batches can resolve it;
- `rare_tail`: insufficient coverage or excessive single-event dominance;
- `statistics_limited`: unresolved at the configured event/time limit;
- `failed`: repeated execution or data-integrity failure.

A configuration completes only when every required detector is `measured`,
`negligible`, or `not_relevant`.

### Near-zero signals

Do not demand relative precision from a result consistent with zero. Use an
absolute upper bound expressed as a configured rate or as a fraction of a
reference detector response. Upstream and neighboring detector-plane rates
must be shown beside the bound to establish where the signal went.

No observed hits alone never proves zero. A result remains `rare_tail` until
its upper bound is sufficiently small or targeted sampling resolves it.

### Operator view

The default status display should fit on one screen:

| Configuration | Running/next | Events | Signal destination | Worst required result | State |
|---|---|---:|---|---|---|

Selecting a configuration should reveal:

- per-detector rates and precision;
- primary/secondary split;
- convergence history by batch;
- signal-flow diagram;
- unstable distributions and dominant events;
- the exact reason for the next action or stopping decision.

Required commands:

```text
campaign start
campaign status
campaign inspect CONFIG
campaign pause
campaign resume
campaign stop
campaign report
```

### Background behavior

- One remoll process and one analysis process at a time.
- Independent seeds and immutable ROOT files for completed batches.
- Low process priority and configurable remoll threads.
- The remoll and ROOT standard streams go directly to bounded files, never
  through a terminal or an in-memory capture.
- Keep one current diagnostic log per active configuration, rotate it at a
  small fixed size, and retain only warnings/errors from successful batches.
- Store one compact run-history table containing seed, event count, start/end
  time, elapsed time, exit status, ROOT path, file size, and analysis status.
- Preserve full logs only for failed batches and cap their number and size.
- Terminal output is one periodically refreshed status view, not an appended
  event or batch transcript.
- Atomic state updates and automatic restart after interruption.
- Optional AC-power, load, temperature, disk-space, and runtime limits.
- A failed batch is quarantined and does not stop the campaign.

### Adaptive sampling

Ordinary batches establish the baseline. The controller then identifies weak
generator-level regions from accumulated `ev`, `part`, and detector-response
data and requests biased remoll batches there. Current proposal coordinates
are thrown angle, azimuth, pre-vertex beam momentum, target depth, and the
outgoing-energy fraction for ep inelastic. Every proposal mixes a targeted
component with a nonzero physical-distribution component and records both PDFs
and their weight ratio. The GPU transport oracle can later improve those
priorities by predicting primary destinations. These tools change where
computation is spent, never what is retained or how convergence is judged.
Periodic ordinary control batches provide the efficiency baseline. A targeted
proposal is retained only when its event variance per unit runtime improves on
that baseline for the detector observable it was intended to resolve;
repeatedly worse proposals are rejected and the controller falls back to
ordinary sampling.

## Build Plan

### Milestone 1: One complete trustworthy batch

- Define and validate `campaign.toml`.
- Establish the authoritative detector-plane ID/group table and configuration
  inventory.
- Run the detector-recording benchmarks, including the final `surfacehits`
  configuration.
- Produce one complete ROOT batch with all required branches and a compact run
  record.

**Gate:** physics output matches a direct remoll run; file size, runtime,
detector IDs, branches, and surface-hit behavior are documented.

### Milestone 2: Accumulating precision and inspection

- Implement event-level rate and `rate * energy` accumulation.
- Accumulate all detector planes, primary/secondary classes, generator cells,
  and weighted kinematic distributions without altering ROOT files.
- Build the one-screen status, complete inventory, and detailed inspection
  views.
- Validate against existing calibration ROOT files, including known
  well-sampled and single-event-dominated cases.

**Gate:** adding batches monotonically updates one reproducible campaign
answer; poor samples remain visibly unresolved and every number can be traced
back to contributing events and ROOT files.

### Milestone 3: Safe unattended campaign

- Implement batch launch, seeds, validation, state transitions, pause/resume,
  failure isolation, and resource limits.
- Implement relative convergence, negligible upper bounds, rare-tail checks,
  distribution stability, and statistics-limited outcomes.
- Implement bounded output, compact run history, recovery, and disk controls.
- Test interruption, failed batches, malformed ROOT files, and multi-day log
  behavior before leaving it unattended.

**Gate:** a stopped or crashed campaign resumes without duplicate work;
near-zero detectors finish only with defensible bounds; terminal and log
memory remain bounded.

### Milestone 4: Adaptive production

- Identify under-sampled generator kinematics after every batch.
- Add remoll generator biasing with recorded sampling probabilities and weight
  corrections.
- Allocate batches according to detector-plane uncertainty and scientific
  priority, including targeted resolution of rare tails.
- Evaluate GPU-guided priorities and roulette/splitting against ordinary
  remoll closure runs.

**Gate:** adaptive and ordinary estimates agree within combined uncertainty,
the adaptive campaign reaches the same precision faster, and no physical
source region is silently assigned zero probability.

## Non-negotiable Rules

- Never treat detector hits as statistically independent events.
- Never declare zero from an empty detector alone.
- Never declare completion from relative error alone.
- Never hide unresolved rare tails inside a combined total.
- Never replace or thin an accepted remoll ROOT batch after it is produced.
- Never let adaptive biasing assign exactly zero probability to physical
  source regions without a proven exclusion.
- Every completed result must retain its inputs, stopping reason, uncertainty,
  and reproducible batch list.
