# Background remoll Campaign

This controller runs one remoll configuration until its required detector
results close or reach their declared budget. Only then does it advance to the
next configuration. Complete ROOT batches are retained without thinning or
rewriting.

## Operator files

- `campaign.toml`: configurations, priorities, precision, and resource limits.
- `detectors.toml`: authoritative precision and diagnostic detector groups.
- `base.mac`: shared remoll geometry, fields, and transport controls.
- `campaign.py`: validation, scheduling, analysis, and status commands.
- `analyze_batch.C`: event-level ROOT accumulation and generator-cell analysis.
- `adaptive_replay/`: the experimental parent-surface replay scanner and all
  of its ROOT helpers. Generated files belong under the ignored `runs/` tree.

The editable inventory expands to 39 configurations: three 11 GeV LH2
production interactions and 36 carbon calibration combinations from three
energies, three foil positions, two sieve states, and two interactions.

## Commands

Run from the remoll repository root:

```bash
python3 _local_additions_archive/background_campaign/campaign.py check
python3 _local_additions_archive/background_campaign/campaign.py inventory
python3 _local_additions_archive/background_campaign/campaign.py status
python3 _local_additions_archive/background_campaign/campaign.py inspect CONFIG
python3 _local_additions_archive/background_campaign/campaign.py run-one CONFIG --dry-run
python3 _local_additions_archive/background_campaign/campaign.py start
python3 _local_additions_archive/background_campaign/campaign.py pause
python3 _local_additions_archive/background_campaign/campaign.py resume
python3 _local_additions_archive/background_campaign/campaign.py stop
```

### Short controller test

`campaign_test.toml` contains only 11 GeV LH2 Moller and writes to the separate
`runs/controller_test` directory. It uses 1,000-event batches and a temporary
10% precision target:

```bash
python3 _local_additions_archive/background_campaign/campaign.py \
  --campaign _local_additions_archive/background_campaign/campaign_test.toml start
```

The current completed exercise reached its 20,000-event test budget. The old
group-only rule would have stopped after 6,000 events, but the tally-reliability
checks correctly kept it running and ultimately reported `statistics_limited`.
These are controller-validation results, not final physics rates or production
precision.

`start` is continuous. It does not rotate through configurations batch by
batch. `pause` takes effect after the active batch; `stop` terminates the
active remoll process and leaves the campaign resumable.

## What closes

For Rings 1-6, every individual ring tile, ShowerMax open, closed, and
transition regions, and every individual ShowerMax plane, the analyzer
uses one generated remoll event as the independent statistical unit. It tracks
crossing rate for all tallies and `rate * kinetic energy` for ShowerMax,
including all hits from the event. Group tallies also report primary and
secondary contributions separately.

A required result closes only after:

- the ordinary-pilot and total-batch minimums are met;
- relative standard error is below the configured target;
- no single event exceeds the configured contribution limit; and
- the latest three batches agree with the earlier estimate within the
  configured stability threshold;
- ordinary-control histories pass the configured VOV, error-scaling, FOM,
  high-score-tail, and next-largest-history checks.

A weak nonzero result can instead close as negligible when its confidence
upper bound is below the configured fraction of the strongest detector-group
result. At least 30 signal events are required for that bound. An empty result
is never declared zero; it remains a rare tail until resolved or explicitly
reported `statistics_limited` at the event/batch budget.

Ordinary and biased batches have different variances. The pooled estimate is
the event-count-weighted mean, while its variance is assembled from each
batch's own within-batch variance. This is valid for the changing adaptive
proposal and avoids treating hits from one shower as independent samples.

## Tally reliability

The analyzer retains the first four event-score moments and the 201 largest
event scores for every required group and tile tally. Reliability is evaluated
from ordinary physical-distribution control batches, independently of the
pooled estimate that uses all correctly weighted batches. The reported suite is
modeled on the MCNP history-score checks:

- relative-error behavior versus history count, expected near `N^-1/2`;
- variance of the variance (VOV), expected to decrease near `N^-1`;
- figure-of-merit stability over cumulative ordinary controls;
- a Hill power-law index from the 201 largest nonzero scores, required to be
  greater than 3 by default; and
- the fractional mean change if the next history equals the largest score seen.

The slopes and thresholds are editable in `[reliability]` in the campaign
file. `inspect CONFIG` shows group diagnostics and tile state counts;
`inspect CONFIG --tiles` prints every individual tile tally. Passing these
checks is evidence of stable sampled history scores, not proof that an unseen
phase-space contribution cannot exist.

Reference: LANL, *An MCNP Primer*, 2024, section 5.6, “Statistical Analysis of
Tally Results”: <https://mcnp.lanl.gov/pdf_files/Book_MonteCarlo_2024_ShultisBahadori_AnMCNPPrimer.pdf>.

## Adaptive batches

After five ordinary pilot batches, the controller finds the unresolved
detector observable with the largest relative uncertainty or event dominance.
It then targets the generator cell with the largest conditional second moment.
Available coordinates are thrown angle, azimuth, pre-vertex beam momentum,
target depth, and, for ep inelastic, outgoing-energy fraction.

At least every fifth post-pilot batch is an ordinary control. For the exact
observable a biased proposal was intended to improve, the controller compares

\[
C=s^2\frac{t}{n},
\]

where `s^2` is its pooled event variance and `t/n` is runtime per event. Lower
`C` means more precision per unit runtime. After three trials of the same
axis/cell/observable proposal, it is rejected if it produces no target
responses or if its cost is at least 1.25 times the ordinary-control cost.
Rejected proposals are not selected again. If no useful biased proposal
remains, sampling automatically returns to the ordinary distribution. All
rejected batches remain valid weighted campaign data.

Every adaptive proposal is a full-support mixture:

\[
q(x)=\epsilon p(x)+(1-\epsilon)g_k(x),
\qquad
w_{\mathrm{bias}}(x)=\frac{p(x)}{q(x)}.
\]

The default `epsilon = 0.20` preserves a physical-distribution component over
the complete phase space. No cell is cut away. Moller electron pairs remain
one remoll event with one shared event weight. The ROOT `ev` record stores the
total bias weight and each component's weight, physical PDF, and sampled PDF.

## Experimental parent-surface replay

`adaptive_replay/adaptive_parent_replay.py` is the single continuation point
for the recent ShowerMax convergence experiment. It analyzes all 62 required
tile/group rate and energy-rate observables, locks onto one unresolved
observable until its complete convergence gates pass, runs exact-state pilots,
adds correctly weighted replay batches to one resumable accumulated ROOT file,
retains ordinary-control batches for the MCNP-style reliability tests, and
leaves the input ROOT file untouched.

The central correction is **parent-first ancestry**. If an important detector
hit is a secondary, its own downstream surface crossing is not a valid biasing
state: replaying that state only resamples transport after the rare interaction
has already happened. Candidate replay states now begin with `hit.mtrid` and
walk farther up the recorded ancestry only as needed. Primary hits may still
use their own upstream crossing. If no parent/ancestor crosses the requested
recorded surface, the source builders fail explicitly instead of silently
substituting the child.

Current checkpoint status:

- The restart/recovery logic, 62-observable convergence accounting,
  original-plus-replay estimator, exact-state pilot learning, full-support
  proposal weighting, and multithreaded external-rate fix are implemented.
- A 1,000-event target-primary bias test was intentionally retired. It put 449
  events into the selected target kinematic neighborhood but produced no hit
  in the sparse ShowerMax tile, demonstrating that the rare secondary must be
  attacked at its responsible parent surface rather than at the original
  target primary.
- The representative sparse gamma was recorded at GEM surfaces, but its causal
  parent lineage was not. `run_twofold_poc.py` now automates recovery from the
  stored event RNG state, records the responsible physical-volume boundary,
  follows parent creation steps in causal order, and builds a proposal from
  physically observed neighboring entries. Integrating that resolver into the
  continuous `adaptive_parent_replay.py` controller remains the next step.

Run the current scanner from the repository root (US carbon 2.2 GeV is the
default input):

```bash
python3 _local_additions_archive/background_campaign/adaptive_replay/\
adaptive_parent_replay.py --check

python3 _local_additions_archive/background_campaign/adaptive_replay/\
adaptive_parent_replay.py --batch-events 1000 --pilot-events 1000 \
  --threads 4 --target-rse 0.05 --max-batches 100
```

The continuous scanner still stops with a clear missing-parent-surface error
when an old file lacks the crossing. Run the two-fold proof below to recover
and test that boundary automatically; it is not yet wired into the scanner's
batch loop. Ancestry/creation-volume summaries live alongside the scanner in
`adaptive_replay/ancestry_attribution.C`.

### Two-fold primary/secondary proof

`adaptive_replay/run_twofold_poc.py` performs the branch decision directly on
an existing ROOT file. It runs all 62 statistics, writes cumulative tile RSE
and Neff versus history count through `diagnose_tile.C`, ranks the histories
dominating the worst result, and classifies the leading hit. A primary result
is explicitly handed off to the campaign generator-neighborhood machinery
(that handoff is not exercised by this secondary-dominated file). For a
secondary, the proof exports its stored random state and a random physical control sample,
deterministically replays those histories with volume interaction recording,
and walks the causal ancestry until it finds a real ancestor entry that occurs
before the selected secondary-production chain.

For the US carbon 2.2 GeV sieve-out file, the worst observable is
`tile_71630/rate`: 29.86% RSE, Neff 11.22, with one history carrying 29.29% of
the total. It is a 1.76 MeV gamma produced by bremsstrahlung in
`twobounce_long`. Exact replay resolves the causal chain as electron 2436 to
gamma 2525 to electron 2546 to gamma 2695; electron 2436 is the first ancestor
with a volume entry before that chain. In an equal-cost 2,000-history pilot
using 68 observed forward electron-entry states, ordinary physical-bank replay
gave 1 useful history (Neff 1.00, RSE 100%), while the full-support targeted
mixture gave 10 useful histories from 42 distinct source states (Neff 4.24,
RSE 48.5%). This is a 10x useful-history gain but not convergence or a final
normalization; the finite observed entry bank remains a pilot approximation.

```bash
python3 _local_additions_archive/background_campaign/adaptive_replay/\
run_twofold_poc.py --bank-events 200 --replay-events 2000 --threads 4
```

## Data and output

Each accepted batch contains the full remoll tree: `ev`, `part`, `hit`, `rate`,
`bm`, and the other standard branches. A compact JSON record stores the exact
macro, seed, proposal, integrity result, elapsed time, alerts, batch moments,
and adaptive-cell moments. Successful verbose logs are deleted after alerts
are extracted. Active and failed logs are capped at 5 MB, and remoll output is
never streamed through the terminal or held indefinitely in memory.

## Validation completed

- The 39-item inventory, target/sieve macros, generator commands, screening,
  detector IDs, and all five proposal command sets pass validation.
- Ordinary and theta-mixture Moller generator rates agreed within `1.53`
  combined standard errors in independent 100,000-event tests.
- Ordinary and outgoing-energy-mixture ep-inelastic rates agreed within `1.02`
  combined standard errors in independent 50,000-event tests.
- Tiny transported batches exercised ordinary, theta, phi, beam-momentum,
  target-depth, and outgoing-energy modes with valid ROOT integrity and stored
  PDF/weight ratios. Validation ROOT files and logs were removed afterward.
