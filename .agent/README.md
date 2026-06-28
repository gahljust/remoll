# remoll Quick State

Read this first. Terse on purpose.

## Goal

ShowerMax rates. Found rare-event rate runaway. Study generator math first.

## Core Meaning

- `rate = effective cross section * effective luminosity / Nthrown`.
- Bias weight multiplies cross section/rate first.
- Code: `src/remollPrimaryGeneratorAction.cc:204-219`.
- Luminosity = beam flux + target areal density/composition.
- Target/radiation: `src/remollBeamTarget.cc`.
- `beamp`: beam momentum at forced vertex, after sampled pre-vertex loss.
- `thcom`: generator-defined angle. Moller: COM angle. C12/Al elastic: thrown lab angle.
- `Q2`: positive `-q_mu q^mu`, stored in MeV^2.
- More rate detail: `.agent/rate information.md`.

## Found Problem

Elastic C12 uses bare Mott/Rutherford pole. Roughly `dSigma/dOmega ~ 1/Q^4`.
Pre-vertex radiation can leave electron near rest. Then Q2 tiny. Weight huge.
Few events control rate. Plain Monte Carlo not converged.

Do not say beam itself mostly low energy. Nominal beam loses energy in sampled
pre-vertex radiator model.

Do not "fix" by deleting a few giant events. Tail continues. Cut changes answer.

## Hard Numbers: 2.2 GeV C12 Optics1, 0.1-5 deg

100M generator-only sample:

- Total sampled rate: `8.04727396480e11 /s`.
- Rate-weight Neff: `50.91`.
- Largest event: `1.09177643008e11 /s` = 13.57% total.
- Largest event: beamp `9.286 MeV`, theta `0.1682 deg`, Q2 `7.4307e-4 MeV2`.
- beamp < 10 MeV: 263 events; 20.71% rate.
- beamp < 20 MeV: 550 events; 22.24% rate.
- beamp < 80 MeV: 2230 events; 23.52% rate.
- Q2 < `9e-4 MeV2`: 8 events; 15.87% rate.
- Keeping Q2 >= `9e-4 MeV2`: `6.77020806368e11 /s`.

Bare numerical quadrature:

- Total: `1.98057488555e12 /s`.
- Doubled nodes: `1.98059338803e12 /s`. Numerically stable.
- Rate below 10 MeV: 65.0%.
- Rate below 80 MeV: 68.9%.
- Rate below Q2 `9e-4 MeV2`: 64.7%.

Meaning: 100M Monte Carlo badly undersampled bare-model low-energy pole.

Wentzel-Moliere screened quadrature:

- Factor: `[Q2/(Q2 + qa2)]^2`.
- Carbon `qa2 = 6.6664e-5 MeV2`; `qa = 8.1648 keV`.
- Total: `8.52069016488e11 /s`.
- Below 10 MeV: 19.4%.
- Below 80 MeV: 27.8%.
- Below Q2 cut: 18.4%.

Screening is candidate fix. Implementation consistent. Not externally benchmarked
yet. Production remoll still unscreened.

## Biased Sampling

Added importance sampling:

- `/remoll/bias/beamp/mode uniform`
- `/remoll/bias/beamp/min ...`
- `/remoll/bias/beamp/max ...`
- `/remoll/bias/thcom/mode uniform`
- `/remoll/bias/thcom/min ...`
- `/remoll/bias/thcom/max ...`
- Weight = physical PDF / sampling PDF.

Code:

- Beam bias: `include/remollBeamTarget.hh`, `src/remollBeamTarget.cc`.
- Angle bias: `include/remollVEventGen.hh`, `src/remollVEventGen.cc`.
- Stored weights: `include/remollEvent.hh`, `src/remollEvent.cc`,
  `include/remolltypes.hh`.
- Generator calls: `src/remollGen*.cc`.

100k bias test: beamp 0-0.1 GeV, angle 0.1-0.4 deg.

- Total weighted rate: `1.26944e12 /s`.
- Neff: `451.3`.
- Lowest 10 MeV: 9.45% throws, 94.32% rate.
- Low Q2 region still weight dominated.

Bias helps sample chosen phase space. Bias cannot repair bad physics model.

## Generator Only

Two paths:

- remoll flag: `/remoll/evgen/generatorOnly true`.
- Fast executable: `genonly_scan.cc`, built as `build/genonly_scan`.

Fast scan does generator math. No Geant4 transport. Compact Float branches:

`xs rate beamp_GeV thcom_deg Q2 W2 A vx vy vz p0_GeV th0 ph0 radlen travelled`

Runner: `_local_additions_archive/current/run_genonly_scan.sh`.
Macro: `_local_additions_archive/current/macros/gen_only.mac`.
100M file: `_local_additions_archive/current/gen_scan.root`.

## Numerical Solver

Not Monte Carlo. Deterministic Gauss-Legendre quadrature.
Integrates vertex depth, radiative loss, angle for fixed beam/target setup.
Production remoll untouched by solver.

Path: `_local_additions_archive/current/rate_solver/`.

Main:

- `c12_rate_solver.cc`: completed elastic-C12 integral.
- `gauss_legendre.hh`: quadrature.
- `models/radiative_loss_model.hh`: copied radiator math.
- `models/c12_models.hh`: C12 models + candidate screening.
- `models/upstream_snapshot/`: exact source snapshots.

Other model math copied and event-tested. Full multidimensional rate integrators
for them are NOT completed.

Validation:

- `verify_model_snapshots.sh`: snapshots vs current remoll byte-for-byte.
- `validate_all_models.sh`: vanilla remoll + genonly + copied math, 100 events.
- Compared stored inputs/outputs event by event: beamp, angle, target Z/A,
  outgoing momentum as needed; recomputed xs, Q2, W2.
- All tested models passed configured precision.
- ep elastic xs not fully checked: hidden radiation-sampling weight not stored.
- Validation proves copy fidelity. Does not prove physical validity.

Screening tests:

- `verify_screening.cc`.
- Bare copy vs remoll arithmetic: max relative difference `4.48e-16`.
- Screened copy vs same arithmetic + factor: `4.48e-16`.
- Screening off = bare bitwise.

## remoll Generator Map

Registry: `src/remollPrimaryGeneratorAction.cc:51-65`.

- `moller`: `src/remollGenMoller.cc`, `include/remollGenMoller.hh`.
  Free target electron. COM angle. Q2 stored as 0 by current code.
- `elastic`: `src/remollGenpElastic.cc`, `include/remollGenpElastic.hh`.
  Proton elastic. Custom radiation sampling. Hard `fE_min = 80 MeV`.
  Low-angle suppression tied to multiple-scattering scale.
- `inelastic`: `src/remollGenpInelastic.cc`, `include/remollGenpInelastic.hh`.
  Christy-Bosted. Fit code: `include/christy_bosted_inelastic.h`.
- `elasticC12`, `quasielasticC12`, `inelasticC12`:
  `src/remollGenC12.cc`, `include/remollGenC12.hh`.
- `elasticAl`, `quasielasticAl`, `inelasticAl`:
  `src/remollGenAl.cc`, `include/remollGenAl.hh`.
- `pion`: `src/remollGenPion.cc`, `include/remollGenPion.hh`. Wiser model.
- `flat`: `src/remollGenFlat.cc`, `include/remollGenFlat.hh`.
- `beam`: `src/remollGenBeam.cc`, `include/remollGenBeam.hh`. No scattering xs.
- `external`: `src/remollGenExternal.cc`, `include/remollGenExternal.hh`.
  Reads ROOT events/weights.
- `pion_LUND`: `src/remollGenLUND.cc`, `include/remollGenLUND.hh`.
  Reads LUND.
- `hyperon`: `src/remollGenHyperon.cc`, `include/remollGenHyperon.hh`.

Shared generator base:

- `src/remollVEventGen.cc`, `include/remollVEventGen.hh`.
- Vertex/radiation/MS: `src/remollBeamTarget.cc`, `include/remollBeamTarget.hh`.
- Rate normalization: `src/remollPrimaryGeneratorAction.cc`.

## Model Audit Status

- C12 elastic: demonstrated low-Q2 bare-Mott problem.
- Al elastic: same bare Coulomb structure. Likely same screening need. Not yet
  integrated/benchmarked.
- C12/Al QE and inelastic: Mott factors + fit domains. Audit needed. Do not apply
  elastic screening blindly.
- ep elastic: already special low-energy/angle machinery. Physics rationale and
  continuity need review.
- ep inelastic: Christy-Bosted domain limits need audit.
- Moller: free-electron approximation. Binding/screening relevance needs audit.
  Do not claim solved merely because normal angle cuts exist.
- Pion: separate Wiser model. Different issue class.

Multiple scattering is Geant4 transport AND remoll has pre-vertex estimate.
It was not proven to cause C12 rate runaway. Main proven chain is radiative loss
-> low beamp -> tiny Q2 -> bare Mott pole.

## ShowerMax Checks

Biased full run: 1M events. Thresholds: e+/e- kinetic >=5 MeV; gamma >=10 MeV.

- Detector 30 visible: 235 events; `1.06509e6 /s`; Neff `1.00`.
- One event = 99.82% of detector-30 rate. Unusable convergence.
- Individual SM planes visible: 16 events; `31.15 /s`; Neff `2.51`.
- Largest plane event = 61.2% rate.

Transport controls:

- 100k forward 3 MeV electrons: 45 det30 events, 3 plane events, zero visible.
- 100k forward 2.197 GeV photons: zero det30/plane hits.

These are narrow checks. Not proof full radiative tail is harmless.

Files/macros:

- `_local_additions_archive/current/showermax_bias_test.root`
- `_local_additions_archive/current/macros/showermax_bias_test.mac`
- `_local_additions_archive/current/showermax_transport_*`

## Evidence Packet

Start:

- `_local_additions_archive/current/rate_evidence/WALKTHROUGH.md`
- `_local_additions_archive/current/rate_evidence/RECEIPTS.md`
- Receipts: `_local_additions_archive/current/rate_evidence/receipts/`.
- Email figures: `_local_additions_archive/current/rate_evidence/figures/`.
- Rebuild figures: `rate_evidence/make_email_figures.sh`.
- Rebuild receipts: `rate_evidence/reproduce.sh`.

Current email figures:

1. Original anomaly.
2. 100M: cumulative rate vs highest-weight event percentage.
3. Biased 100k: same metric, separate figure.
4. Solver/remoll xs parity.
5. Three totals: MC, bare quadrature, screened quadrature.

## Research Notes

- `_local_additions_archive/claudes_ideas.md`
- `_local_additions_archive/gpt_deep-research-report.md`
- `_local_additions_archive/gemini's_ideas.txt`
- `_local_additions_archive/MOLLER_TDR-Final.md`

Agent reports are leads. Verify primary sources before physics claim.

## Next Work

1. External benchmark Wentzel-Moliere screening scale/formula against trusted
   elastic electron-atom calculation/data (e.g. ELSEPA literature).
2. Implement screened C12 locally only after benchmark.
3. Repeat generator integral + biased detector tracking.
4. Audit Al elastic next.
5. Audit QE/inelastic/Moller/ep separately. No universal blind factor.
6. Keep old ROOT runs as bare-model evidence. Do not call them final rates.

## Safety

- `_local_additions_archive/` intentionally not tracked by remoll git.
- It has separate local history.
- Do not delete large ROOT files casually.
- Do not change production remoll formulas while doing solver-only tests.
- Worktree may be dirty. Never undo user changes.
