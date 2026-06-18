# CLAUDE.md — ShowerMax analysis quick reference

Local ShowerMax/`A_PV` study on remoll output. See `.agent/working guide.md` (layout,
rules) and `.agent/rate information.md` (rate/Neff details). Work lives in
`_local_additions_archive/` (own git, gitignored from main repo); active area is
`_local_additions_archive/current/`. Never edit/delete `*.root` (read-only data).

## Reading ROOT files HERE (no ROOT/uproot/reroot available)
The sandbox has no ROOT, no uproot (no internet), and the built `reroot` is a macOS
binary that won't run on the Linux VM. Read ROOT directly with a small numpy parser:
- Scan the file for basket keys: find `b"\x07TBasket"`; the key header starts 34 bytes
  before it (64-bit seeks, fVersion>=1000). Validate `fSeekKey == keystart`.
- Decompress each basket: payload = `raw[ks+fKeylen : ks+fNbytes]`; if compressed it's
  a 9-byte ROOT header + zlib stream (`zlib.decompress(comp[9:])[:fObjlen]`).
- Simple leaf (double/bool/int): the decompressed buffer is the values, big-endian
  (`>f8`, `u1`, `>i4`). Concatenate baskets sorted by `ks` = entry order.
- Variable (vector) branch e.g. `hit.det`: read extended TBasket header right after the
  title (`ver h, bufsize i, nevbufsize i, nevbuf i, last i`). Buffer splits at
  `last-fKeylen`: `[data][entry-offset array]`. Offsets = `off[1:]` (drop leading count);
  per-event element span = `body[(off[k]-kl)//size : (off[k+1]-kl)//size]`.
- Use `mmap` for multi-GB files; guard rare per-event length desyncs (skip mismatched).
Working scripts are in the scratch outputs dir (sector_check / cmp / cmpcut patterns).

## Key data facts (verified, easy to get wrong)
- **`thcom` = the THROWN lab scattering angle**, NOT center-of-mass. It ranges exactly
  over `thmin..thmax` (set by `/remoll/evgen/thmin|thmax`). Do NOT convert it. (remoll
  C12 gen samples lab `th` flat in cos, then `SetThCoM(th)`.)
- **`mtrid == 0` ⇒ primary** (Geant4 ParentID 0). `hit.mtrid = track->GetParentID()`
  (`remollGenericDetector.cc`). Primaries have empty `creator_process_name`.
- ShowerMax dets: single broad plane = **30**; 28 stack modules = **70030..72730**
  (step 100). Region sets by det ID: open `{70030,70430,70830,71230,71630,72030,72430}`,
  closed `{70230,70630,71030,71430,71830,72230,72630}`, rest = transition.
- **Det-ID order is scrambled vs physical position** (`70030`=`singledet_15`, etc.;
  `mollerParallel.gdml:755`). Trust per-module RATES, distrust geometric LABELS. Real
  azimuth: read each plane's `(x,y)` from its placement in `mollerParallel.gdml`.
- `rate` = per-event weight (events/s), already /nthrown. SM rate = Σ over selected hits
  of rate. Per-event weight `W = rate * (selected hits in event)`.

## Analysis findings / methods
- **Precision metric:** `Neff = (ΣW)²/ΣW²`, precision `= 1/√Neff`, on per-EVENT W.
  Exact for independent events. BUT Neff LIES under heavy-tailed weights — it only sees
  sampled weights, so two runs can show good Neff yet disagree. Cross-check with the
  chunk/convergence view (sum rate in blocks of N events, INCLUDING empty blocks; the
  old code dropped empties — a bug).
- **Heavy tail:** elastic dσ/dΩ ~ 1/sin⁴(θ/2) diverges at small angle, so a few tiny-θ
  events dominate the rate and wreck Neff/convergence. Throwing more events makes the
  rate climb, not converge (saw 11 GeV rate 4.6e7→6.8e7 from 100k→1M). FIX = θ cut, not
  more events.
- **θ acceptance band** (where events reach ShowerMax): take `thcom` percentiles (1–99)
  of events with an SM hit. Same for sieve in and out (sieve only masks position, not
  angle). Lower edge ~0.30° (raw thcom) is the real cut; cutting there took 100k 11 GeV
  Neff 449→2237 (4.7%→2.1%). Upper edge runs to the throw cap (~1.25). Band shifts with
  target vertex (US/MS/DS) and energy → measure per config from a sieve-out run.
- **Sieve out vs in:** band identical; YIELD hugely different — sieve-out passes ~25–35%
  of throw vs ~0.1–0.4% sieve-in (holes are a geometric position mask). Use sieve-out for
  efficient high-yield runs; sieve-in stays low-yield.
- **7-fold (septant) symmetry:** sieve-IN breaks it badly (hole mask → factors of 20+
  between modules). sieve-OUT is ~symmetric: 7 coil sectors agree to ~9% (was χ²/ndf~19,
  so a real ~3–5% residual). So one-septant×7 is viable for sieve-out to ~few-%, NOT for
  sieve-in. φ control = `/remoll/evgen/phmin|phmax`, correctly rate-normalized
  (`remollGenC12.cc`), but a single φ window biases the full ring (signal spread over all
  7 sectors) — keep full 2π unless doing a per-sector scheme.
- **Phase-space recycling:** record particles at an SM virtual plane and restart from
  there (`remollGenExternal`, reads a remoll ROOT file at a detid w/ zOffset). Replay =
  exact, handles all secondary origins for free, but capped at source stats. Fit+resample
  (KDE/GMM/flows) only to oversample; must fit the JOINT density per particle-type
  component (not independent 1D marginals), and can't exceed source info in rare tails.

## Gotchas
- The new `carbon_US_*_1M.root` sieve-out runs have a pathological weight tail (a couple
  events = ~20%+ of total rate, odd `beamp`), likely beam-loss events near the dσ/dΩ
  divergence. 1M there is worse than θ-cut 100k. Cap/understand those before trusting 1M.
- Acceptance windows in `current/macros/neff_carbon_acceptance_window.mac` were derived
  with the wrong thcom→lab conversion; correct lower edge is ~0.30° (raw thcom). Redo.
- Make plots/region rates from per-MODULE det sums; the region split is correct in data
  (open ≫ transition ≫ closed for primaries) — equal-region results come from buggy
  summaries, not the data.
