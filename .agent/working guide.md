# Agent Working Guide

This repository is `/Users/justingahley/G4/remoll`. The work here is a local ShowerMax/GEM analysis effort for MOLLER ShowerMax `A_PV` extraction, using remoll simulation output and local ROOT/C++ analysis macros.

## Ground Rules

- Never move, rename, edit, overwrite, or delete ROOT files. Treat `*.root` files as read-only data products.
- Never modify base remoll source, geometry, macros, or build files unless Justin explicitly asks for that exact change.
- Do not compile, run long tests, run remoll, or create extra output files unless Justin explicitly asks for it.
- Do not assume this is a normal autonomous coding workspace. Justin asks for specific help and usually wants exactly that scope only.
- Prefer explanations and references to code/files/lines before making changes.
- When asked to make files, keep them local to the requested project/helper area unless told otherwise.
- Use C++/ROOT `.C` macros for analysis work here unless Justin explicitly asks for another language.
- Keep generated study files under `_local_additions_archive/showermax-gem-analysis/` or another explicitly requested local additions area.

## Project Areas

- Main local project scaffold:
  - `_local_additions_archive/showermax-gem-analysis/`
- Local analysis helpers:
  - `_local_additions_archive/plane_analysis/`
  - `_local_additions_archive/run_helpers/`
- Hidden agent notes:
  - `.agent/`
- Rate reference:
  - `.agent/rate information.md`
- Relevant project docs:
  - `_local_additions_archive/showermax-gem-analysis/README.md`
  - `_local_additions_archive/showermax-gem-analysis/docs/run_matrix.md`
  - `_local_additions_archive/showermax-gem-analysis/docs/open_questions.md`
  - `_local_additions_archive/showermax-gem-analysis/docs/deliverables.md`
  - `_local_additions_archive/plane_analysis/README.md`
- MOLLER TDR copy:
  - `_local_additions_archive/MOLLER_TDR-Final.md`
- remoll rate notes from prior work:
  - `guide/rateCalc.md`

## Current Physics/Analysis Focus

The main study is ShowerMax virtual-plane hit content and later ShowerMax signal/background corrections. The current chain is:

1. Run remoll production/calibration samples.
2. Read virtual-plane hits near ShowerMax.
3. Use the remoll `rate` branch as the per-hit event-rate weight.
4. Build rate-weighted hit maps by process, detector, particle type, energy, and position.
5. Later convert those hits to PE or PE-equivalent signal using a ShowerMax response model.
6. Compute process fractions and correction terms for ShowerMax `A_PV`.

The rate details, including code references and dimensional analysis, are in `.agent/rate information.md`.

## ShowerMax Detector Notes

The ShowerMax detector system is in:

- `geometry/showermax/showerMaxDetectorSystem.gdml`

The sensitive ShowerMax stack is the quartz/tungsten stack:

- `solid_quartz` is defined at `geometry/showermax/showerMaxDetectorSystem.gdml:19`.
- `solid_tungsten` is defined at `geometry/showermax/showerMaxDetectorSystem.gdml:27`.
- `solid_stack_region` is defined at `geometry/showermax/showerMaxDetectorSystem.gdml:160`.
- Individual detectors are placed as `singledet_01` through `singledet_28` in `showerMaxMother` at `geometry/showermax/showerMaxDetectorSystem.gdml:11204`.
- Detector placements run through `geometry/showermax/showerMaxDetectorSystem.gdml:11207` to `geometry/showermax/showerMaxDetectorSystem.gdml:11346`.

There are 28 ShowerMax detector modules.

## ShowerMax Virtual Planes

Important correction: the main ShowerMax virtual plane detector ID is `30`, not `70`.

The single legacy/global ShowerMax virtual plane is in:

- `geometry/mollerParallel.gdml:1047`
- It uses `midVirtualPlane_solid`.
- Its `DetNo` is `30` at `geometry/mollerParallel.gdml:1051`.
- Its position reference is `showerMaxVirtualPlane_pos` at `geometry/mollerParallel.gdml:1815`.
- That position is defined in `geometry/positions.xml:75` to `geometry/positions.xml:76` as `z="23920.00-90.00"`.

Older local notes may mention ShowerMax virtual plane `70`. Treat that as stale unless Justin explicitly says to use it. The current main ShowerMax virtual plane is `30`.

## Individual ShowerMax Stack Virtual Planes

There are also 28 individual ShowerMax stack virtual planes, one per detector module. These were added in the parallel geometry so each ShowerMax module can be studied separately instead of relying only on the single broad plane.

They are in:

- `geometry/mollerParallel.gdml`

How they were made:

- Two z constants were added for the staggered detector stack front locations:
  - `showerMaxStackFrontPlaneZLow` at `geometry/mollerParallel.gdml:13`
  - `showerMaxStackFrontPlaneZHigh` at `geometry/mollerParallel.gdml:14`
- A stack-plane box solid was added:
  - `Showermax_stack_VirtualPlane_solid` at `geometry/mollerParallel.gdml:109`
  - size: `x="181.61" y="314.69" z="0.5"` mm
- 28 logical virtual-plane volumes were added:
  - `showerMaxStackVirtualPlane01_log` through `showerMaxStackVirtualPlane28_log`
  - detector IDs run from `70030` through `72730`
  - definitions start at `geometry/mollerParallel.gdml:545`
  - final definition reaches `geometry/mollerParallel.gdml:740`
- 28 physical placements were added:
  - `showerMaxStackVirtualPlane01_phys` through `showerMaxStackVirtualPlane28_phys`
  - placements start at `geometry/mollerParallel.gdml:1591`
  - placements end at `geometry/mollerParallel.gdml:1756`

The individual planes follow the ShowerMax detector ring positions and rotations, with the same alternating low/high z staggering as the real detector modules.

Important mapping note from the file comment:

- `geometry/mollerParallel.gdml:545` to `geometry/mollerParallel.gdml:546` says the stack planes follow the legacy proxy ID order:
  - `70030` at `singledet_15`
  - `71430` at `singledet_01`
  - `72130` at `singledet_08`

So do not assume `70030` means `singledet_01`. Check the placement and mapping before interpreting detector-by-detector plots.

## Rate To Signal Reminder

At the virtual-plane level, the current proxy is:

```text
hit contribution = rate
```

For real ShowerMax signal, the needed form is:

```text
hit contribution = rate * Response(pid, energy, detector_id, hit_position, direction)
```

The response should become PE or PE-equivalent signal. Until that response model is supplied, `sum(rate)` at the virtual planes is only a rate-weighted hit proxy, not the final ShowerMax detector signal.

## Practical Analysis Conventions

- Use `hit.det` to select virtual planes.
- Use `hit.pid`, hit momentum/energy, hit position, and detector ID for response studies.
- Use one `rate` weight per selected virtual-plane hit for current ShowerMax hit-rate precision studies.
- Do not multiply by unrelated particle-array entries unless reproducing an older compatibility macro on purpose.
- For independent ROOT files generated with the same settings, average physical rate estimates between files rather than summing them as if they were different processes.
- For different physics processes, compare or combine their weighted signals after each process has been normalized through its own remoll `rate`.

## Existing Local Scripts To Know

- `_local_additions_archive/showermax-gem-analysis/scripts/calculate_rate_neff.C`
  - Older compatibility-style Neff macro, made to match `simple_showermax_xy.C` behavior.
- `_local_additions_archive/showermax-gem-analysis/scripts/calculate_production_precision.C`
  - Cleaner production precision macro using one `rate` per selected ShowerMax virtual-plane hit.
- `_local_additions_archive/showermax-gem-analysis/scripts/run_30_showermax_gem_simulations.sh`
  - Simple back-to-back remoll runner script. Do not run it unless Justin asks.

## Before Answering Or Editing

If asked about geometry, first inspect the relevant GDML/XML lines. If asked about ROOT analysis, first inspect the exact macro in question. If asked to change something, keep the edit tightly scoped and state exactly which file will be touched before touching it.
