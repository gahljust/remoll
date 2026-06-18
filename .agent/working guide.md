# Agent Working Guide

Repo root: `/Users/justingahley/G4/remoll`.

## What To Treat As Fixed

- Treat `*.root` as read-only data. Do not edit, regenerate, or delete them unless explicitly asked.
- Do not change remoll source, geometry, macros, or build files unless the request is specifically about those files.
- Do not run remoll, long builds, or long tests unless asked.
- Use C++/ROOT `.C` macros for analysis work unless the user asks for another language.

## Current Layout

- Active local analysis work lives in `_local_additions_archive/current/`.
- Older snapshots are in `_local_additions_archive/old/` and should not be edited.
- Local analysis docs live in `_local_additions_archive/current/README.md`.
- Rate notes live in `.agent/rate information.md`.
- remoll rate background also exists in `guide/rateCalc.md`.

## Current Analysis Focus

- The main work is ShowerMax virtual-plane hit studies for `A_PV` extraction.
- Use plain `rate` for current displays and summaries.
- Treat `rate * energy` as a placeholder for a future PE/signal model, not the final signal.

## Useful Current Files

- `_local_additions_archive/current/make_showermax_acceptance_tree.C`
- `_local_additions_archive/current/sview.C`
- `_local_additions_archive/current/view_showermax_acceptance_tree.C`
- `_local_additions_archive/current/analysis_help/`
- `_local_additions_archive/current/neff_phase_space/`

## Practical Rules

- Keep edits local and narrowly scoped.
- Prefer references to exact files and lines before changing behavior.
- Do not assume older path names are still valid; check under `current/` first.
- If asked about rate weighting, use `.agent/rate information.md` as the source of truth.
- Make all responses very short and concise unless asked differently.
