# remoll Rate Information

This note records what we worked out about the `rate` branch in remoll ROOT files and how we are using it for the ShowerMax virtual-plane studies.

## Short Meaning

For normal remoll event generators such as `moller`, `elastic`, `inelastic`, `elasticAl`, `inelasticAl`, `elasticC12`, and `pion`, the ROOT `rate` branch is an event weight with units of approximately events/second. It converts a forced/generated Monte Carlo event into its physical contribution at the chosen beam current, target, material, and generator phase-space settings.

The key formula in remoll is:

```cpp
fEvent->fRate = fEvent->fEffXs * fBeamTarg.GetEffLumin(sampling_type) / nthrown;
```

Code reference:

- `src/remollPrimaryGeneratorAction.cc:205-212`
- `src/remollRunAction.cc:67-68`
- `src/remollIO.cc:167-172`

`nthrown` is the number of events requested for that remoll job. This means each event's `rate` is already normalized by the number of thrown events in that job.

## Dimensional Meaning

The physical structure is:

```text
rate = effective cross section x effective luminosity / number thrown
```

For target-based generators:

```text
effective luminosity = beam electrons per second x target scattering centers per area
```

The code path is:

```cpp
return fBeamCurrent / (e_SI*coulomb) * fEffectiveMaterialLength;
```

Code reference:

- `src/remollBeamTarget.cc:79-84`

`fBeamCurrent / electron charge` gives beam electrons per second.

`fEffectiveMaterialLength` is not a plain geometric length. It is an areal number density factor. The relevant code is:

```cpp
effective_length = z_half_length*2.0 * material->GetDensity();
fEffectiveMaterialLength =
  (total_effective_length/effective_length) *
  effective_length * Avogadro/masssum;
```

Code reference:

- `src/remollBeamTarget.cc:305-310`
- `src/remollBeamTarget.cc:368-386`

Dimensionally:

```text
2*z_half_length x density = mass / area
(mass / area) / molar_mass = mole / area
(mole / area) x Avogadro = scattering centers / area
```

Then:

```text
cross section [area] x beam electrons/s x scattering centers/area
= interactions/s
```

The factor:

```cpp
(total_effective_length/effective_length) * effective_length
```

algebraically reduces to `total_effective_length` for the selected volume, but written this way it exposes the intended Monte Carlo sampling weight times the sampled material thickness.

## Target Sampling

For normal generators, the sampling type defaults to active target volume:

- `src/remollVEventGen.cc:48`

The target volumes are discovered from GDML auxiliary tags:

- `src/remollDetectorConstruction.cc:640-708`
- `include/remollBeamTarget.hh:57-70`

For the LH2 production target macro:

```text
/remoll/target/mother LH2
/remoll/target/volume LH2
```

Reference:

- `macros/target/LH2.mac:1-3`

The current target ladder has LH2, upstream aluminum window, and downstream aluminum window sampling names:

- `geometry/target/targetLadder.gdml:58-75`
- `geometry/target/targetLadder.gdml:98`

For aluminum-window studies we used:

```text
/remoll/target/volume USAl
/remoll/target/volume DSAl
```

## Cross Sections

Geant4 does not decide the primary Moller/elastic/inelastic cross-section weight for these forced remoll generated events. remoll's event generators calculate an effective cross section, then Geant4 transports the particles.

Examples:

Moller:

- `src/remollGenMoller.cc:33-44`
- analytic QED-like Moller cross section, then `SetEffCrossSection(...)`

ep elastic:

- `src/remollGenpElastic.cc:187-223`
- `src/remollGenpElastic.cc:242`
- Mott-like cross section and form factors, then `SetEffCrossSection(...)`

ep inelastic:

- `src/remollGenpInelastic.cc:35-44`
- uses proton/neutron inelastic models, then `SetEffCrossSection(...)`

C12:

- `src/remollGenC12.cc:68-104`
- generator phase space and nuclear model, then `SetEffCrossSection(...)`

Al:

- `src/remollGenAl.cc:68-102`
- generator phase space and aluminum model, then `SetEffCrossSection(...)`

Pions:

- `src/remollGenPion.cc:74-124`
- Wiser parameterization, phase space, then `SetEffCrossSection(...)`

## What `sum(rate)` Means

If two simulations have the same generator, beam current, target, phase-space cuts, and geometry, then the event-level sum:

```text
sum over generated events of rate
```

should converge to the same physical generated-process rate as event count increases.

It should not grow forever with more events because each event's rate is divided by `nthrown`.

For a selected detector condition, such as ShowerMax virtual-plane hits:

```text
sum selected hits of rate
```

estimates the selected hit-rate proxy for that detector condition.

If an event has multiple selected hits, summing once per hit counts hit multiplicity. That is useful for virtual-plane hit-rate studies, but it is not the same as counting selected events.

## Rate Across Multiple ROOT Jobs

The remoll guide says not to simply add final rate estimates from separate jobs:

- `guide/rateCalc.md:32-40`

Reason: each job's `rate` values were normalized by that job's `Nthrown`. If each job is an independent run with the same settings and the same number of requested events, then each job estimates the same physical rate. The final physical rate estimate should be averaged:

```text
R_final = (R_1 + R_2 + ... + R_Njobs) / Njobs
```

But the statistical precision improves because independent samples add statistics:

```text
Neff_total = Neff_1 + Neff_2 + ... + Neff_Njobs
```

So five independent 20k jobs should have about the same precision as one 100k job, if all settings match and all jobs completed cleanly.

If a job did not complete the requested number of events, its rate normalization is wrong for the number of events actually present. The guide recommends discarding incomplete jobs:

- `guide/rateCalc.md:36-38`

## Neff Precision

For weighted samples:

```text
Neff = (sum w)^2 / sum(w^2)
precision fraction = 1 / sqrt(Neff)
precision percent = 100 / sqrt(Neff)
```

For our ShowerMax virtual-plane precision pilot, the intended weight is:

```text
w = rate
```

with one weight per selected ShowerMax virtual-plane hit.

Later, when detector response is available, the intended weight becomes:

```text
w = rate x ShowerMax response
```

## How We Used Rate For ShowerMax

For each production process `i`, the virtual-plane hit-rate proxy is:

```text
R_i^SM = sum over selected ShowerMax virtual-plane hits from process i of rate
```

With detector response:

```text
S_i^SM = sum over selected hits from process i of rate x response(hit)
```

Then the ShowerMax process fraction is:

```text
f_i^SM = S_i^SM / sum_j S_j^SM
```

Before response is available, we use:

```text
f_i^SM approx R_i^SM / sum_j R_j^SM
```

The response-weighted process asymmetry should be:

```text
<A_i>_SM = sum(rate x response x A_i) / sum(rate x response)
```

The mixed measured ShowerMax asymmetry is approximately:

```text
A_meas^SM = sum_i f_i^SM x <A_i>_SM
```

This is why `rate` is essential: it puts separately generated samples, such as Moller, ep elastic, ep inelastic, aluminum, pions, and pi0/gamma proxy, onto the same physical scale.

Equal generated event counts do not mean equal physical contribution.

## Important Warning About `calculate_rate_neff.C`

The earlier compatibility macro:

- `_local_additions_archive/showermax-gem-analysis/scripts/calculate_rate_neff.C`

was written to mimic `simple_showermax_xy.C` behavior. That compatibility calculation adds rate once per selected ShowerMax hit times every accepted `part.p` entry:

```text
selected_sum_rate = sum_events rate_event x N_selected_hits_event x N_good_part_entries_event
```

Relevant code:

- `_local_additions_archive/showermax-gem-analysis/scripts/calculate_rate_neff.C:185-194`
- `_local_additions_archive/plane_analysis/simple_showermax_xy.C:100-138`

That was useful for reproducing the old `simple_showermax_xy.C` Neff/Akin number, but it is not the clean official ShowerMax physical hit-rate definition.

For official production precision studies, use one `rate` weight per selected virtual-plane hit, as implemented in:

- `_local_additions_archive/showermax-gem-analysis/scripts/calculate_production_precision.C`

Relevant code:

- `_local_additions_archive/showermax-gem-analysis/scripts/calculate_production_precision.C:176-200`

## Process Names Used

remoll generator names:

```text
Moller:           /remoll/evgen/set moller
ep elastic:       /remoll/evgen/set elastic
ep inelastic:     /remoll/evgen/set inelastic
Al elastic:       /remoll/evgen/set elasticAl
Al quasi-elastic: /remoll/evgen/set quasielasticAl
Al inelastic:     /remoll/evgen/set inelasticAl
pion:             /remoll/evgen/set pion
```

Pion subtype examples:

```text
/remoll/evgen/set pion
/remoll/evgen/pion/settype pi-
```

```text
/remoll/evgen/set pion
/remoll/evgen/pion/settype pi0
```

Generator map reference:

- `src/remollPrimaryGeneratorAction.cc:50-64`

