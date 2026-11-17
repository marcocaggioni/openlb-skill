# Run report — single deformable cell in shear flow (RBC test)

Date: 2026-11-17. Machine: this VM (4 CPUs, ~7.7 GB RAM, serial CPU).
OpenLB 1.9r0, commit `4b911c38`. Case: `examples/rbc-shear/`.

## What was run

OpenLB's shipped benchmark `fsi/deformableParticleInShearFlow3d`
(Krüger et al. 2011), copied unchanged into `examples/rbc-shear/`
(`rbc-shear.cpp`, `parametersMeshes.xml`, `parametersPositions.xml`,
`sph_ico_1280.msh`) with a run/post wrapper. Single spherical
capsule (icosphere, 1280 faces, radius 5 lattice units) at the centre
of a 50³ domain, shear in the x–y plane, Re = 0.02, G = 0.01,
9000 lattice steps to κt = 120, 16 snapshots (every κt = 8).

## Commands

```bash
cd examples/rbc-shear && make        # builds the case (first build ~142 s)
python3 rbc_run.py --run             # 359.6 s wall, 3.13 MLUPS average
python3 rbc_post.py                  # D/G vs κt, validation gate
```

## Result — PASSED

Steady D/G (mean of the last 3 snapshots, κt = 104–120) = **7.019**
vs analytic **6.25** → **+12.3%**, inside the 15% gate
(`rbc_post.py` exits 0). The shipped Δx = 1 run reaches ~6.8 (read
off the shipped figure); this run is +3% above that, within
figure-reading error. The curve shape matches the shipped figure:
steep rise, plateau from κt ≈ 60.

| κt | D | D/G | a | b | c | V |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.0000 | 0.000 | 4.981 | 4.981 | 4.981 | 519.1 |
| 8 | 0.0085 | 0.849 | 5.023 | 4.981 | 4.939 | 519.1 |
| 16 | 0.0307 | 3.068 | 5.135 | 4.983 | 4.830 | 519.1 |
| 24 | 0.0473 | 4.734 | 5.219 | 4.988 | 4.747 | 519.1 |
| 32 | 0.0573 | 5.735 | 5.268 | 4.995 | 4.696 | 519.1 |
| 40 | 0.0630 | 6.301 | 5.294 | 5.002 | 4.666 | 519.1 |
| 48 | 0.0662 | 6.618 | 5.308 | 5.008 | 4.649 | 519.1 |
| 56 | 0.0680 | 6.795 | 5.315 | 5.013 | 4.638 | 519.1 |
| 64 | 0.0689 | 6.894 | 5.318 | 5.016 | 4.632 | 519.1 |
| 72 | 0.0695 | 6.950 | 5.320 | 5.018 | 4.629 | 519.1 |
| 80 | 0.0698 | 6.982 | 5.321 | 5.019 | 4.627 | 519.1 |
| 88 | 0.0700 | 7.000 | 5.322 | 5.020 | 4.625 | 519.1 |
| 96 | 0.0701 | 7.011 | 5.322 | 5.021 | 4.625 | 519.1 |
| 104 | 0.0702 | 7.017 | 5.322 | 5.021 | 4.624 | 519.1 |
| 112 | 0.0702 | 7.020 | 5.322 | 5.021 | 4.624 | 519.1 |
| 120 | 0.0702 | 7.021 | 5.322 | 5.022 | 4.624 | 519.1 |

(Full table in `deformation.csv`; plot in `deformation.png`.)

## Checks and caveats

- **D definition**: Taylor D = (a−c)/(a+c) from the principal
  semi-axes of the volume inertia tensor — the same definition as
  OpenLB's own `writeParticleInertia` (GSL path). No GSL is
  installed in this VM, so `rbc_post.py` mirrors that computation
  in numpy from the written node positions (`.vtp` files). The
  post-processing was validated before the run on a synthetic
  sphere (D = 0.00000, V exact) and a 1.1/0.9 stretched ellipsoid
  (D = 0.1000 exactly). The run's own κt = 0 snapshot gives
  D = 0.0000, a = b = c = 4.981 — the icosphere's equivalent
  ellipsoid radius (faceted mesh, not exactly 5.000).
- **Volume**: constant at 519.1 through the whole run (no volume
  drift), even though kV = 0 in this benchmark.
- **Why +12.3% off analytic?** The shipped figure shows the same
  overshoot for both resolutions (+9% at Δx = 1, +2% at Δx = 0.5),
  so it is a property of the benchmark at G = 0.01, not a failure
  of this run. The 15% gate was chosen for exactly this reason
  (a 10% gate would fail a correct reproduction). A finer grid
  (Δx = 0.5, i.e. N = 100, r = 10) is the obvious next check if a
  tighter match is wanted — not run here (cost ~8× the 3D grid).
- **Scope**: spherical capsule, not a biconcave RBC; Skalak membrane
  with no bending and no volume/area conservation (kB = kV = kA = 0).
  That is the benchmark's setup. A biconcave RBC shape and
  physiological membrane parameters are follow-ups, not claimed.
- **Suspension**: not attempted — the single-cell test passes, so a
  small suspension is now a legitimate next step, not yet done.
- **GPU**: not planned for this test (CPU here first, per the brief).

## Runtime

Build ~142 s (first build; OpenLB core already built once).
Run 359.6 s wall for 9000 steps on a 50³ grid (3.13 MLUPS).
Post-processing <2 s.
