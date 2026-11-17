# Single deformable cell in shear flow — OpenLB immersed boundary

A single deformable capsule in shear flow, using OpenLB 1.9's
immersed-boundary model for resolved deformable particles. This is
OpenLB's shipped benchmark `fsi/deformableParticleInShearFlow3d`
(Krüger et al. 2011), copied here unchanged with a run/post wrapper:
the validation is the Taylor deformation parameter D vs normalized
time κt, against the analytic steady value.

## What is here

| File | What it does |
|---|---|
| `case.yaml` | The contract: geometry, membrane, shear, validation reference. |
| `rbc-shear.cpp` + `Makefile` | The shipped OpenLB example, copied unchanged (serial, no MPI). |
| `parametersMeshes.xml` | Membrane mesh + physics (Skalak kS = kAlpha = 0.0111111). |
| `parametersPositions.xml` | Cell position (centre of the 50³ domain). |
| `sph_ico_1280.msh` | Icosphere mesh (1280 faces, radius 5 lattice units). |
| `rbc_run.py` | Runs the case in `runs/single/`. Dry-run by default. |
| `rbc_post.py` | Node positions → D/G vs κt + validation gate. **Tested here.** |

## Scope — read this first

The "cell" is a **spherical capsule** (icosphere mesh), not a
biconcave red blood cell. That is what Krüger's benchmark
validates, and it is the honest first step: single cell, known
reference, before any suspension or RBC shape. The membrane is a
Skalak membrane (shear + area dilatation, kS = kAlpha = 0.0111111
lattice units) with **no** bending rigidity and **no** volume/area
conservation (kB = kV = kA = 0) — the benchmark's parameters, not a
physiological RBC model. A biconcave shape and RBC mechanics are
follow-ups, not claimed here.

## Setup

Domain 50×50×50 lattice units, shear in the x–y plane (top wall
+x, bottom wall −x, periodic in x and y). Cell: sphere radius 5,
centre (25, 25, 25). Fluid: Newtonian BGK with Guo forcing,
Re = 0.02, relaxation time 1.0. Dimensionless shear rate (capillary
number) G = η·γ̇·r/k_s = 0.01. Run length κt = 120 (9000 lattice
steps), 16 snapshots (every κt = 8).

## Run it

```bash
bash ../../scripts/check-setup.sh
make                         # builds the case (OpenLB core already built once)
python3 rbc_run.py --run     # runs in runs/single/
python3 rbc_post.py          # deformation.csv + deformation.png, validation gate
```

## Validation

D = (a−c)/(a+c) from the principal semi-axes (a ≥ b ≥ c) of the
volume inertia tensor — the same definition as OpenLB's own
`writeParticleInertia` (GSL path; no GSL is installed here, so
`rbc_post.py` mirrors that computation in numpy from the written
node positions — validated on a synthetic sphere, D = 0, and a
1.1/0.9 stretched ellipsoid, D = 0.100). Reference: Krüger,
Varnik & Raabe (2011), Comput. Math. Appl. 61(12), 3485–3505,
doi:10.1016/j.camwa.2010.03.057, and OpenLB's shipped validation
figure: D/G vs κt, analytic steady value D/G = 6.25 at κt = 120;
the shipped runs reach ~6.8 (Δx = 1) and ~6.4 (Δx = 0.5). Gate:
steady D/G (mean of the last 3 snapshots) within 15% of 6.25 —
the shipped Δx = 1 run is itself +9% off analytic, so a 10% gate
would fail a correct reproduction.

## Status

**PASSED, 2026-11-17** (see `RUN_REPORT.md`). Steady D/G = 7.019 vs
analytic 6.25 (+12.3%, inside the 15% gate); shipped Δx = 1 run
reaches ~6.8. Run: 359.6 s wall, 9000 steps, 50³ grid, 3.13 MLUPS.
Volume constant at 519.1 through the run. Next steps (not done):
finer grid (Δx = 0.5) for a tighter match, biconcave RBC shape,
small suspension.
