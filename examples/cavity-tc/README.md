# Lid-driven cavity — OpenLB, three-component (TC) fluid

Square cavity, top lid moving, other walls fixed. Same TC fluid and
eta-cap regularization as `examples/couette-tc/`. Built to make a
video of the shear bands / yielding — the fluid yields near the lid
while the slow corners stay unyielded — and to run on the home PC
(`matteo`, WSL Ubuntu, RTX 3090); verified here on CPU first.

## What is here

| File | What it does |
|---|---|
| `case.yaml` | The contract: geometry, TC fluid, lid speed, lattice. |
| `cavity.py` | Physical → lattice parameters + Ghia reference. **Tested here.** |
| `cavity.cpp` + `Makefile` | OpenLB case (2D cavity, GNU Make build). Same TC dynamics as Couette-TC. |
| `cavity_run.py` | Runs the Newtonian validation first, then the TC run. Dry-run by default. |
| `cavity_post.py` | Ghia check + steady-state + plug statistics. **Tested here.** |
| `make_video.py` | Renders the 3-panel video (speed, viscosity, vorticity). |

## Setup

L = 0.1 m square, U_lid = 1.0 m/s, N = 128 intervals per side
(dx = 0.781 mm, 129×129 lattice points). TC fluid σ_y = 70 Pa,
γ̇_c = 100 1/s, η_bg = 0.05 Pa·s, ρ = 1000 kg/m³, same η_max cap
(7.02e4 Pa·s at γ̇_reg = 1e-3 1/s) as Couette-TC. Characteristic
rate U_lid/L = 10 1/s → τ = 92.6 Pa at the lid, above the 70 Pa
yield stress (yielded), while the slow bottom corners stay
unyielded (plug). Effective Re ≈ 10.8. Run time 5.0 s (≈ 4.6
diffusion times, ≈ 50 lid turnovers), snapshots every 0.25 s.

## Run it

```bash
bash ../../scripts/check-setup.sh
python3 cavity.py            # lattice params + Ghia reference — must PASS
make                         # builds the case (OpenLB core already built once)
python3 cavity_run.py --run  # Newtonian Re=100 check, then the TC run
python3 cavity_post.py       # Ghia check + steady-state stats
python3 make_video.py        # cavity-tc.mp4 + cavity-steady.png
```

On `matteo` (WSL Ubuntu, RTX 3090) the same commands run as-is —
see "Running on the Windows RTX 3090 machine" in the repo README.
The v0.1 case is serial CPU and needs no CUDA; the GPU notes and
the RTX 3090 runtime estimate are in `RUN_REPORT.md`.

## Status

- `cavity.cpp`: **compiled and run** against OpenLB 1.9r0
  (commit 4b911c3, 2026-11-17) — see `RUN_REPORT.md`. Geometry and
  lid boundary follow OpenLB's own `laminar/cavity2d` example
  (one fix was needed: the lid rename needs overlap ≥ 2, and the
  mesh over [0, L] has N+1 lattice points per side).
- **Newtonian validation PASSED: Ghia Re = 100 centreline relative
  L2 error 2.84%** (tolerance 10%), 147,456 steps, ~160 s.
- **TC run converged**: mean-KE drift 0.00009% between the 80–90%
  and 90–100% windows; 455,325 steps, 1315 s on this VM's CPU
  (5.76 MLUPS). Steady plug (η > 1000 Pa·s) = 38.7% of fluid
  cells, capped (η = η_max) = 25.1%. Video `cavity-tc.mp4`
  (3.1 s, 1500×500) and `cavity-steady.png` in this folder.
  The steady field shows the yielded lid layer, a rotating
  unyielded plug at the vortex core, and the stagnant bottom plug.
- 3D: NOT attempted — 2D was chosen deliberately (the task: 3D
  only if cheap; a 3D TC cavity at this resolution would be
  ~100x the cost for a video that 2D already delivers).

## Results (steady state, t = 5 s)

| Quantity | Newtonian (Re = 100) | TC |
|---|---|---|
| Ghia centreline rel. L2 | 2.84% (PASSED) | — (no analytic reference) |
| Mean KE (m²/s²) | 0.03272 | 0.01012 |
| Max speed (m/s) | 0.952 | 0.926 |
| Energy drift (80–90% vs 90–100%) | 0.000% | 0.00009% |
| Density drift over the run | −0.98% | −0.25% |
| Plug, η > 1000 Pa·s | — | 38.7% of cells |
| Capped, η = η_max | — | 25.1% of cells |

Density drifts are slow and roughly linear in time (−0.20%/s
Newtonian, −0.05%/s TC) through the InterpolatedVelocity walls;
they do not disturb the steady velocity field (KE drift ~0,
Ghia 2.84%).
