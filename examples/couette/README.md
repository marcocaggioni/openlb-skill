# Couette example — OpenLB

Fluid in a Couette rheometer: inner cylinder rotates, outer cylinder
fixed, torque on the bob gives the flow curve. Same geometry and
fluid as the M-Star example (`marcocaggioni/mstar-skill`), so the
solvers can be compared point for point.

## What is here

| File | What it does |
|---|---|
| `case.yaml` | The contract: geometry, fluid, sweep, lattice. Confirm before running. |
| `lattice.py` | Physical → lattice parameters (dx, dt, omega, Mach check). **Tested here.** |
| `couette.cpp` + `Makefile` | OpenLB case (2D annulus, GNU Make build). See status below. |
| `couette_sweep.py` | Runs the sweep, one run per shear rate. Dry-run by default. |
| `couette_post.py` | Torque → flow curve + Newtonian validation gate. **Tested on synthetic torque.** |
| `make_video.py` | Renders the spin-up and sweep videos from the field snapshots. |

## Defaults

R_i = 0.020 m, R_o = 0.025 m, L = 0.060 m, Newtonian η = 1.0 Pa·s,
ρ = 1000 kg/m³, γ̇ = 0.1 → 100 s⁻¹ log-spaced, 10 cells across the
5 mm gap, ω = 1.0, 1 s per point, steady window = last 20%.
2D cross-section first (torque per unit length); 3D is a follow-up.

## Run it

```bash
bash ../../scripts/check-setup.sh
uv run python lattice.py          # dx, dt, omega, Mach check — must PASS
make                              # builds the case (first build compiles OpenLB core)
uv run python couette_sweep.py --run
uv run python couette_post.py     # flow_curve.csv + flow_curve.png, validation gate
uv run python make_video.py       # couette-gamma10.mp4 + couette-sweep.mp4
```

## Status (v0.1)

- `lattice.py`: verified — dx 0.5 mm, dt 4.17e-5 s, ω = 1.0,
  u_lb,max = 0.030 at γ̇ = 100 s⁻¹, 24,000 steps per point.
- `couette_post.py`: verified on synthetic torque — recovers
  η_app ≈ 1.0 Pa·s, validation PASSED.
- `couette.cpp`: **compiled and run** against OpenLB 1.9r0
  (commit 4b911c3, 2026-11-17) — see `RUN_REPORT.md`. Moving inner
  wall = Bouzidi interpolated bounce-back with imposed velocity
  (`setBouzidiBoundary<BouzidiVelocityPostProcessor>` +
  `setBouzidiVelocity`, the aorta3d pattern), updated during the
  ramp; torque = momentum exchange at inner-wall links (the same
  formula as OpenLB's `SuperLatticePhysBoundaryForce2D`).
  **Newtonian validation PASSED: max analytic torque error 0.94%,
  η_app = 1.0093 Pa·s** (true 1.0) at all 7 sweep points; sweep
  32.5 s. Videos `couette-gamma10.mp4` / `couette-sweep.mp4` and
  `flow_curve.csv/png` in this folder.
- Still unverified: non-Newtonian fluids (no model implemented yet
  — Newtonian only), 3D (2D cross-section here), finer resolution,
  MPI/GPU.

## Analytic check

γ̇_i = 2·Ω·R_o²/(R_o²−R_i²), τ_i = T/(2π·R_i²·L),
T_analytic = 4π·L·η·Ω·R_i²·R_o²/(R_o²−R_i²).
At γ̇ = 10 s⁻¹: Ω = 1.8 rad/s (17.2 rpm), T = 1.508e-3 N·m, τ = 10 Pa.
