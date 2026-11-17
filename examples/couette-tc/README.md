# Couette example — OpenLB, three-component (TC) fluid

Fluid in a Couette rheometer: inner cylinder rotates, outer cylinder
fixed, torque on the bob gives the flow curve. Same geometry and sweep
labels as the Newtonian example (`examples/couette/`), so the solvers
can be compared point for point — but the fluid is the TC model
(Marco, Nov 17, voice), not Newtonian.

## What is here

| File | What it does |
|---|---|
| `case.yaml` | The contract: geometry, TC fluid, sweep, lattice, regularization. |
| `lattice.py` | Physical → lattice parameters, PER SWEEP POINT (dx, dt, omega, Mach check, yield radius). **Tested here.** |
| `couette.cpp` + `Makefile` | OpenLB case (2D annulus, GNU Make build). TC rheology as a local relaxation frequency. See status below. |
| `couette_sweep.py` | Runs the sweep, one run per shear rate. Dry-run by default. |
| `couette_post.py` | Torque → flow curve + TC validation gate (10%). **Tested here.** |
| `make_video.py` | Renders the spin-up and sweep videos from the field snapshots. |

## The TC fluid

τ = σ_y + σ_y·√(γ̇/γ̇_c) + η_bg·γ̇, with σ_y = 70 Pa,
γ̇_c = 100 1/s, η_bg = 0.05 Pa·s (50 cP), ρ = 1000 kg/m³.
R_i = 0.020 m, R_o = 0.025 m, L = 0.060 m, γ̇ = 0.1 → 100 s⁻¹
log-spaced (labels use the Newtonian inner-wall-rate convention,
Ω = γ̇·(R_o²−R_i²)/(2·R_o²)). 2D cross-section, torque per unit
length; 10 cells across the 5 mm gap.

**Regularization (required).** The TC viscosity
η(γ̇) = σ_y/γ̇ + σ_y/√(γ̇·γ̇_c) + η_bg diverges as γ̇ → 0, so it is
capped: η = min(η(γ̇), η_max) with η_max = η(γ̇_reg),
γ̇_reg = 1e-3 1/s, i.e. η_max ≈ 7.02e4 Pa·s. The cap is the
unyielded plug: at γ̇ ≤ 10 the analytic solution has a genuine plug
attached to the outer wall (R_y < R_o; see the table below).

**Validation is torque vs analytic torque.** The sweep is validated
against the analytic TC torque at each Ω (bisection on the
τ(r) = T/(2πLr²) integral, `tc-analytic.py` math copied into
`lattice.py` and `couette.cpp`) — NOT by comparing η_app to the TC
law at the label rate, because the true inner-wall rate exceeds the
label rate (the yield stress carries most of the stress).

## Run it

```bash
bash ../../scripts/check-setup.sh
python3 lattice.py          # per-point dx, dt, omega, Mach check — must PASS
make                        # builds the case (OpenLB core already built once)
python3 couette_sweep.py --run
python3 couette_post.py     # flow_curve.csv + flow_curve.png, validation gate
python3 make_video.py       # couette-tc-gamma10.mp4 + couette-tc-sweep.mp4
```

## Status (v0.1)

- `lattice.py`: verified — per-point dt from η at the label rate at
  ω = 1.0 (dt 5.77e-8 s at γ̇ = 0.1 → 2.87e-5 s at γ̇ = 100),
  24,000 steps per point, u_lb,max = 0.021 at γ̇ = 100. Local ω at
  the inner wall 1.17–1.82 (in [0.8, 1.9]); in the capped plug
  ω ≈ 0.02 by design (overdamped solid-like plug, not an error).
- `couette.cpp`: **compiled and run** against OpenLB 1.9r0
  (commit 4b911c3, 2026-11-17) — see `RUN_REPORT.md`. TC rheology
  via OpenLB's own generic non-Newtonian hook
  (`visco::OmegaFromCell<collision::BGK, tc::ThreeComponent>`,
  the Carreau-Yasuda/Casson pattern in
  `src/dynamics/nonNewtonianBGKdynamics.h`): each cell's shear rate
  comes from its non-equilibrium stress, the TC viscosity (capped)
  sets the cell ω every step. Moving inner wall = Bouzidi
  interpolated bounce-back with imposed velocity (as in the
  Newtonian case); torque = momentum exchange at inner-wall links.
- **TC validation PASSED: max analytic torque error 1.14%**
  (tolerance 10%), at γ̇ = 0.1 where the yielded layer is only
  ~2 cells wide. Full sweep 74.8 s on the re-run (~11 s per point —
  the scalar TC collision path runs slower per step than Newtonian
  BGK, ~11–22 MLUPS here vs ~52; first run took 170.8 s under
  heavier host load). Density drift ≤ 0.03% at γ̇ ≤ 30, +0.34% at
  γ̇ = 100 — no XLB-style leak (see RUN_REPORT.md). Videos `couette-tc-gamma10.mp4` /
  `couette-tc-sweep.mp4` and `flow_curve.csv/png` in this folder.
  The γ̇ = 10 frame shows the physics clearly: flow confined to
  the yielded ring near the bob, the outer plug still.
- Still unverified: 3D (2D cross-section here, and 3D is NOT in
  scope for this case), finer resolution (the 1.14% worst point
  should shrink — the yielded layer is 2 cells wide there),
  MPI/GPU.

## Results

| γ̇ (1/s) | Ω (rad/s) | T sim (N·m) | T analytic (N·m) | rel. err | R_y (mm) | η_app (Pa·s) |
|---|---|---|---|---|---|---|
| 0.1 | 0.018 | 1.18031e-02 | 1.16706e-02 | 1.14% | 21.03 | 782.72 |
| 0.3 | 0.054 | 1.22633e-02 | 1.21841e-02 | 0.65% | 21.49 | 271.08 |
| 1 | 0.18 | 1.30863e-02 | 1.30377e-02 | 0.37% | 22.23 | 86.78 |
| 3 | 0.54 | 1.42634e-02 | 1.42313e-02 | 0.23% | 23.22 | 31.53 |
| 10 | 1.8 | 1.63004e-02 | 1.62797e-02 | 0.13% | 24.84 | 10.81 |
| 30 | 5.4 | 1.93341e-02 | 1.93138e-02 | 0.11% | 25.00 | 4.27 |
| 100 | 18.0 | 2.53342e-02 | 2.52340e-02 | 0.40% | 25.00 | 1.68 |

R_y = 25.00 mm means fully yielded (no plug). At γ̇ = 10:
Ω = 1.8 rad/s (17.2 rpm), T = 1.630e-2 N·m, τ_i = 108.1 Pa,
true inner-wall rate ≈ 27 s⁻¹ (vs the 10 s⁻¹ label).
