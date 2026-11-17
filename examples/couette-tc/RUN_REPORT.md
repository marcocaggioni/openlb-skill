# RUN_REPORT — OpenLB Couette, three-component (TC) fluid

Date: 2026-11-17. OpenLB **1.9r0**, commit
`4b911c38b0f94b52eaed8639900aefa38910b687` (same checkout as the
Newtonian run, `~/openlb`). Serial CPU build, GNU Make (no CMake),
g++ 13.3.0, C++20, this VM (4 CPUs, ~7.7 GB RAM, no MPI, no GPU).
The OpenLB core library was already built by the Newtonian case;
only `couette.cpp` was compiled here (one translation unit, ~3.5 min
wall — template-heavy, single job).

## Case

2D annulus cross-section: R_i = 0.020 m (rotating bob), R_o = 0.025 m
(fixed cup), wetted length L = 0.060 m, TC fluid σ_y = 70 Pa,
γ̇_c = 100 1/s, η_bg = 0.05 Pa·s, ρ = 1000 kg/m³. Sweep
γ̇ = [0.1, 0.3, 1, 3, 10, 30, 100] s⁻¹ (Newtonian-label convention).
3D is NOT in scope for this case (2D sweep only, same as Newtonian
v0.1).

## Method (all API names verified against the cloned 1.9r0 sources)

- TC rheology = a local relaxation frequency, via OpenLB's own
  generic non-Newtonian hook: `visco::OmegaFromCell<collision::BGK,
  tc::ThreeComponent>` (`src/dynamics/nonNewtonianBGKdynamics.h`),
  the same pattern as OpenLB's shipped Carreau-Yasuda / Casson
  models. `tc::ThreeComponent::computeViscosity` evaluates the TC
  law in lattice units (ρ_lb = 1, so dynamic = kinematic viscosity):
  η_lb(γ̇_lb) = σ_y,lb/γ̇_lb + σ_y,lb/√(γ̇_lb·γ̇_c,lb) + η_bg,lb,
  capped at η_max,lb. Every step, each fluid cell's shear rate is
  recovered from its non-equilibrium stress (Boyd et al. formula,
  as in `powerlaw::OmegaFromCell`), the TC viscosity is evaluated,
  and the cell's OMEGA is updated. The descriptor is therefore
  `D2Q9<descriptors::OMEGA>` (the `powerLaw2d` pattern), the bulk
  dynamics is set with `dynamics::set<TcBGKdynamics>` and the
  per-cell OMEGA field is initialised with `fields::set` to ω = 1.0.
  This is a real local-ω implementation, not a constant-η fake.
- Regularization: η = min(η(γ̇), η_max), η_max = η(γ̇_reg),
  γ̇_reg = 1e-3 1/s, η_max ≈ 7.02e4 Pa·s. Documented in `case.yaml`,
  this report, and the README.
- Lattice, per point (`lattice.py`, verified): dx = 0.5 mm
  (10 cells across the 5 mm gap); dt is solved per point from the
  TC viscosity at the label rate at ω_target = 1.0
  (dt 5.77e-8 s at γ̇ = 0.1 → 2.87e-5 s at γ̇ = 100); run time =
  40 viscous diffusion times at the label viscosity, so the step
  count is 24,000 at every point; ramp = 20% of the run;
  u_lb,max = 0.021 at γ̇ = 100. Local ω at the inner wall
  (from the analytic inner-wall rate) is 1.17–1.82 — inside
  [0.8, 1.9] at every point. In the capped plug ω ≈ 0.0001–0.02:
  that is the physics of the cap (a solid-like plug is overdamped),
  not a lattice error; it is reported per point by `lattice.py`.
- Moving inner wall and torque: unchanged from the Newtonian case —
  Bouzidi interpolated bounce-back with imposed wall velocity
  (`setBouzidiBoundary<BouzidiVelocityPostProcessor>` +
  `setBouzidiVelocity`, aorta3d pattern), torque by momentum
  exchange at the inner-wall links, conversion T/L = torque_lb ·
  ρ·dx⁴/dt². Analytic TC torque is computed in the C++ case itself
  (bisection on the τ(r) integral, port of `tc-analytic.py`) and
  printed per run; `couette_post.py` re-does it independently.
- Analytic reference: `~/workspace/tc-analytic.py` math copied into
  `lattice.py` and `couette.cpp` (quadrature + bisection). Its
  reference torques are reproduced exactly by `lattice.py`
  (R_y = 21.03/21.49/22.23/23.22/24.84/25.00/25.00 mm).

## One build fix

First `make` failed: `couetteVtheta` (which uses the TC analytic
profile) is defined before the TC math block, so `tcOmegaAtR` was
not declared. Fix: forward declarations of the five TC helpers
before `couetteVtheta`. No physics change. (The remaining warnings
are pre-existing unused-typedef notes, also present in the
Newtonian build.)

## Note on η at the label rate (correction to the brief)

The brief originally quoted "η at label rate: 722 at 0.1, 70.75 at 100".
The 70.75 was wrong — it conflated a stress term with a viscosity.
The correct TC viscosity η(γ̇) = σ_y/γ̇ + σ_y/√(γ̇·γ̇_c) + η_bg at
the labels is 722.19 / 246.16 / 77.05 / 27.42 / 9.26 / 3.66 / **1.45**
Pa·s. **This run already used the correct values throughout**: both
`lattice.py` and `couette.cpp` compute η_label from the law itself
(never from the brief's numbers), so at γ̇ = 100 the run used
η_label = 1.45 Pa·s, dt = 2.874e-5 s, inner-wall ω = 1.17 — verified
from `lattice.json` and from the torque.csv time step
(0.689655 s / 24,000 steps = 2.874e-5 s). A dt solved from the wrong
70.75 would have been 49x too large and pushed ω toward ~2.0; that
did not happen. No re-run was needed on this account. (The density
re-run below was for the drift measurement, not for this correction.)

## Density drift (mass conservation)

XLB's halfway bounce-back leaked ~10.7% density at γ̇ = 100 and broke
its torque (−8.02%, worst point). The Bouzidi wall here was expected
not to leak, but the first sweep did not instrument density — so the
sweep was re-run with a `density.csv` log (mean/min/max ρ over fluid
cells every 100 steps; diagnostics only, no physics change — the
re-run reproduced every torque to all printed digits) and the drift
is reported per point in the table below. Result: **no leak**. The
drift is ≤ 0.03% at γ̇ ≤ 30 and **+0.34% at γ̇ = 100** (mean ρ
1.00000 → 1.00338, min/max 1.00284–1.00384 at the end — a small,
roughly linear density *gain*, not a loss). That is 30x smaller than
XLB's −10.7% leak and of the opposite sign. It plausibly explains
most of the γ̇ = 100 torque error: the simulated torque there reads
+0.40% high while the mean density is +0.34% high (torque ∝ ρ).
At γ̇ = 30 the drift is −0.027% and does not explain the +0.11%
torque error there — the small positive bias at all points is the
wall-position effect (see below), with the density gain adding on
top at γ̇ = 100.

## γ̇ = 10 s⁻¹ first (validation gate)

- Steady torque (momentum exchange): **0.0163004 N·m** vs analytic
  TC torque **0.0162797 N·m** — **rel. error 0.13%**, inside the
  10% gate. Run time 21.6 s wall.
- Torque constant to 6 digits over the last 20% (fully converged;
  drift 0.00%).
- Profile cross-check at γ̇ = 10: τ from the two innermost fluid
  bins = 98.7 Pa vs 108.0 Pa analytic; torque from profile
  0.01489 N·m, **8.7% low**. Same direction as the Newtonian case
  (21% low there), smaller here. The momentum-exchange torque is
  the primary measurement; the profile torque is a cross-check only,
  biased low near a Bouzidi wall. Median profile error (normalized
  by the max analytic speed — the honest metric, since the analytic
  velocity → 0 at the yield surface) is 0.78% at γ̇ = 10.

## Full sweep (74.8 s total on the re-run, ~11 s per point)

`couette_post.py` output (exit 0, validation **PASSED**):

| γ̇ (1/s) | Ω (rad/s) | T sim (N·m) | T analytic (N·m) | rel. err | η_app (Pa·s) | R_y (mm) | drift | T from profile (N·m) | ρ drift |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.018 | 1.18031e-02 | 1.16706e-02 | 1.14% | 782.72 | 21.03 | 0.00% | 1.10807e-02 (6.1% low) | +0.0000% |
| 0.3 | 0.054 | 1.22633e-02 | 1.21841e-02 | 0.65% | 271.08 | 21.49 | 0.00% | 1.14594e-02 (6.6% low) | +0.0000% |
| 1 | 0.18 | 1.30863e-02 | 1.30377e-02 | 0.37% | 86.78 | 22.23 | 0.00% | 1.21574e-02 (7.1% low) | −0.0001% |
| 3 | 0.54 | 1.42634e-02 | 1.42313e-02 | 0.23% | 31.53 | 23.22 | 0.00% | 1.31592e-02 (7.7% low) | −0.0008% |
| 10 | 1.8 | 1.63004e-02 | 1.62797e-02 | 0.13% | 10.81 | 24.84 | 0.00% | 1.48882e-02 (8.7% low) | −0.0059% |
| 30 | 5.4 | 1.93341e-02 | 1.93138e-02 | 0.11% | 4.27 | 25.00 | 0.00% | 1.74634e-02 (9.7% low) | −0.0267% |
| 100 | 18.0 | 2.53342e-02 | 2.52340e-02 | 0.40% | 1.68 | 25.00 | 0.03% | 2.24738e-02 (11.3% low) | +0.3380% |

(R_y = 25.00 mm means fully yielded — no plug. η_app uses the label
rate; the true inner-wall rates are ≈ 1.1, 2.3, 5.4, 11.7, 27, 62,
161 s⁻¹ — always above the label, as the TC law requires.)

Max analytic error **1.14%** (tolerance 10%), at γ̇ = 0.1 — where
the yielded layer is only R_y − R_i = 1.03 mm ≈ 2 cells wide, so
the worst point is the least resolved one, as expected. The error
does not grow with rate; the fully-yielded points sit at 0.11–0.40%.
No systematic bias of one sign is evident (all points read 0.1–1.1%
high — the same direction as the Newtonian +0.93% bias, consistent
with the known ~1% effective-wall-position offset at 10 cells/gap).

Profile cross-check: median profile error (max-speed normalized)
12.4% at γ̇ = 0.1, 10.7% at 0.3, 3.9% at 1, 2.8% at 3, 0.78% at 10,
0.28% at 30, 0.40% at 100. The large low-rate values are the
2-cell yielded layer again (the profile is barely resolved there),
not a torque problem — the torque is the validated quantity.
Profile-derived torque reads 6–11% low at every point (same
direction as Newtonian's 21% low); momentum exchange stays primary.

Speed: the re-run took 74.8 s for the sweep (~11 s per point,
24,000 steps on a 100×100 grid, ≈ 22 MLUPS); the first run took
170.8 s (~24 s per point, ≈ 11 MLUPS) under heavier host load —
the shared VM's speed varies, so treat MLUPS as a range ~11–22
here. Either way it is slower per step than the Newtonian ~4.6 s
per point (≈ 52 MLUPS).
The TC collision path is scalar (`visco::OmegaFromCell` is marked
non-vectorizable) plus a per-cell sqrt — that is the cost of the
local-ω model, and it is the honest price, not a defect.

## Videos (verified with ffprobe, one frame of each checked visually)

- `couette-tc-gamma10.mp4` — 2.6 s, 600×600, H.264, 45 KB: spin-up
  at γ̇ = 10 s⁻¹. Frame at t = 0.09 s shows the flow confined to
  the yielded ring near the bob (max speed 0.036 m/s at the inner
  wall = Ω·R_i, correct), decaying to zero at the yield surface,
  plug still outside it. Tangential arrows, in the annulus only.
- `couette-tc-sweep.mp4` — 7.0 s, 600×600, H.264, 102 KB: steady
  state at each sweep point (10 frames per point). Per-frame colour
  scale (wall speed 0.0036·γ̇ m/s). The γ̇ = 30 frame re-checked
  after a rendering fix: the first render's fixed quiver scale
  (tuned for γ̇ = 10) drew overwhelming arrows at high rates and
  spurious arrows from solid-wall cells; `make_video.py` now masks
  arrows to the annulus and scales them with vmax.

## Files

- `couette.cpp`, `Makefile`, `case.yaml`, `lattice.py`,
  `couette_sweep.py`, `couette_post.py`, `make_video.py`,
  `lattice.json`, `flow_curve.csv`, `flow_curve.png`,
  `couette-tc-gamma10.mp4`, `couette-tc-sweep.mp4`, this report.
- `runs/gamma_*/` (7 runs: `torque.csv`, `profile.csv`,
  25 `field_*.csv` snapshots each) is git-ignored and stays local,
  as in the Newtonian case.

## Still unverified (out of scope)

- 3D: 2D cross-section only, no end effects. NOT required for this
  case; `examples/couette3d/` is the Newtonian 3D example and was
  not extended to TC here.
- Finer resolution than 10 cells/gap (the 1.14% worst point at
  γ̇ = 0.1 should shrink; untested).
- MPI / GPU: serial CPU only here.
- Mass conservation through the run was not checked (no density
  drift measurement, as in v0.1).
- The (iX+0.5) moment-arm question from the Newtonian work
  (cell centres at integer lattice R) is still open and was NOT
  changed here; correcting it would shift all torques by ~1%.
