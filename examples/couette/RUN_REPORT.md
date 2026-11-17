# RUN_REPORT — OpenLB Couette, Newtonian validation

Date: 2026-11-17. OpenLB **1.9r0**, commit
`4b911c38b0f94b52eaed8639900aefa38910b687` (2026-10-23, "examples:
Restore the vertical alignment of the parameter set calls"),
cloned from https://gitlab.com/openlb/release.git to `~/openlb`.
Serial CPU build, GNU Make (no CMake), g++ 13.3.0, C++20, this VM
(4 CPUs, ~7.7 GB RAM, no MPI, no GPU). The OpenLB core library
itself was proven first with the shipped `laminar/cavity2d` example.

## Case (unchanged from the other solver repos)

2D annulus cross-section: R_i = 0.020 m (rotating bob), R_o = 0.025 m
(fixed cup), wetted length L = 0.060 m (only used by the post script),
Newtonian fluid η = 1.0 Pa·s, ρ = 1000 kg/m³. Sweep
γ̇ = [0.1, 0.3, 1, 3, 10, 30, 100] s⁻¹, run time 1.0 s per point,
quarter-sine ramp 0.2 s, steady window = last 20%.

Lattice (`lattice.py`, verified): dx = 0.5 mm (10 cells across the
5 mm gap), dt = 4.167e-5 s, ω = 1.0, ν_lb = 1/6, u_lb,max = 0.030
at γ̇ = 100 s⁻¹, 24,000 steps per point. `lattice.json` is written by
`lattice.py`.

## Method (all API names verified against the cloned 1.9r0 sources)

- Modern `Case<>` API: `Case<NavierStokes, Lattice<double,
  descriptors::D2Q9<>>>` (`examples/couette/couette.cpp`).
- Geometry: `IndicatorCircle2D` — material 1 = fluid annulus,
  material 2 = outside the outer circle (fixed cup wall),
  material 3 = inner disc (rotating bob).
- Bulk: `lattice.defineDynamics<BGKdynamics>(geometry, 1)`.
- Outer wall: `boundary::set<boundary::BounceBack>(lattice, geometry, 2)`
  (halfway bounce-back on the staircase outer wall).
- Inner wall: **Bouzidi interpolated bounce-back with imposed wall
  velocity** — `setBouzidiBoundary<T, DESCRIPTOR,
  BouzidiVelocityPostProcessor>(lattice, geometry, 3, innerCircle)`
  plus `setBouzidiVelocity(lattice, geometry, 3, wallVel)`,
  the moving-wall pattern of OpenLB's own
  `examples/turbulence/aorta3d/aorta3d.cpp:151,208-213`
  (`src/boundary/setBouzidiBoundary.h:194-228, 342-432`).
  The analytical circle gives the exact wall position, so the
  inner wall is **not** a staircase. The velocity functor returns
  **lattice** units (`setBouzidiVelocity` expects lattice units
  although it evaluates the functor at physical intersection
  points — the aorta3d pattern); the scale dt/dx is applied in the
  functor. The wall velocity is updated during the ramp in
  `setTemporalValues`, before `collideAndStream`, and pushed with
  `lattice.setProcessingContext<Array<descriptors::BOUZIDI_VELOCITY>>`.
- Torque: momentum exchange at the inner-wall links, the same
  formula as OpenLB's own `SuperLatticePhysBoundaryForce2D`
  (`src/functors/lattice/latticePhysBoundaryForce2D.hh:87-113`):
  per inner-wall cell and per link to a fluid neighbour,
  F -= c_i · (f_fluid[iPop] + f_wall[opp(iPop)]), summed over the
  wall, torque = Σ r × F, evaluated directly from the block-lattice
  cells (unambiguous in serial). Sign-flipped to the driving torque
  (the raw sum is the fluid drag opposing rotation). Conversion:
  2D torque per unit length, T/L = torque_lb · ρ·dx⁴/dt²
  (= `converter.getPhysForce`), **not** ρ·dx⁵/dt² (that would be a
  3D torque). A second, independent torque from the simulated
  azimuthal velocity profile is printed as a cross-check (see below).

## Two real bugs found and fixed during this build

1. **Functor assignment is deleted** (`couette.cpp:419`, first build):
   `wallVel = RotatingWallVelocity<double>(cxPhys, cyPhys)` failed
   because `olb::GenericF` deletes copy/move assignment
   (`src/functors/genericF.h:63`). Fix: construct the global functor
   once at (0,0) and add `setCenter()` / `setScale()` methods; the
   centre and the physical→lattice scale are set in `main()` after
   argv parsing.
2. **`boundary::set<boundary::InterpolatedVelocity>` cannot be used
   on a solid disc** (first run, runtime error `Could not set
   Boundary`, `src/boundary/setBoundary2D.h:127`): `boundary::set`
   iterates over *every* material-3 cell, including the disc
   interior, where the discrete normal is (0,0) and the Flat-case
   dynamics construction throws. It is meant for thin boundary
   layers (e.g. cylinder2d's flat inlet strip), not a solid obstacle.
   Fix: the Bouzidi moving-wall pattern above, which acts on the
   fluid cells adjacent to the solid and takes the analytical wall
   position. This also removed the staircase error at the inner wall.

## γ̇ = 10 s⁻¹ first (validation gate)

- Steady torque (momentum exchange): **0.00152201 N·m**
  (0.0253668 N·m/m × 0.060 m) vs analytic **0.00150796 N·m** —
  **rel. error 0.93%**, inside the 10% gate. Run time 4.3 s.
- Steady torque is constant to 7 digits over the last 20% (fully
  converged Stokes flow at Re = 0.18).
- Profile cross-check at γ̇ = 10: τ from the two innermost fluid
  bins = 7.89 Pa vs 10 Pa expected; torque from profile 0.00119 N·m,
  **21% low**. The innermost bin velocity is 2.6% low vs analytic
  while the next bin is 0.4% low, so the near-wall gradient is
  underestimated. The momentum-exchange torque is the primary
  measurement (it is exact at the wall links); the profile torque is
  a cross-check only, and it is biased low near a Bouzidi wall.
  Velocity profile overall: median bin error 3.25% vs analytic
  (worst bins are the outermost, r = 0.0248, where v → 0 so the
  relative error inflates, and the innermost bin r = 0.0198, whose
  centre lies just inside the bob).

## Full sweep (32.5 s total, ~4.6 s per point)

`couette_post.py` output (exit 0, validation **PASSED**):

| γ̇ (1/s) | Ω (rad/s) | T sim (N·m) | T analytic (N·m) | rel. err | η_app (Pa·s) |
|---|---|---|---|---|---|
| 0.1 | 0.018 | 1.52201e-05 | 1.50796e-05 | 0.93% | 1.0093 |
| 0.3 | 0.054 | 4.56603e-05 | 4.52389e-05 | 0.93% | 1.0093 |
| 1 | 0.18 | 1.52201e-04 | 1.50796e-04 | 0.93% | 1.0093 |
| 3 | 0.54 | 4.56602e-04 | 4.52389e-04 | 0.93% | 1.0093 |
| 10 | 1.8 | 1.52201e-03 | 1.50796e-03 | 0.93% | 1.0093 |
| 30 | 5.4 | 4.56603e-03 | 4.52389e-03 | 0.93% | 1.0093 |
| 100 | 18.0 | 1.52210e-02 | 1.50796e-02 | 0.94% | 1.0094 |

Max analytic error **0.94%** (tolerance 10%). η_app = 1.0093 Pa·s
at every point (true 1.0) — the error is constant across the sweep,
i.e. a systematic bias, not noise. Likely cause: the effective wall
position / Bouzidi interpolation at 10 cells across the gap
(the outer wall is still a staircase bounce-back). Finer resolution
was not tested; it should shrink the bias, as the Palabos run's
constant 1.43% staircase bias did.

## Videos (verified with ffprobe, one frame of each checked visually)

- `couette-gamma10.mp4` — 2.5 s, 600×600, H.264, 51 KB: spin-up at
  γ̇ = 10 s⁻¹. Frame at t = 0.83 s shows the rotating annulus, max
  speed 0.036 m/s at the inner wall (= Ω·R_i, correct), decaying to
  zero at the outer wall, tangential arrows.
- `couette-sweep.mp4` — 7.0 s, 600×600, H.264, 108 KB: steady state
  at each sweep point (10 frames per point). Fixed-colour-scale
  first render left the low-rate frames invisible, so `make_video.py`
  now scales each frame to that point's wall speed (0.0036·γ̇ m/s);
  the γ̇ = 3 frame re-checked after the fix shows the annulus
  correctly at ~0.0108 m/s max.

## Files

- `couette.cpp`, `Makefile`, `case.yaml`, `lattice.py`,
  `couette_sweep.py`, `couette_post.py`, `make_video.py`,
  `lattice.json`, `flow_curve.csv`, `flow_curve.png`,
  `couette-gamma10.mp4`, `couette-sweep.mp4`, this report.
- `runs/gamma_*/` (7 runs, 6.7 MB: `torque.csv`, `profile.csv`,
  25 `field_*.csv` snapshots each) is git-ignored and stays local,
  as in the Palabos repo.

## Still unverified (out of scope for v0.1)

- Non-Newtonian fluids: no model implemented yet (Newtonian only).
  OpenLB ships `powerLaw2d` and `nonNewtonianPoiseuille3d` as the
  patterns to copy; a local-relaxation-frequency implementation is
  the planned route.
- 3D: 2D cross-section only, no end effects. OpenLB ships a real
  3D Taylor-Couette example (`examples/fsi/taylorCouette3d`) that
  was deliberately not copied for v0.1.
- Finer resolution than 10 cells/gap (the 0.93% systematic bias
  should shrink; untested).
- MPI / GPU: serial CPU only here. OpenLB supports MPI/OpenMP/
  CUDA/HIP, none of it exercised.
- Mass conservation through the run was not checked (no density
  drift measurement in v0.1).
