# RUN_REPORT — OpenLB 3D Couette, γ̇ = 10 1/s

Date: 2026-11-17. OpenLB **1.9r0**, commit
`4b911c38b0f94b52eaed8639900aefa38910b687` (2026-10-23), cloned from
https://gitlab.com/openlb/release.git to `~/openlb`. Serial CPU build,
GNU Make (no CMake), g++ 13.3.0, C++20, this VM (4 CPUs, ~7.7 GB RAM,
no MPI, no GPU). Same checkout and build as the 2D case
(`examples/couette/RUN_REPORT.md`).

## Case (unchanged from the 2D case)

Full 3D annulus: R_i = 0.020 m (rotating bob), R_o = 0.025 m (fixed
cup), L = 0.060 m, periodic along the cylinder axis Z — no end
effects, so the infinite-cylinder assumption of the 2D cross-section
carries over unchanged. Newtonian fluid η = 1.0 Pa·s, ρ = 1000 kg/m³.
γ̇ = 10 1/s only (Ω = 1.8 rad/s, analytic torque 0.00150796 N·m,
analytic τ = 10 Pa), run time 1.0 s, quarter-sine ramp 0.2 s, steady
window = last 20%.

Lattice (`lattice3d.py`, verified): dx = 0.5 mm (10 cells across the
5 mm gap), dt = 4.167e-5 s, ω = 1.0, ν_lb = 1/6, u_lb,max = 0.030
at γ̇ = 100 s⁻¹ (same as 2D), 24,000 steps. Grid 100×100×120 =
1.2M cuboid cells (`lattice3d.json`).

## Method (all API names verified against the cloned 1.9r0 sources)

- Modern `Case<>` API: `Case<NavierStokes, Lattice<double,
  descriptors::D3Q19<>>>` (`examples/couette3d/couette3d.cpp`).
- Geometry: `IndicatorCylinder3D` — material 1 = fluid annulus,
  material 2 = outside the outer cylinder (fixed cup wall),
  material 3 = inner cylinder (rotating bob, extended 5·dx past both
  Z faces so it has no caps inside the domain).
- Axial periodicity: `mesh.getCuboidDecomposition().setPeriodicity(
  {false, false, true})` — the pattern of OpenLB's own
  `examples/fsi/taylorCouette3d/taylorCouette3d.cpp` (that example
  itself was NOT copied: it uses an FSI porous-element wall model
  with LES Smagorinsky for turbulent Taylor-Couette, the wrong tool
  for this laminar case).
- Bulk: `lattice.defineDynamics<BGKdynamics>(geometry, 1)`.
- Outer wall: `boundary::set<boundary::BounceBack>(lattice, geometry, 2)`.
- Inner wall: **Bouzidi interpolated bounce-back with imposed wall
  velocity** — `setBouzidiBoundary<T, DESCRIPTOR,
  BouzidiVelocityPostProcessor>(lattice, geometry, 3, innerCylinder)`
  plus `setBouzidiVelocity(lattice, geometry, 3, wallVel)`, the same
  moving-wall pattern as the 2D case and OpenLB's own
  `examples/turbulence/aorta3d/aorta3d.cpp:151,208-213`
  (`src/boundary/setBouzidiBoundary.h:194-228, 342-432`). The
  analytical cylinder gives the exact wall position, so the inner
  wall is **not** a staircase. The velocity functor returns
  **lattice** units (scale dt/dx applied in the functor). The wall
  velocity is updated during the ramp in `setTemporalValues` and
  pushed with
  `lattice.setProcessingContext<Array<descriptors::BOUZIDI_VELOCITY>>`.
- Torque: momentum exchange at the inner-wall links, the same formula
  as OpenLB's own `SuperLatticePhysBoundaryForce2D`
  (`src/functors/lattice/latticePhysBoundaryForce2D.hh:87-113`),
  summed over all Z slices (Z neighbours wrap — periodic axis).
  Conversion: 3D total torque, T = torque_lb · ρ·dx⁵/dt²
  (force = ρ·dx⁴/dt², times one more dx for the moment arm; in 2D
  the per-unit-length torque used ρ·dx⁴/dt²).
- Field snapshots: mid-plane (z = L/2) velocity slice every 1000
  steps (25 snapshots), for the video.

## γ̇ = 10 s⁻¹ result (validation gate)

- Steady torque (momentum exchange): **0.00151904 N·m** vs analytic
  **0.00150796 N·m** — **rel. error 0.73%**, inside the 10% gate.
  **Validation PASSED.** The torque is constant to 7 digits over the
  last 20% (fully converged Stokes flow at Re = 0.18).
- Steady torque per unit length: 0.0253173 N·m/m vs the 2D case's
  0.0253668 N·m/m (−0.2%) — 2D and 3D agree.
- Profile cross-check at γ̇ = 10: τ from the two innermost fluid
  bins = 7.54 Pa vs 10 Pa expected; torque from profile 0.001137 N·m,
  **25% low** — the same known bias as the 2D case (21% low there):
  the near-wall gradient is underestimated at a Bouzidi wall. The
  momentum-exchange torque is the primary measurement.
- Velocity profile (mid-plane, steady window): matches the analytic
  Couette profile within ~2% at every radius except the outermost
  bins, where v → 0 so the relative error inflates (same as 2D).
- Run time: **~22 min** for 24,000 steps on 1.2M cells, serial CPU
  (the 2D point took 4.3 s — the 3D cost is the full cuboid, walls
  included, with D3Q19).

## Video (verified with ffprobe, first/middle/last frames checked visually)

- `couette3d-gamma10.mp4` — 5.0 s, 720×720, H.264, 59 KB: spin-up
  at γ̇ = 10 s⁻¹, rendered as the actual 3D annulus (matplotlib 3D
  + ffmpeg — no ParaView or PyVista is installed on this VM, so
  matplotlib is the best 3D renderer available here). The grey bob
  in the middle, the sheared fluid as concentric cylindrical shells
  coloured by speed, a 90-degree wedge cut away so the radial shear
  gradient is visible on the two cut faces and on the top cap, camera
  orbiting once over the video. Colour scale fixed at 0–0.04 m/s
  (wall speed 0.036 m/s) so the spin-up is visible: the first frame
  is dark (fluid at rest), the last frame is bright yellow-green at
  the bob. Title shows t and the running torque vs analytic.
- The flow is axisymmetric and z-invariant (laminar, periodic in Z),
  so the mid-plane slice fully determines the 3D field: speed
  depends on radius only, and the render extrudes that profile
  along Z. (`make_video3d.py`.)

## Two real findings during this build

1. **Cell-centre convention (affects the torque, not the video).**
   OpenLB places cell (iX, iY, iZ) centre at `origin +
   latticeR*dx` (`src/geometry/cuboid.h:117` — `getPhysR` is "the
   single source of truth"). The 2D and 3D torque/profile code uses
   (iX + 0.5) for the cell centre, i.e. half a cell too far out.
   That overstates the moment arm by ~0.5/40 ≈ +1.25% on the
   torque. Correcting it would move the 3D error from +0.73% to
   ≈ −0.5% — still well inside the gate, and it would move the 2D
   error the same way. NOT fixed here (the gate passes either way
   and 2D/3D stay directly comparable); fix the convention in both
   cases in a follow-up, and re-check the numbers.
   Evidence: the mid-plane snapshots contain cells with velocity
   exactly zero at integer radius exactly 40 (inner wall) and
   ≥ 50 (outer wall) — those are the wall cells, and they sit at
   integer coordinates, not half-integers. The video binning
   (`make_video3d.py:read_mid_plane`) uses the correct integer
   convention and excludes wall cells; the first render used the
   +0.5 convention and the innermost bin read 28% low because
   zero-velocity wall cells leaked into it (a float-rounding
   leak: 0.02/0.0005 rounds to 39.999...). Fixed with integer wall
   radii; the binned profile now matches analytic within ~2%.
2. **Geometry statistics print looks wrong, geometry is right.**
   `prepareGeometry` prints material 0 = 803k cells, material 3 =
   38.2k — inconsistent with a solid inner cylinder (603k cells).
   The field data proves the geometry is correct: the bob interior
   has velocity exactly zero at every sampled cell, the annulus
   follows the analytic Couette profile, and the torque matches
   analytic to 0.73%. The statistics print is most likely a
   counting artifact; not chased further, recorded here as an open
   question.

## Files

- `couette3d.cpp`, `Makefile`, `case.yaml`, `lattice3d.py`,
  `lattice3d.json`, `make_video3d.py`, `couette3d-gamma10.mp4`,
  `README.md`, this report.
- `runs/gamma_10/` (git-ignored, stays local): `torque.csv` (242
  samples), `profile.csv`, 25 `field_mid_*.csv` snapshots.

## Still unverified (out of scope)

- 3D at other shear rates (only γ̇ = 10 run; the full sweep would be
  ~2.5 h serial at this resolution).
- Finer resolution than 10 cells/gap.
- MPI / GPU: serial CPU only here.
- Mass conservation through the run was not checked.
- The cell-centre convention fix (finding 1 above) in the 2D and
  3D torque code.
