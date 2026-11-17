# OpenLB — short guide for this repo

Sources: https://www.openlb.net, https://gitlab.com/openlb/release,
the user guide in `doc/` of the checkout. Anything not in those
sources is marked as unverified in this repo — do not treat it as fact.

## What OpenLB is

Open-source lattice Boltzmann framework from the Karlsruhe Institute
of Technology (KIT), GPLv2, native interface in C++20. 130+ example
cases ship with it, including a real Taylor-Couette example
(`examples/fsi/taylorCouette3d`, 3D, FSI wall model — the pattern
this repo's 2D case is adapted from). Parallelism via MPI, OpenMP,
CUDA and HIP. It is a **library**, not an application: you write a
C++ program against it and compile it with GNU Make. There is no
GUI, no case file format, no built-in post-processor — output is
VTK (ParaView) plus whatever your program writes.

## Install / first build

```bash
git clone https://gitlab.com/openlb/release.git
cd openlb/examples/laminar/cavity2d
make
```

Requires a C++20 compiler (gcc/clang) and GNU Make — those are the
only mandatory dependencies. MPI is optional (parallel runs); this
repo's v0.1 path is serial. Use plain `make` (one job) on small
machines: parallel builds of the template-heavy core can exhaust
memory. The first build compiles the core library (`libolbcore.a`)
and the embedded dependencies (zlib, tinyxml2) — a few minutes, once.

## The modern Case API (v1.9)

Examples follow one pattern (`examples/forBeginners/cylinder2d` is
the cleanest to copy):

```cpp
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D2Q9<>>>;
// Mesh from an IndicatorCuboid2D + physDeltaX
// Geometry: SuperGeometry material numbers via geometry.rename(...)
//   with IndicatorCircle2D / IndicatorCuboid2D
// Lattice: lattice.defineDynamics<BGKdynamics>(geometry, 1)
//   boundary::set<boundary::BounceBack>(lattice, geometry, 2)
//   setBouzidiBoundary<T,DESCRIPTOR,BouzidiVelocityPostProcessor>(
//     lattice, geometry, 3, innerCircle)   // moving wall, aorta3d pattern
//   lattice.setParameter<descriptors::OMEGA>(omega)
// Velocity on a moving Bouzidi wall: setBouzidiVelocity(lattice,
//   geometry, mat, functor) — the functor returns LATTICE velocity
//   (evaluated at physical wall-intersection points); update it
//   during a run from setTemporalValues, before collideAndStream.
//   Do NOT use boundary::set<boundary::InterpolatedVelocity> on a
//   solid obstacle: it iterates over every material cell including
//   the interior, where the discrete normal is (0,0), and throws
//   "Could not set Boundary".
// Time loop: lattice.collideAndStream() per step.
```

Unit conversion: `UnitConverter(dx, dt, charL, charU, physViscosity,
physDensity)` — or the resolution-based converters used by the
examples. `converter.getLatticeRelaxationFrequency()` should equal
your omega; check it at runtime. Force conversion is
`rho * dx^4 / dt^2` (converter.getPhysForce).

## Lattice units — the one thing to get right

OpenLB works in lattice units internally. Convert physical → lattice
before setting anything, and convert forces/torques back afterwards:

```
dx = characteristic length / resolution
dt from omega:  omega = 1 / (3*nu_lb + 0.5),  nu_lb = nu * dt / dx^2
u_lb = u_phys * dt / dx          keep < 0.1
rho_lb = 1.0 (by convention)
force conversion = rho * dx^4 / dt^2
```

`examples/couette/lattice.py` does this for the Couette case.
Choosing `omega` near 1 is a good default; near 2 the scheme goes
unstable, near 0.5 it is over-damped.

## Boundaries

- No-slip walls: `boundary::set<boundary::BounceBack>`.
- Velocity boundaries on flat inlets/outlets:
  `boundary::set<boundary::InterpolatedVelocity>` plus
  `momenta::setVelocity` with a functor (constant or spatially
  varying). This works on thin boundary layers, NOT on a solid
  obstacle — on a solid disc it throws "Could not set Boundary"
  (the discrete normal is (0,0) in the interior).
- Curved / moving walls: the shipped `setBouzidiBoundary` (used for
  the cylinder in `forBeginners/cylinder2d`); for a moving wall use
  the `BouzidiVelocityPostProcessor` operator plus
  `setBouzidiVelocity` with a functor returning lattice velocity
  (the `turbulence/aorta3d` pattern) — this is what this repo's
  Couette case uses, and the analytical indicator gives the exact
  wall position instead of a staircase.

## Torque

There is no built-in torque output. The tool OpenLB provides is
`SuperLatticePhysBoundaryForce2D` — pointwise physical force on a
boundary material by momentum exchange
(`src/functors/lattice/latticePhysBoundaryForce2D.hh`). This repo
evaluates the same formula directly from the block-lattice cells
and sums r x F over the inner-wall nodes; a profile-derived torque
is computed as a cross-check.

## Non-Newtonian fluids

There is no rheology menu. The pattern to copy is OpenLB's own
`powerLaw2d` example (and `nonNewtonianPoiseuille3d`): a local
relaxation frequency from the local shear rate. Herschel–Bulkley
needs the same plus a yield stress and a viscosity cap at low shear
rate. Implement Newtonian first, validate it, then add the model —
never the other way round.

## Limits (v0.1)

- The Couette case is 2D; the shipped Taylor-Couette example is 3D
  with a turbulent wall model and FSI machinery — deliberately not
  copied, to keep v0.1 laminar and simple.
- Post-processing is your own code (this repo's `couette_post.py`)
  or ParaView on the VTK output.
- Serial CPU only in v0.1; MPI/CUDA/HIP untested here (OpenLB
  supports all three — that is a genuine advantage over Palabos,
  whose GPU support is pre-release).
- OpenLB is GPLv2: fine for a public repo, check before using it
  inside anything you would not publish.
