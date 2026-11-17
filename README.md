# openlb-skill

Skill teaching an agent to run open-source OpenLB lattice Boltzmann
simulations from an English discussion — starting with a Couette
rheometer flow curve. **Public repo. No license needed.**

OpenLB (Karlsruhe Institute of Technology) is a free LBM framework
for general-purpose CFD, GPLv2, C++20, with 130+ shipped example
cases and MPI/OpenMP/CUDA/HIP parallelism. This is the open-source
counterpart of the private `mstar-skill` repo: same workflow, same
Couette example, different solver — no license seat, no GUI, no
`.msb` file. Instead you write a small C++ case and compile it with
GNU Make. It is the fourth solver in the series after Palabos, XLB
and OpenFOAM, and it ships a real Taylor-Couette example
(`examples/fsi/taylorCouette3d`) that this repo's 2D case is adapted
from.

## Quickstart

```bash
git clone https://github.com/marcocaggioni/openlb-skill.git
cd openlb-skill
uv sync              # Python side only (numpy, matplotlib, pyyaml)
bash scripts/install-openlb.sh   # OpenLB C++ — no PyPI package exists for it
bash scripts/check-setup.sh
copilot
```

First prompt:

```
/openlb-simulation run the Couette example in examples/couette/case.yaml
```

Or one go: `copilot -i "Use the /openlb-simulation skill to run the Couette example in examples/couette/case.yaml"`. Copilot must be started from the repo root; if the skill was added mid-session, run `/skills reload`. See `QUICKSTART.md`.

## The workflow

1. **Discuss** in English — geometry, fluid + rheology, sweep, what is measured.
2. **Write `case.yaml`**, show it, get a yes. Never skip this, never invent a number.
3. **Lattice parameters** — `lattice.py` converts physical → lattice units and checks the low-Mach condition. This is the step M-Star hides.
4. **Build** the C++ case (GNU Make). Dry-run by default.
5. **Sweep / run** — one run per shear rate, torque CSV per run.
6. **Post-process** torque → flow curve, with a mandatory Newtonian validation gate (analytic Couette torque, ~10%).
7. **Report** — plot + CSV + assumptions + validation + unverified items.

## The Couette example

Fluid in a Couette rheometer, inner cylinder rotating, outer fixed,
torque on the bob → flow curve. Defaults (labelled, to confirm):
R_i = 0.020 m, R_o = 0.025 m, L = 0.060 m, Newtonian η = 1.0 Pa·s,
ρ = 1000 kg/m³, γ̇ = 0.1 → 100 s⁻¹ log-spaced, 10 cells across the
gap. 2D cross-section first (torque per unit length); 3D is a
follow-up. Same numbers as the M-Star example, so the solvers can
be compared point for point.

## Solver comparison

Same Couette case, same validation, five solvers. Full write-up: [docs/comparison.md](docs/comparison.md).

| Solver | Method | Repo | Validation (max error) | η_app | Sweep video |
|---|---|---|---|---|---|
| M-Star | LBM, GPU, licensed | [mstar-skill](https://github.com/marcocaggioni/mstar-skill) (private) | not yet run | — | — |
| Palabos | LBM, C++/MPI, AGPLv3 | [palabos-skill](https://github.com/marcocaggioni/palabos-skill) | 1.43% | 1.0143 | [video](https://github.com/marcocaggioni/palabos-skill/blob/main/examples/couette/couette-sweep.mp4) |
| XLB | LBM, Python/JAX | [xlb-skill](https://github.com/marcocaggioni/xlb-skill) | 1.16% | 0.988–1.000 | [video](https://github.com/marcocaggioni/xlb-skill/blob/main/examples/couette/couette-sweep.mp4) |
| OpenFOAM | Finite-volume, GPL | [openfoam-skill](https://github.com/marcocaggioni/openfoam-skill) | 0.03% | 1.0003 | [video](https://github.com/marcocaggioni/openfoam-skill/blob/main/examples/couette/couette-sweep.mp4) |
| OpenLB | LBM, C++/CUDA, GPLv2 | [openlb-skill](https://github.com/marcocaggioni/openlb-skill) | 0.94% | 1.0093 | [video](https://github.com/marcocaggioni/openlb-skill/blob/main/examples/couette/couette-sweep.mp4) |

## Status — v0.1

- `lattice.py` and `couette_post.py`: **tested in this VM** (lattice
  check passes; post-processing recovers η ≈ 1.0 Pa·s on synthetic
  torque, validation PASSED).
- `couette.cpp`: **compiled and run** against OpenLB 1.9r0
  (commit 4b911c3) on 2026-11-17. Moving wall = Bouzidi
  interpolated bounce-back with imposed velocity
  (`setBouzidiBoundary<BouzidiVelocityPostProcessor>` +
  `setBouzidiVelocity`, aorta3d pattern); torque = momentum
  exchange at inner-wall links. Full sweep validated:
  η_app = 1.0093 Pa·s (true 1.0), analytic torque error 0.94%
  max — see `examples/couette/RUN_REPORT.md`,
  `flow_curve.csv/png`, and the videos
  `couette-gamma10.mp4` / `couette-sweep.mp4`. Sweep took 32.5 s.
- Still open: non-Newtonian fluids and 3D (2D cross-section here).

## Files

- `skills/openlb-simulation/SKILL.md` — the skill, mirrored at `.github/skills/openlb-simulation/SKILL.md`
- `docs/openlb-guide.md` — install, build, lattice units, non-Newtonian notes, limits
- `examples/couette/` — the worked example
- `scripts/check-setup.sh` — readiness check (copilot, uv, g++, make, OpenLB)

OpenLB: https://www.openlb.net — https://gitlab.com/openlb/release (GPLv2).
This repo's own code: MIT.
