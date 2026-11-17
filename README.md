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

## The Couette example — three-component fluid

Same geometry and sweep as the Newtonian Couette example, but the
fluid is the three-component (TC) model: τ = σ_y + σ_y·√(γ̇/γ̇_c) +
η_bg·γ̇, σ_y = 70 Pa, γ̇_c = 100 1/s, η_bg = 0.05 Pa·s. Implemented
as a local relaxation frequency via OpenLB's own generic
non-Newtonian hook (`visco::OmegaFromCell<collision::BGK,
tc::ThreeComponent>` — the Carreau-Yasuda/Casson pattern), with the
viscosity capped at η_max ≈ 7.02e4 Pa·s (γ̇_reg = 1e-3 1/s) because
the TC viscosity diverges as γ̇ → 0. The sweep is validated as
torque vs the analytic TC torque at each Ω (never η_app vs the TC
law at the label rate — the true inner-wall rate exceeds the label).
See `examples/couette-tc/`.

**TC validation PASSED: max analytic torque error 1.14%**
(tolerance 10%), at γ̇ = 0.1 where the yielded layer is only ~2
cells wide; other points 0.11–0.65%. At γ̇ ≤ 10 there is a genuine
unyielded plug attached to the outer wall (R_y = 21.03 → 24.84 mm);
γ̇ = 30 and 100 are fully yielded. Density drift ≤ 0.03% at
γ̇ ≤ 30, +0.34% at γ̇ = 100 — no leak of the kind XLB's halfway
bounce-back showed (−10.7%). Full sweep 74.8 s on the re-run.
Profile-derived torque reads 6–11% low (same direction as the
Newtonian 21% low); momentum exchange stays primary.

## The lid-driven cavity — three-component fluid

Square cavity L = 0.1 m, lid at U = 1.0 m/s, same TC fluid — built
to make a video of the shear bands / yielding, and to run on the
home PC (see below). The fluid yields near the lid (τ ≈ 93 Pa >
σ_y) while the slow corners stay unyielded. See
`examples/cavity-tc/`.

**Newtonian validation PASSED: Ghia Re = 100 centreline relative
L2 error 2.84%** (same binary, `--newtonian`). **TC run
converged** (mean-KE drift 0.00009%): steady plug (η > 1000 Pa·s)
= 38.7% of cells, capped = 25.1%; 455,325 steps, 1315 s on this
VM's CPU. Video `cavity-tc.mp4` shows speed, log-viscosity with
the η = 1000 Pa·s yield surface, and vorticity. 2D only, by
choice — 3D was not attempted.

## Solver comparison

Same Couette case, same validation, five solvers. Full write-up: [docs/comparison.md](docs/comparison.md).

| Solver | Method | Repo | Validation (max error) | η_app | Sweep video |
|---|---|---|---|---|---|
| M-Star | LBM, GPU, licensed | [mstar-skill](https://github.com/marcocaggioni/mstar-skill) (private) | not yet run | — | — |
| Palabos | LBM, C++/MPI, AGPLv3 | [palabos-skill](https://github.com/marcocaggioni/palabos-skill) | 1.43% | 1.0143 | [video](https://github.com/marcocaggioni/palabos-skill/blob/main/examples/couette/couette-sweep.mp4) |
| XLB | LBM, Python/JAX | [xlb-skill](https://github.com/marcocaggioni/xlb-skill) | 1.16% | 0.988–1.000 | [video](https://github.com/marcocaggioni/xlb-skill/blob/main/examples/couette/couette-sweep.mp4) |
| OpenFOAM | Finite-volume, GPL | [openfoam-skill](https://github.com/marcocaggioni/openfoam-skill) | 0.03% | 1.0003 | [video](https://github.com/marcocaggioni/openfoam-skill/blob/main/examples/couette/couette-sweep.mp4) |
| OpenLB | LBM, C++/CUDA, GPLv2 | [openlb-skill](https://github.com/marcocaggioni/openlb-skill) | 0.94% | 1.0093 | [video](https://github.com/marcocaggioni/openlb-skill/blob/main/examples/couette/couette-sweep.mp4) |

## Status

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
- 3D Couette (`examples/couette3d/`): **compiled and run** at γ̇ = 10 s⁻¹
  only — torque 0.00151904 N·m vs analytic 0.00150796 (+0.73%),
  validation PASSED; 100×100×120 grid, 24,000 steps, ~22 min serial.
  Video: [couette3d-gamma10.mp4](https://github.com/marcocaggioni/openlb-skill/blob/main/examples/couette3d/couette3d-gamma10.mp4)
  (5 s, 720×720, orbiting camera, wedge cutaway). See
  `examples/couette3d/RUN_REPORT.md`.
- Couette-TC (`examples/couette-tc/`): **compiled and run** against
  OpenLB 1.9r0 on 2026-11-17. TC rheology as a local relaxation
  frequency (OpenLB's `visco::OmegaFromCell` hook, capped viscosity);
  per-point dt from η at the label rate, ω = 1.0, 24,000 steps per
  point. **Validation PASSED: max analytic TC torque error 1.14%**,
  density drift ≤ 0.03% (γ̇ ≤ 30) / +0.34% (γ̇ = 100). Sweep 74.8 s.
  See `examples/couette-tc/RUN_REPORT.md`, `flow_curve.csv/png`,
  `couette-tc-gamma10.mp4` / `couette-tc-sweep.mp4`.
- Lid-driven cavity, TC fluid (`examples/cavity-tc/`): **compiled
  and run** on 2026-11-17. Newtonian Re = 100 check vs Ghia et al.
  (1982) **PASSED at 2.84%**; TC run converged (KE drift 0.00009%),
  plug 38.7%, 1315 s CPU. Video: [cavity-tc.mp4](https://github.com/marcocaggioni/openlb-skill/blob/main/examples/cavity-tc/cavity-tc.mp4).
  RTX 3090 runtime estimate ≈ 10–60 s of solver time — see
  `examples/cavity-tc/RUN_REPORT.md` for the basis and caveats.
- Still open: the 3D Couette sweep (γ̇ = 10 only so far), finer
  resolution, MPI/GPU (serial CPU only here — the GPU run is the
  `matteo` machine's job, see below).

## Files

- `skills/openlb-simulation/SKILL.md` — the skill, mirrored at `.github/skills/openlb-simulation/SKILL.md`
- `docs/openlb-guide.md` — install, build, lattice units, non-Newtonian notes, limits
- `examples/couette/` — the worked example
- `examples/couette-tc/` — the same Couette case with the TC fluid (validated 1.14%)
- `examples/cavity-tc/` — lid-driven cavity with the TC fluid (Ghia-validated setup, video)
- `scripts/check-setup.sh` — readiness check (copilot, uv, g++, make, OpenLB)

OpenLB: https://www.openlb.net — https://gitlab.com/openlb/release (GPLv2).
This repo's own code: MIT.
## Running on the Windows RTX 3090 machine

The serial CPU cases in this repo also run on Marco's home PC `matteo`
(Windows 11, NVIDIA GeForce RTX 3090 24 GB, driver 610.88, on
Tailscale) — via WSL2 Ubuntu. Nothing in the v0.1 cases needs the GPU
or CUDA: the serial CPU build needs only a C++20 compiler and GNU
Make (OpenLB's own documented requirements: "Compiler: C++20
compliant compiler (recent versions of GCC, Clang, Intel ICX).
Build System: GNU Make. GPU (Optional): NVIDIA CUDA 12.4+").
The 3090 matters for the bigger runs (finer grids, 3D, the cavity
video at higher resolution) — see the GPU notes at the end.

### WSL2 Ubuntu setup (one time)

```powershell
# PowerShell (Admin) on matteo — if WSL was never enabled:
wsl --install -d Ubuntu
# then reboot if asked, and create the Ubuntu user account.
# If Ubuntu is already installed but stopped:
wsl -d Ubuntu
```

Inside Ubuntu:

```bash
sudo apt update && sudo apt install -y g++ make git python3 python3-pip
# Python side (same as this VM): numpy, matplotlib, pyyaml
pip install numpy matplotlib pyyaml   # or: pip install uv && uv sync (repo root)
git clone https://gitlab.com/openlb/release.git ~/openlb
git clone https://github.com/marcocaggioni/openlb-skill.git ~/openlb-skill
cd ~/openlb-skill/examples/couette-tc
make            # plain make, one job; first build compiles the OpenLB core
python3 lattice.py
python3 couette_sweep.py --run
python3 couette_post.py
python3 make_video.py
```

What is already known about `matteo` (verified Nov 17, 2026): WSL2
is present with Ubuntu installed (currently stopped, startable); git
2.39.2 is present on the Windows side; no OpenLB directory exists
yet; no CMake and no nvcc (CUDA toolkit) are in the Windows PATH —
neither is needed for the serial CPU build above. The Windows-side
gaps do not block the WSL path: the compiler and Make live inside
WSL, installed by the `apt` line.

**To confirm on the machine** (not verifiable from here): that
`wsl -d Ubuntu` starts cleanly, that `g++ --version` and `make
--version` answer inside WSL after the `apt install`, and that the
first `make` in `examples/couette-tc/` completes (a few minutes for
the OpenLB core on first build). If `apt` has no network inside WSL,
check the WSL DNS / mirrored-networking setting before anything else.

### GPU (CUDA) build under WSL — only for the bigger runs

- The Windows 610.88 driver already exposes the GPU to WSL2 (NVIDIA's
  CUDA-on-WSL passthrough: no GPU driver install inside WSL — the
  driver comes from Windows). **To confirm on the machine:** run
  `nvidia-smi` inside WSL; it should list the RTX 3090.
- For a GPU build you additionally need the CUDA toolkit *inside*
  WSL (for `nvcc`): install CUDA 12.4+ per NVIDIA's WSL guide
  (`cuda-toolkit` from NVIDIA's Ubuntu repo — do NOT install a
  driver inside WSL). Then build OpenLB with the GPU platform:
  in `~/openlb/config.mk` set `PLATFORMS := CPU_SISD GPU_CUDA` (see
  `config/gpu_only.mk` for the template) and `make` the case the
  same way. The case code does not change — OpenLB's platform layer
  moves the same lattice to the GPU.
- Caveat: the TC cases use OpenLB's `visco::OmegaFromCell` path,
  which is marked non-vectorizable on CPU; its GPU throughput has
  NOT been measured — expect it to be slower than a pure-BGK GPU
  kernel, and measure it on the first GPU run rather than assuming
  the textbook MLUPS figure. The cavity runtime estimate in
  `examples/cavity-tc/RUN_REPORT.md` (measured CPU MLUPS here,
  extrapolated) is the basis for the first GPU run's expectations.

### Native Windows build (fallback, NOT the recommended path)

If WSL turns out unusable, a native build needs, at minimum: Visual
Studio 2022 with the C++ workload (MSVC, C++20), CMake (not in PATH
today — install it), the CUDA toolkit 12.4+ (for a GPU build; for a
CPU-only build, MSVC + CMake suffice), and a Git clone of OpenLB
built with CMake (OpenLB's `CMakeLists.txt` route) instead of GNU
Make. That is a strictly bigger install than the WSL route above
and has not been tried for this repo — prefer WSL.
