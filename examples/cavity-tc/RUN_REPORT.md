# RUN_REPORT — OpenLB lid-driven cavity, three-component (TC) fluid

Date: 2026-11-17. OpenLB **1.9r0**, commit
`4b911c38b0f94b52eaed8639900aefa38910b687` (same checkout as the
Couette runs, `~/openlb`). Serial CPU build, GNU Make (no CMake),
g++ 13.3.0, C++20, this VM (4 CPUs, ~7.7 GB RAM, no MPI, no GPU).

## Case

2D square cavity L = 0.100 m, top lid moving at U_lid = 1.0 m/s,
other walls fixed. TC fluid σ_y = 70 Pa, γ̇_c = 100 1/s,
η_bg = 0.05 Pa·s, ρ = 1000 kg/m³, same η_max cap as Couette-TC
(7.02e4 Pa·s at γ̇_reg = 1e-3 1/s). Characteristic rate U_lid/L =
10 1/s → τ_char = 92.6 Pa > σ_y (yielded at the lid), slow bottom
corners unyielded (plug). Effective Re = U·L/ν(η_char) ≈ 10.8.
Run time 5.0 s (≈ 4.6 diffusion times L²/ν_char = 1.08 s, ≈ 50 lid
turnovers L/U_lid = 0.1 s), snapshots every 0.25 s. 2D only, by
choice: 3D was not attempted (the task allowed skipping it if not
cheap; a 3D TC cavity at this resolution would cost ~100x for a
video 2D already delivers).

## Method

- TC dynamics: the SAME `tc::ThreeComponent` model and
  `visco::OmegaFromCell` hook as `examples/couette-tc/` (no new
  rheology code). Descriptor `D2Q9<descriptors::OMEGA>`.
- Lattice (`cavity.py`, verified): N = 128 intervals per side,
  dx = 0.781 mm; dt solved from η_char = η(10 1/s) = 9.264 Pa·s at
  ω = 1.0 → dt = 1.098e-5 s, u_lb = 0.0141, 455,325 steps. Local ω
  ≈ 1.0 near the lid, ~0.0003 in the capped plug (same documented
  exception as Couette-TC).
- Geometry / lid boundary follow OpenLB's own `laminar/cavity2d`
  example (IndicatorCuboid2D lid, `boundary::InterpolatedVelocity`
  on walls and lid, lid velocity via `momenta::setVelocity`).
- Per-cell shear rate and viscosity are dumped with every field
  snapshot (recovered from the non-equilibrium stress, the same
  formula `visco::OmegaFromCell` uses), so the video can show the
  yield pattern, not just the velocity.
- Validation: no analytic TC cavity solution exists, so it is
  twofold — (1) the TC dynamics is the code validated in Couette-TC
  (1.14% max torque error); (2) the cavity setup itself is validated
  by a Newtonian run of the same binary (`--newtonian --eta 1.0`,
  Re = 100) against Ghia, Ghia & Shin (1982) centreline velocities
  (table copied from the values embedded in OpenLB's own
  `laminar/cavity2d/cavity2d.cpp`).

## Two setup bugs found and fixed

1. **No lid (first Newtonian run: mean KE = 0 exactly).** The lid
   rename `geometry.rename(2, 3, 1, lid)` printed
   `[SuperGeometry2D] error rename only implemented for overlap>=2`
   and silently did nothing — the case had set `OVERLAP = 1`
   (copied from Couette, where the 3-arg rename works). Fix:
   `OVERLAP = 2`. The run then drives the flow. Lesson: a rename
   that prints an error and continues is easy to miss — the KE = 0
   check caught it.
2. **Mesh over [0, L] has N+1 lattice points per side** (both
   endpoints included: 129×129 at N = 128; fluid 127×127, walls 1
   cell thick). The C++ loops and the video grid now use
   `grid_n = N+1`, with the lid row at index N. (This also matches
   cavity2d, whose Ghia y-coordinates run 0..128.)
3. **Newtonian check at ω = 1.0 gives u_lb = 0.13 (Ma = 0.23)** —
   too compressible for a Ghia comparison. Fix: the Newtonian run
   uses `--omega 1.5` (dt = 3.391e-5 s, u_lb = 0.0434, Ma = 0.075),
   recorded in `case.yaml`. The TC run stays at ω = 1.0
   (u_lb = 0.0141).

## Newtonian validation run (Re = 100)

- 147,456 steps, dt = 3.391e-5 s, ω = 1.5, ~160 s wall
  (≈ 15.3 MLUPS, 129×129 grid).
- Centreline u_x(y) at x = L/2 vs Ghia et al. (1982) Re = 100:
  **relative L2 error 2.84%** (16 points, lid point excluded) —
  **PASSED** (tolerance 10%). Point values track Ghia throughout
  (e.g. y/L = 0.500: −0.207 vs −0.206; y/L = 0.977: 0.817 vs 0.841).
- Mean KE 80–90% vs 90–100%: drift 0.000% (fully steady).
  Final max speed 0.9524 m/s (lid 1.0).
- Density drift −0.98% over the 5 s run, roughly linear in time
  (−0.20%/s) through the InterpolatedVelocity walls. It does not
  disturb the steady velocity field (KE drift 0, Ghia 2.84%),
  but it is real and reported, not hidden.

## TC run

- 455,325 steps, dt = 1.098e-5 s, ω = 1.0, **1315.3 s wall
  (5.76 MLUPS)**. The per-step cost is ~2.6x the Newtonian cavity
  (15.3 MLUPS) — the scalar TC collision path, same ratio as
  Couette-TC vs Newtonian Couette. (Absolute MLUPS varies with
  host load on this shared VM; Couette-TC measured 11–22 MLUPS on
  its re-run.)
- Mean KE 80–90% vs 90–100%: 0.0101194 → 0.0101194 m²/s²,
  **rel. drift 0.00009%** — fully steady. Final max speed
  0.9264 m/s.
- Density drift −0.25% over the run (−0.05%/s, linear), 4x smaller
  than the Newtonian run's.
- Steady field (t = 5 s, 16,129 fluid cells): **plug
  (η > 1000 Pa·s) = 38.7%** of cells, **capped (η = η_max) =
  25.1%**. Max speed 0.9265 m/s.
- What the steady field shows (see `cavity-steady.png`): a thin
  yielded layer under the lid (η ~ 1–10 Pa·s), an unyielded plug
  at the vortex core that rotates as a solid body (the yellow
  high-η blob in the middle — expected for a viscoplastic cavity:
  the vortex core is below yield), a yielded band between that
  core plug and the stagnant bottom plug, and the bottom corners
  fully capped. The η = 1000 Pa·s contour (white) marks the
  visible yield surface. The horizontal banding inside the core
  plug is grid-scale variation of the recovered shear rate where η
  sits at the cap (γ̇ ~ 1e-3 1/s, ω ~ 0.0003) — shown as computed,
  not smoothed or interpreted further.

## Videos (verified with ffprobe, frames checked visually)

- `cavity-tc.mp4` — 3.1 s, 1500×500, H.264, 91 KB: 22 frames
  (t = 0 → 5 s), three panels per frame — speed (with quiver
  arrows), log₁₀ viscosity with the η = 1000 Pa·s yield surface
  contoured in white, vorticity — plus title, colourbars, axis
  labels and the TC parameters on screen. Frame at t = 1.75 s
  already shows the steady yield pattern; the pattern is stable
  to t = 5 s.
- `cavity-steady.png` — the same three panels for the last field,
  for the README / report.

## RTX 3090 runtime estimate (for the `matteo` run — ESTIMATE, not a measurement)

Basis: the TC cavity run measured here is 455,325 steps × 16,641
lattice cells = 7.58e9 lattice updates, done in 1315 s wall at
5.76 MLUPS on this VM's serial CPU. The GPU figure to compare
with: OpenLB's 2D BGK kernel measures 711 MLUPS unoptimized /
1077 ± 10 MLUPS optimized on an RTX A5000 (98.9% of roofline;
Kummerländer et al., via the OpenLB roofline study), and ~2,500
MLUPS on an A100 in the Moreau et al. barn-airflow study (Results
in Engineering). The RTX 3090 (sm_86, ~936 GB/s) sits above the
A5000 (~768 GB/s), so a pure-BGK 2D ballpark of ~700–1100 MLUPS
is grounded — but two derates apply to THIS case: (1) the TC
local-ω path is scalar/non-vectorized on CPU and its GPU
throughput is unmeasured — assume 300–700 MLUPS effective;
(2) the grid is only 16,641 cells, far too small to saturate a
GPU — kernel-launch overhead (~5–10 µs/step × 455k steps ≈ 2–5 s)
sets a floor regardless of MLUPS. So: **solver time ≈ 10–60 s**
(7.58e9 updates at 300–700 MLUPS ≈ 11–25 s, plus the launch
floor), i.e. roughly 20–130x faster than the 1315 s here — not
the ~200x the peak-MLUPS ratio alone would suggest. The same
extrapolation puts the Couette-TC sweep (1.68e9 updates) at ~5–30 s
of solver time on the 3090. One-time cost not included: the CUDA
build of the OpenLB core + case under WSL (unmeasured here;
allow 10–30 min) — and the serial CPU build needs no CUDA at all
(see the repo README's WSL section). Verify by timing the first
GPU run on `matteo`; if it lands outside 10–60 s, the launch
floor and the TC-path derate are the first suspects.

## Files

- `cavity.cpp`, `Makefile`, `case.yaml`, `cavity.py`,
  `cavity_run.py`, `cavity_post.py`, `make_video.py`,
  `lattice.json`, `cavity_summary.json`, `cavity-steady.png`,
  `cavity-tc.mp4`, this report.
- `runs/newtonian/` and `runs/tc/` (`energy.csv`,
  `centerline.csv`, `field_*.csv`, `meta.json`) are git-ignored
  and stay local, as in the Couette cases.

## Still unverified (out of scope)

- 3D cavity (2D only, by choice).
- Finer resolution than N = 128 (the yield surface is ~2–3 cells
  wide in places; untested).
- MPI / GPU: serial CPU only here — the GPU run is `matteo`'s
  job, with the estimate above.
- The slow linear density drift (−0.20%/s Newtonian, −0.05%/s TC)
  through the InterpolatedVelocity walls was not chased further;
  a longer run would drift further, so a production run should
  log density and, if needed, re-zero or use a mass-conserving
  wall treatment.
