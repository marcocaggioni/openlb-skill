# Couette rheometer, 3D — OpenLB example

Same geometry and fluid as the 2D case (`../couette/case.yaml`):
R_i = 0.020 m rotating bob, R_o = 0.025 m fixed cup, L = 0.060 m,
Newtonian η = 1.0 Pa·s, ρ = 1000 kg/m³ — but the full 3D annulus
(D3Q19) instead of a 2D cross-section, periodic along the cylinder
axis Z, so there are no end effects and the infinite-cylinder
assumption of the 2D case carries over unchanged.

This example runs γ̇ = 10 1/s only (Ω = 1.8 rad/s, analytic torque
0.00150796 N·m). See `RUN_REPORT.md` for the result.

## Files

- `case.yaml` — the case contract (same numbers as the 2D case).
- `lattice3d.py` — physical → lattice parameters; writes `lattice3d.json`.
  dx = 0.5 mm (10 cells across the gap), dt = 4.167e-5 s, 24,000 steps,
  grid 100×100×120 = 1.2M cuboid cells.
- `couette3d.cpp` — the case. D3Q19/BGK, Bouzidi interpolated
  bounce-back moving inner wall (same pattern as the 2D case),
  bounce-back outer wall, torque by momentum exchange at the
  inner-wall links (3D conversion ρ·dx⁵/dt²), mid-plane velocity
  slices every 1000 steps for the video.
- `make_video3d.py` — renders `couette3d-gamma10.mp4`: the 3D annulus
  with a wedge cut away, shells coloured by speed, camera orbiting
  once over the spin-up (matplotlib 3D + ffmpeg).
- `RUN_REPORT.md` — run report.

## Build and run (serial CPU)

```bash
cd examples/couette3d
uv run python lattice3d.py   # writes lattice3d.json, checks u_lb < 0.1
make                          # GNU Make, OpenLB at ~/openlb
./couette3d --gamma 10 --out runs/gamma_10
uv run python make_video3d.py
```

`runs/` is git-ignored and stays local.
