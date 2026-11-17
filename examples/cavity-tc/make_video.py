#!/usr/bin/env python3
"""Render the cavity video from the field snapshots (runs/tc/field_*.csv).

Outputs (in examples/cavity-tc/):
  cavity-tc.mp4        — 3-panel video: speed | log10(viscosity) |
                         vorticity, with the yield surface (eta = 1000
                         Pa.s) contoured on the viscosity panel.
  cavity-steady.png    — the same 3 panels for the last (steady) field,
                         for the README / report.

Field CSVs carry (iX, iY, ux, uy, gamma, eta) for fluid cells only on
a 128x128 grid; the lid row (top) is set to U_lid and all other
non-fluid cells to zero when the full grid is rebuilt. Vorticity is
computed in Python from the velocity grid (numpy gradient).

The yield indicator: log10(eta) spans ~0 (eta ~ 1 Pa.s, yielded near
the lid) to ~4.85 (eta_max = 7.02e4 Pa.s, the capped plug). The white
contour marks eta = 1000 Pa.s (gamma ~ 0.07 1/s) — the visible yield
surface between the flowing layer and the plug.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
RUN = HERE / "runs" / "tc"
ETA_CONTOUR = 1000.0  # Pa.s — the plotted yield surface


def read_field(path, meta):
    N = meta["grid_n"]  # lattice points per side: N intervals + endpoints
    dx = meta["dx_m"]
    UX = np.zeros((N, N))
    UY = np.zeros((N, N))
    ETA = np.full((N, N), np.nan)
    GAMMA = np.full((N, N), np.nan)
    fluid = np.zeros((N, N), dtype=bool)
    with open(path) as fh:
        for row in csv.DictReader(fh):
            iX, iY = int(row["iX"]), int(row["iY"])
            UX[iY, iX] = float(row["ux_m_s"])
            UY[iY, iX] = float(row["uy_m_s"])
            ETA[iY, iX] = float(row["eta_Pa_s"])
            GAMMA[iY, iX] = float(row["gamma_1_s"])
            fluid[iY, iX] = True
    # Lid row (top, index N-1 of the grid) moves at U_lid; other
    # non-fluid cells (walls) are zero.
    UX[N - 1, :] = meta["U_lid_m_s"]
    ETA[N - 1, :] = np.nanmin(ETA)
    return UX, UY, ETA, GAMMA, fluid


def render_frame(path, out_png, meta, t):
    N = meta["grid_n"]
    dx = meta["dx_m"]
    L = meta["L_m"]
    U = meta["U_lid_m_s"]
    UX, UY, ETA, GAMMA, fluid = read_field(path, meta)
    speed = np.hypot(UX, UY)
    # Vorticity (1/s) from the velocity grid.
    dudy, dudx = np.gradient(UY, dx, dx)
    dvdx, dvdy = np.gradient(UX, dx, dx)
    vort = dvdx - dudy
    x = (np.arange(N) + 0.5) * dx
    y = (np.arange(N) + 0.5) * dx
    X, Y = np.meshgrid(x, y)
    logeta = np.log10(np.where(np.isfinite(ETA), ETA, np.nan))

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"Lid-driven cavity — TC fluid: L = {L} m, U_lid = {U} m/s, "
        f"$\\sigma_y$ = {meta['sigma_y_Pa']:.0f} Pa, "
        f"$\\dot\\gamma_c$ = {meta['gamma_dot_c_1_s']:.0f} 1/s, "
        f"$\\eta_{{bg}}$ = {meta['eta_bg_Pa_s']} Pa·s   |   t = {t:.2f} s",
        fontsize=11)
    im0 = ax[0].imshow(speed, extent=[0, L, 0, L], origin="lower",
                       cmap="viridis", vmin=0, vmax=U)
    ax[0].quiver(X[::8, ::8], Y[::8, ::8], UX[::8, ::8], UY[::8, ::8],
                 color="white", scale=25, width=0.004)
    ax[0].set_title("speed (m/s)")
    fig.colorbar(im0, ax=ax[0], label="speed (m/s)", shrink=0.85)
    im1 = ax[1].imshow(logeta, extent=[0, L, 0, L], origin="lower",
                       cmap="plasma", vmin=0, vmax=5)
    ax[1].contour(X, Y, np.where(np.isfinite(ETA), ETA, np.nan),
                  levels=[ETA_CONTOUR], colors="white", linewidths=1.5)
    ax[1].set_title(f"log10 viscosity (Pa·s) — white: η = {ETA_CONTOUR:.0f} (yield surface)")
    fig.colorbar(im1, ax=ax[1], label="log10 η (Pa·s)", shrink=0.85)
    vmax_v = 60.0
    im2 = ax[2].imshow(vort, extent=[0, L, 0, L], origin="lower",
                       cmap="RdBu_r", vmin=-vmax_v, vmax=vmax_v)
    ax[2].set_title("vorticity (1/s)")
    fig.colorbar(im2, ax=ax[2], label="vorticity (1/s)", shrink=0.85)
    for a in ax:
        a.set_xlabel("x (m)")
        a.set_ylabel("y (m)")
        a.set_aspect("equal")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def main():
    meta = json.loads((RUN / "meta.json").read_text())
    fields = sorted(RUN.glob("field_*.csv"),
                    key=lambda p: int(p.stem.split("_")[1]))
    snap_dt = meta["snap_interval_s"]
    run_time = meta["run_time_s"]
    tmp = Path("/tmp/cavity-tc-frames")
    tmp.mkdir(parents=True, exist_ok=True)
    for k, f in enumerate(fields):
        # The last field is written after the loop at t = run_time, not
        # at (n_fields-1)*snap_dt.
        t = min(k * snap_dt, run_time)
        render_frame(f, tmp / f"frame_{k:03d}.png", meta, t)
        print(f"frame {k}: t={t:.2f} s", flush=True)
    # Steady figure (last frame) for the README / report.
    render_frame(fields[-1], HERE / "cavity-steady.png", meta, run_time)
    subprocess.run(["ffmpeg", "-y", "-v", "error", "-framerate", "7",
                    "-i", str(tmp / "frame_%03d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p",
                    str(HERE / "cavity-tc.mp4")], check=True)
    print("wrote", HERE / "cavity-tc.mp4", "and cavity-steady.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
