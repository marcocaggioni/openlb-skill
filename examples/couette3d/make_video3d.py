#!/usr/bin/env python3
"""Render a 3D video from the OpenLB 3D Couette run (runs/gamma_10).

Output (in examples/couette3d/):
  couette3d-gamma10.mp4 — spin-up at gamma = 10 1/s, rendered as the
  actual 3D annulus: the grey bob in the middle, the sheared fluid as
  concentric cylindrical shells coloured by speed, a 90-degree wedge
  cut away so the radial shear gradient is visible on the two cut
  faces and on the top cap. The camera slowly orbits once over the
  z-invariant (laminar, periodic in Z), so the mid-plane slice
  (field_mid_*.csv) fully determines the 3D field: speed depends on
  radius only, and the render extrudes that profile along Z.

Frames are rendered with matplotlib 3D (Agg) and encoded with ffmpeg
(libx264). Requires the gamma = 10 run to have finished first.

Colour scale is fixed (vmax = 0.04 m/s, just above the wall speed
0.036 m/s) so the spin-up is visible: early frames are dark, the
steady frame is bright at the bob.
"""
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
RI, RO, L = 0.020, 0.025, 0.060
VMAX = 0.04
WEDGE = 270.0  # degrees of the annulus shown; 90-degree wedge cut away
T_ANALYTIC = 0.00150796


def read_mid_plane(path):
    """Mid-plane slice -> (radii_m, speed_m_s) binned by radius, 1 cell wide.

    Cell centres are at integer lattice coordinates: OpenLB's
    Cuboid::getPhysR(latticeR) = origin + latticeR*dx
    (src/geometry/cuboid.h:117), so a cell (iX, iY) sits at radius
    hypot(iX - 50, iY - 50) lattice units from the annulus centre
    (50, 50). Cells with radius exactly 40 or 50 are wall cells
    (material 3 / 2, zero velocity) and are excluded — keeping them
    would drag the innermost/outermost bins down. Only fluid cells
    (40 < r < 50) enter the bins.
    """
    lat = json.loads((HERE / "lattice3d.json").read_text())
    dx = lat["dx_m"]
    c = lat["nx"] // 2
    # Integer lattice radii of the walls (RI/dx = 40, RO/dx = 50 exactly;
    # do NOT compute them as floats — 0.02/0.0005 rounds to 39.999...,
    # which would let the r = 40 wall cells, velocity exactly zero,
    # leak into the innermost bin).
    r_lo, r_hi = int(round(RI / dx)), int(round(RO / dx))
    r_bins = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            r_lat = math.hypot(int(row["iX"]) - c, int(row["iY"]) - c)
            if r_lat <= r_lo or r_lat >= r_hi:
                continue
            sp = math.hypot(float(row["ux_m_s"]), float(row["uy_m_s"]))
            b = int(r_lat)
            r_bins.setdefault(b, []).append(sp)
    radii = np.array([(b + 0.5) * dx for b in sorted(r_bins)])
    speeds = np.array([float(np.mean(r_bins[b])) for b in sorted(r_bins)])
    return radii, speeds


def read_torque():
    ts, tq = [], []
    with open(HERE / "runs" / "gamma_10" / "torque.csv") as fh:
        for row in csv.DictReader(fh):
            ts.append(float(row["time_s"]))
            tq.append(float(row["torque_N_m"]))
    return np.array(ts), np.array(tq)


def surf_cylinder(ax, radius, z0, z1, theta_max_deg, color, alpha,
                  n_theta=72, n_z=13):
    """Cylindrical shell of fixed radius; color is a single RGBA or
    a (n_theta, n_z) array of RGBA face colours."""
    th = np.linspace(0, np.deg2rad(theta_max_deg), n_theta)
    z = np.linspace(z0, z1, n_z)
    TH, Z = np.meshgrid(th, z)
    X = radius * np.cos(TH)
    Y = radius * np.sin(TH)
    if isinstance(color, np.ndarray) and color.ndim == 3:
        fc = color
    else:
        fc = np.empty(TH.shape + (4,), dtype=float)
        fc[:] = color
    ax.plot_surface(X, Y, Z, facecolors=fc, shade=False, alpha=alpha,
                    edgecolor="none", antialiased=False)


def surf_radial_face(ax, theta_deg, radii, speeds, norm, cmap):
    """Flat cut face at fixed azimuth: rectangle r in [Ri, Ro], z in [0, L],
    coloured by the speed profile."""
    th = np.deg2rad(theta_deg)
    z = np.linspace(0, L, 13)
    R, Z = np.meshgrid(radii, z)
    X = R * np.cos(th)
    Y = R * np.sin(th)
    fc = cmap(norm(speeds))[np.newaxis, :, :].repeat(len(z), axis=0)
    # meshgrid gives (z, r); facecolors must match (n_z, n_r, 4)
    ax.plot_surface(X, Y, Z, facecolors=fc, shade=False, alpha=1.0,
                    edgecolor="none", antialiased=False)


def surf_top_cap(ax, radii, speeds, norm, cmap, theta_max_deg):
    """Annulus disc at z = L, coloured by the speed profile."""
    th = np.linspace(0, np.deg2rad(theta_max_deg), 72)
    R, TH = np.meshgrid(radii, th)
    X = R * np.cos(TH)
    Y = R * np.sin(TH)
    Z = np.full_like(X, L)
    fc = cmap(norm(speeds))[np.newaxis, :, :].repeat(len(th), axis=0)
    ax.plot_surface(X, Y, Z, facecolors=fc, shade=False, alpha=1.0,
                    edgecolor="none", antialiased=False)


def render_frame(radii, speeds, t, torque, azim, out_png):
    norm = Normalize(vmin=0, vmax=VMAX)
    cmap = cm.viridis
    fig = plt.figure(figsize=(7.2, 7.2))
    ax = fig.add_subplot(111, projection="3d")
    # --- solid bob (inner cylinder), full 360 degrees, light grey ---
    surf_cylinder(ax, RI, 0, L, 360, (0.62, 0.62, 0.62, 1.0), 1.0)
    # --- fluid shells at the bin radii, wedge cut away ---
    for r, v in zip(radii, speeds):
        if r <= RI or r >= RO:
            continue
        surf_cylinder(ax, r, 0, L, WEDGE, cmap(norm(v)), 1.0,
                      n_theta=64, n_z=9)
    # --- cut faces at both edges of the wedge ---
    surf_radial_face(ax, 0.0, radii, speeds, norm, cmap)
    surf_radial_face(ax, WEDGE, radii, speeds, norm, cmap)
    # --- top cap ---
    surf_top_cap(ax, radii, speeds, norm, cmap, WEDGE)
    ax.set_xlim(-0.027, 0.027)
    ax.set_ylim(-0.027, 0.027)
    ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, L / (2 * RO)))
    ax.set_axis_off()
    ax.view_init(elev=18, azim=azim)
    ax.set_title(f"OpenLB 3D Couette, $\\dot\\gamma$ = 10 1/s, t = {t:.2f} s\n"
                 f"torque = {torque:.5f} N·m (analytic {T_ANALYTIC:.5f})",
                 fontsize=11)
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(mappable, ax=ax, label="speed (m/s)", shrink=0.65, pad=0.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def main():
    lat = json.loads((HERE / "lattice3d.json").read_text())
    dt = lat["dt_s"]
    t_tq, v_tq = read_torque()
    fields = sorted((HERE / "runs" / "gamma_10").glob("field_mid_*.csv"),
                    key=lambda p: int(p.stem.split("_")[2]))
    print(f"{len(fields)} snapshots, torque samples: {len(t_tq)}")
    tmp = Path("/tmp/couette3d-frames")
    tmp.mkdir(parents=True, exist_ok=True)
    # One full camera orbit over the video; each snapshot held 2 frames.
    n_frames = 2 * len(fields)
    for k, f in enumerate(fields):
        idx = int(f.stem.split("_")[2])
        t = idx * 1000 * dt
        radii, speeds = read_mid_plane(f)
        torque = float(np.interp(t, t_tq, v_tq)) if len(t_tq) else 0.0
        torque += 0.0  # -0.0 (first sample, exactly zero torque) -> 0.0
        for j in range(2):
            frame = 2 * k + j
            azim = -60 + 360.0 * frame / max(n_frames - 1, 1)
            render_frame(radii, speeds, t, torque, azim,
                         tmp / f"frame_{frame:03d}.png")
        print(f"snapshot {k+1}/{len(fields)}  t={t:.3f} s  "
              f"max speed={speeds.max():.4f} m/s  torque={torque:.5f} N·m")
    out = HERE / "couette3d-gamma10.mp4"
    subprocess.run(["ffmpeg", "-y", "-v", "error", "-framerate", "10",
                    "-i", str(tmp / "frame_%03d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", str(out)],
                   check=True)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
