#!/usr/bin/env python3
"""Post-process the RBC single-cell shear run.

Reads the membrane node positions from the .vtp files the shipped
example writes (runs/<mode>/tmp/vtkData/VTKParticles/Particles_t<iT>.vtp)
and computes the Taylor deformation parameter exactly the way
OpenLB's own writeParticleInertia does (GSL path — no GSL installed
here, so this mirrors it in numpy):

  * volume inertia tensor per face (computeInertiaTensorVolume):
        pv = faceCentre - centre, sp = pv . n_hat
        I_xx += area/5 * ((pv.pv)*sp - pv.x^2*sp),  I_xy += area/5 * (-pv.x*pv.y*sp), ...
  * diagonalise (numpy eigh), eigenvalues ascending e0 <= e1 <= e2
  * principal semi-axes (a largest, c smallest):
        a = sqrt(2.5*(e1+e2-e0)/V), b = sqrt(2.5*(e0+e2-e1)/V),
        c = sqrt(2.5*(e0+e1-e2)/V)
  * D = (a-c)/(a+c); plotted as D/G vs kappa*t, kappa = shearRate/G.

Validation gate: steady D/G (mean over the last 3 snapshots) within
15% of the analytic 6.25 (Krüger et al. 2011; OpenLB's shipped
figure). Exits non-zero on failure. Writes deformation.csv and
deformation.png in examples/rbc-shear/.

Usage: python3 rbc_post.py [--runs runs]
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).parent
G = 0.01                 # dimensionless shear rate (case.yaml)
KAPPA = 0.0133333333333  # shearRate/G per lattice step (case.yaml)
D_OVER_G_ANALYTIC = 6.25
TOL = 0.15


def read_vtp(path):
    """Return (nodes (n,3), faces (m,3)) from an ASCII .vtp file."""
    text = path.read_text()
    m = re.search(r'Name="Points"[^>]*>(.*?)</DataArray>', text, re.S)
    if not m:
        raise ValueError(f"no Points array in {path}")
    pts = np.array([float(x) for x in m.group(1).split()], dtype=float)
    nodes = pts.reshape(-1, 3)
    m = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', text, re.S)
    if not m:
        raise ValueError(f"no connectivity array in {path}")
    conn = np.array([int(x) for x in m.group(1).split()], dtype=int)
    faces = conn.reshape(-1, 3)
    return nodes, faces


def deformation(nodes, faces):
    """(D, a, b, c, volume) from node positions + face connectivity."""
    v0 = nodes[faces[:, 0]]
    v1 = nodes[faces[:, 1]]
    v2 = nodes[faces[:, 2]]
    centres = (v0 + v1 + v2) / 3.0
    cross = np.cross(v1 - v0, v2 - v0)
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    normals = cross / (2.0 * areas)[:, None]  # unit normals
    centre = nodes.mean(axis=0)
    pv = centres - centre
    sp = np.einsum("ij,ij->i", pv, normals)
    # Volume by the divergence theorem (outward normals).
    vol = float(np.sum(np.einsum("ij,ij->i", centres, normals) * areas) / 3.0)
    if vol < 0:
        vol = -vol
        sp = -sp
    w = areas / 5.0 * sp
    pv2 = np.einsum("ij,ij->i", pv, pv)
    Ixx = np.sum(w * (pv2 - pv[:, 0] ** 2))
    Iyy = np.sum(w * (pv2 - pv[:, 1] ** 2))
    Izz = np.sum(w * (pv2 - pv[:, 2] ** 2))
    Ixy = np.sum(w * (-pv[:, 0] * pv[:, 1]))
    Ixz = np.sum(w * (-pv[:, 0] * pv[:, 2]))
    Iyz = np.sum(w * (-pv[:, 1] * pv[:, 2]))
    I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    evals = np.linalg.eigvalsh(I)
    e0, e1, e2 = sorted(evals)
    a = math.sqrt(max(0.0, 2.5 * (e1 + e2 - e0) / vol))
    b = math.sqrt(max(0.0, 2.5 * (e0 + e2 - e1) / vol))
    c = math.sqrt(max(0.0, 2.5 * (e0 + e1 - e2) / vol))
    D = (a - c) / (a + c) if (a + c) > 0 else float("nan")
    return D, a, b, c, vol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--validation-tol", type=float, default=TOL)
    args = ap.parse_args()
    vtp_dir = Path(args.runs) / "single" / "tmp" / "vtkData" / "VTKParticles"
    files = sorted(vtp_dir.glob("Particles_t*.vtp"),
                   key=lambda p: int(p.stem.split("_t")[1]))
    if not files:
        print(f"no .vtp files in {vtp_dir} — run the case first")
        return 1
    rows = []
    for f in files:
        iT = int(f.stem.split("_t")[1])
        nodes, faces = read_vtp(f)
        D, a, b, c, vol = deformation(nodes, faces)
        kt = iT * KAPPA
        rows.append({"iT": iT, "kappa_t": kt, "D": D, "D_over_G": D / G,
                     "a": a, "b": b, "c": c, "volume": vol})
        print(f"kappa*t={kt:6.1f}  D={D:.4f}  D/G={D/G:6.3f}  "
              f"a={a:.3f} b={b:.3f} c={c:.3f} V={vol:.1f}")
    steady = np.mean([r["D_over_G"] for r in rows[-3:]])
    err = abs(steady - D_OVER_G_ANALYTIC) / D_OVER_G_ANALYTIC
    passed = err <= args.validation_tol
    print(f"steady D/G (last 3 snapshots) = {steady:.3f} vs analytic "
          f"{D_OVER_G_ANALYTIC} (err {err*100:.1f}%, tol "
          f"{args.validation_tol*100:.0f}%) — {'PASSED' if passed else 'FAILED'}")
    print(f"shipped dx=1 run reaches ~6.8; this run: {steady:.3f}")

    with open(HERE / "deformation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["iT", "kappa_t", "D", "D_over_G",
                                           "a", "b", "c", "volume"])
        w.writeheader()
        w.writerows(rows)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([r["kappa_t"] for r in rows], [r["D_over_G"] for r in rows],
            "o-", label="OpenLB (this run, dx=1)")
    ax.axhline(D_OVER_G_ANALYTIC, color="k", ls="--",
               label=f"analytic steady D/G = {D_OVER_G_ANALYTIC}")
    ax.axhline(6.8, color="gray", ls=":", label="shipped dx=1 run (~6.8)")
    ax.set_xlabel("normalized time κt")
    ax.set_ylabel("scaled deformation D/G")
    ax.set_title("Single deformable capsule in shear flow (Krüger et al. 2011)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "deformation.png", dpi=150)
    print("wrote deformation.csv and deformation.png")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
