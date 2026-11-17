#!/usr/bin/env python3
"""Post-process the cavity runs.

Newtonian validation run (runs/newtonian/): compare the centreline
u_x(y) at x = L/2 with Ghia, Ghia & Shin (1982) Re = 100 (table in
cavity.py) and report the relative L2 error. Gate: < 10% (same spirit
as the Couette gates). Exits non-zero on failure.

TC run (runs/tc/): check steady state from energy.csv (mean-KE drift
between the 80–90% and 90–100% windows), and report the steady
field: plug fraction (cells with eta > 1000 Pa.s — the visible yield
surface; eta = 1000 corresponds to gamma ~ 0.07 1/s), capped fraction
(eta == eta_max), max speed, and mean KE. Writes cavity_summary.json
and cavity-steady.png (the 3-panel steady figure, same as the video's
last frame).

Usage: python3 cavity_post.py [--runs runs]
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from cavity import GHIA_RE100_U, GHIA_RE100_Y, load_case  # noqa: E402

HERE = Path(__file__).parent
ETA_PLUG = 1000.0  # Pa.s; the yield-surface threshold used in the video


def read_centerline(path):
    ys, us = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            ys.append(float(row["y_m"]))
            us.append(float(row["ux_m_s"]))
    return np.array(ys), np.array(us)


def read_energy(path):
    t, ke, umax, rho = [], [], [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t.append(float(row["time_s"]))
            ke.append(float(row["mean_ke_m2_s2"]))
            umax.append(float(row["max_speed_m_s"]))
            rho.append(float(row.get("mean_rho", 1.0)))
    return np.array(t), np.array(ke), np.array(umax), np.array(rho)


def ghia_check(run_dir, L, U):
    ys, us = read_centerline(run_dir / "centerline.csv")
    y_frac = np.array(GHIA_RE100_Y) / 128.0
    u_ghia = np.array(GHIA_RE100_U) * U
    # Simulated u_x interpolated onto Ghia's y positions (drop the
    # lid point y/L = 1, where u is imposed = U by construction).
    u_sim = np.interp(y_frac, ys / L, us)
    m = y_frac < 1.0
    err = np.sqrt(np.mean((u_sim[m] - u_ghia[m]) ** 2)) / np.sqrt(
        np.mean(u_ghia[m] ** 2))
    print(f"Ghia Re=100 centreline: relative L2 error = {err*100:.2f}% "
          f"({m.sum()} points, lid point excluded)")
    for yf, ug, usi in zip(y_frac[m], u_ghia[m], u_sim[m]):
        print(f"  y/L={yf:.3f}  Ghia u/U={ug/U:7.3f}  sim u/U={usi/U:7.3f}")
    return err


def energy_drift(run_dir, U):
    t, ke, umax, rho = read_energy(run_dir / "energy.csv")
    n = len(ke)
    k1 = ke[int(0.8 * n):int(0.9 * n)].mean()
    k2 = ke[int(0.9 * n):].mean()
    drift = abs(k2 - k1) / k1 if k1 else float("nan")
    rho_drift = (rho[-1] - rho[0]) / rho[0] if rho[0] else float("nan")
    print(f"energy: mean KE 80-90% = {k1:.6e}, 90-100% = {k2:.6e} "
          f"m2/s2, rel drift = {drift*100:.3f}%, "
          f"final max speed = {umax[-1]:.4f} m/s (lid {U:.1f}), "
          f"density drift = {rho_drift*100:+.4f}%")
    return drift, k2, umax[-1], rho_drift


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    args = ap.parse_args()
    runs = Path(args.runs)
    case = load_case()
    L = case["geometry"]["L"]
    U = case["flow"]["U_lid_m_s"]
    summary = {}

    newt = runs / "newtonian"
    if (newt / "centerline.csv").exists():
        print("== Newtonian validation run (Re = 100) ==")
        drift, _, _, rho_drift = energy_drift(newt, U)
        err = ghia_check(newt, L, U)
        passed = err <= 0.10
        print(f"Newtonian cavity validation: {'PASSED' if passed else 'FAILED'} "
              f"(tol 10%)")
        summary["newtonian"] = {"ghia_rel_L2": float(err),
                                "energy_drift": float(drift),
                                "density_drift": float(rho_drift),
                                "passed": bool(passed)}
    else:
        print("Newtonian run not found — run cavity_run.py --newtonian first")
        summary["newtonian"] = None

    tc = runs / "tc"
    if (tc / "meta.json").exists():
        print("== TC run ==")
        drift, ke_final, umax_final, rho_drift = energy_drift(tc, U)
        meta = json.loads((tc / "meta.json").read_text())
        last = sorted(tc.glob("field_*.csv"),
                      key=lambda p: int(p.stem.split("_")[1]))[-1]
        etas, gammas, speeds = [], [], []
        with open(last) as fh:
            for row in csv.DictReader(fh):
                etas.append(float(row["eta_Pa_s"]))
                gammas.append(float(row["gamma_1_s"]))
                speeds.append(math.hypot(float(row["ux_m_s"]),
                                         float(row["uy_m_s"])))
        etas = np.array(etas)
        plug_frac = float((etas > ETA_PLUG).mean())
        capped_frac = float((etas >= meta["eta_max_Pa_s"] * 0.999).mean())
        print(f"steady field {last.name}: {len(etas)} fluid cells, "
              f"plug (eta>{ETA_PLUG:.0f} Pa.s) = {plug_frac*100:.1f}%, "
              f"capped (eta=eta_max) = {capped_frac*100:.1f}%, "
              f"max speed = {max(speeds):.4f} m/s")
        summary["tc"] = {"energy_drift": float(drift),
                         "density_drift": float(rho_drift),
                         "mean_KE_final": float(ke_final),
                         "max_speed_final": float(umax_final),
                         "plug_fraction": plug_frac,
                         "capped_fraction": capped_frac,
                         "eta_plug_threshold_Pa_s": ETA_PLUG,
                         "last_field": last.name}
    else:
        print("TC run not found — run cavity_run.py first")
        summary["tc"] = None

    (HERE / "cavity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote cavity_summary.json")
    # Fail only on a run that is present and failed its check.
    failed = False
    if summary["newtonian"] is not None and not summary["newtonian"]["passed"]:
        failed = True
    if summary["tc"] is not None and summary["tc"]["energy_drift"] > 0.05:
        print("WARNING: TC run energy drift > 5% — not converged")
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
