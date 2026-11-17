#!/usr/bin/env python3
"""Couette-TC post-processing: torque.csv -> flow curve + TC validation.

Input per sweep point: runs/gamma_<value>/torque.csv with columns
    time_s, torque_N_m_per_m
(torque per unit length; the 2D case. Multiply by L for 3D torque.)

Conversions (same as the Newtonian Couette example):
    gamma_i = 2*Omega*Ro^2 / (Ro^2 - Ri^2)      Omega from gamma via case.yaml
    tau_i   = T / (2*pi*Ri^2*L),  T = torque_per_m * L
    eta_app = tau_i / gamma_i   (label rate — NOT the TC inner-wall rate)

Validation gate: the simulated torque must match the ANALYTIC TC
torque at each Omega (tc-analytic.py math, copied into lattice.py)
within 10% (default). Exits non-zero on failure. The comparison is
torque vs analytic torque — never eta_app vs the TC law at the label
rate, because the TC inner-wall rate exceeds the label rate.

Also, per point: a profile cross-check from runs/gamma_*/profile.csv
(r_m, vtheta_m_s, vtheta_analytic_m_s) — median relative error of the
simulated azimuthal velocity vs the analytic TC profile over the
yielded bins, and a profile-derived torque from the two innermost
fluid bins (shear rate x TC viscosity there), for comparison with the
momentum-exchange torque (primary).

Usage: python3 couette_post.py [--runs runs] [--validation-tol 0.10]
Writes flow_curve.csv and flow_curve.png in examples/couette-tc/.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from lattice import load_case, tc_eta, tc_gamma_of_tau, tc_torque_of_omega  # noqa: E402

HERE = Path(__file__).parent


def read_torque(path):
    times, torques = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            times.append(float(row["time_s"]))
            torques.append(float(row["torque_N_m_per_m"]))
    return times, torques


def steady_mean(times, torques, frac):
    n = len(times)
    k = max(1, int(n * frac))
    window = torques[-k:]
    mean = sum(window) / len(window)
    sd = (sum((x - mean) ** 2 for x in window) / len(window)) ** 0.5
    # drift check: first vs second half of the steady window
    h = max(1, len(window) // 2)
    m1 = sum(window[:h]) / h
    m2 = sum(window[h:]) / (len(window) - h)
    drift = abs(m2 - m1) / mean if mean else float("nan")
    return mean, sd, drift


def read_density(path):
    """Mean density drift over the run: (last-first)/first."""
    rhos = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rhos.append(float(row["mean_rho"]))
    if not rhos:
        return float("nan")
    return (rhos[-1] - rhos[0]) / rhos[0] if rhos[0] else float("nan")


def read_profile(path):
    rs, vs, va = [], [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rs.append(float(row["r_m"]))
            vs.append(float(row["vtheta_m_s"]))
            va.append(float(row["vtheta_analytic_m_s"]))
    return rs, vs, va


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--validation-tol", type=float, default=0.10)
    args = ap.parse_args()

    case = load_case()
    g, f, s = case["geometry"], case["fluid"], case["sweep"]
    Ri, Ro, L = g["R_i"], g["R_o"], g["L"]
    sy, gc, ebg = f["sigma_y_Pa"], f["gamma_dot_c_1_s"], f["eta_bg_Pa_s"]
    eta_max = f["eta_max_Pa_s"]
    frac = s["steady_window_fraction"]
    dx = (Ro - Ri) / case["lattice"]["cells_across_gap"]
    rows = []
    for gamma in s["shear_rate_1_s"]:
        path = Path(args.runs) / f"gamma_{gamma}" / "torque.csv"
        if not path.exists():
            print(f"missing {path} — run the sweep first")
            return 1
        times, torques = read_torque(path)
        T_per_m, sd, drift = steady_mean(times, torques, frac)
        dens_path = Path(args.runs) / f"gamma_{gamma}" / "density.csv"
        rho_drift = read_density(dens_path) if dens_path.exists() else float("nan")
        T = T_per_m * L
        Omega = gamma * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro)
        tau = T / (2.0 * math.pi * Ri * Ri * L)
        T_an = tc_torque_of_omega(Omega, Ri, Ro, L, sy, gc, ebg)
        err = abs(T - T_an) / T_an if T_an else float("nan")
        r_y = math.sqrt(T_an / (2.0 * math.pi * L * sy))
        # profile cross-check
        prof_path = Path(args.runs) / f"gamma_{gamma}" / "profile.csv"
        prof_err, prof_torque = float("nan"), float("nan")
        if prof_path.exists():
            rs, vs, va = read_profile(prof_path)
            # Median profile error NORMALIZED BY THE MAX ANALYTIC SPEED:
            # pointwise |v-a|/a inflates near the yield surface, where
            # the analytic velocity -> 0, so it is not the honest metric.
            a_max = max(va) if va else 0.0
            rels = sorted(abs(v - a) / a_max for r, v, a in zip(rs, vs, va)
                          if r < min(r_y, Ro) and a_max > 0)
            prof_err = rels[len(rels) // 2] if rels else float("nan")
            # profile torque from the two innermost fluid bins
            # (nearest-radius lookup: float keys from C++ vs Python
            # differ in the last bit, so no exact dict match)
            def v_at(r_target):
                j = min(range(len(rs)), key=lambda i: abs(rs[i] - r_target))
                return vs[j] if abs(rs[j] - r_target) < dx / 2 else None
            r0, r1 = (int(Ri / dx) + 1.5) * dx, (int(Ri / dx) + 2.5) * dx
            v0 = v_at(r0)
            v1 = v_at(r1)
            if v0 is not None and v1 is not None:
                gamma_prof = -r0 * ((v1 / r1) - (v0 / r0)) / (r1 - r0)
                tau_prof = tc_eta(gamma_prof, sy, gc, ebg, eta_max) * gamma_prof
                prof_torque = tau_prof * 2 * math.pi * Ri * Ri * L
        rows.append({"gamma": gamma, "omega": Omega, "tau": tau,
                     "eta_app": tau / gamma, "T": T, "T_analytic": T_an,
                     "rel_err": err, "torque_sd": sd * L, "drift": drift,
                     "R_y_mm": min(r_y, Ro) * 1000,
                     "yielded_throughout": r_y >= Ro,
                     "prof_err": prof_err, "prof_torque": prof_torque,
                     "rho_drift": rho_drift})
        print(f"gamma={gamma:>6} 1/s  tau={tau:9.3f} Pa  eta_app={tau/gamma:8.3f} "
              f"Pa.s  T={T:.6e} N.m  analytic err={err*100:6.2f}%  "
              f"R_y={min(r_y, Ro)*1000:5.2f} mm  drift={drift*100:.2f}%  "
              f"prof_med_err={prof_err*100:5.2f}%  "
              f"T_prof={prof_torque:.6e} N.m  "
              f"rho_drift={rho_drift*100:+.4f}%")

    max_err = max(r["rel_err"] for r in rows)
    passed = max_err <= args.validation_tol
    print(f"TC validation: max analytic torque error {max_err*100:.2f}% "
          f"(tol {args.validation_tol*100:.0f}%) — {'PASSED' if passed else 'FAILED'}")

    with open(HERE / "flow_curve.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gamma_1_s", "omega_rad_s", "tau_Pa",
                                           "eta_app_Pa_s", "torque_N_m",
                                           "torque_analytic_N_m", "rel_err",
                                           "R_yield_mm", "yielded_throughout",
                                           "torque_sd_N_m", "steady_drift",
                                           "profile_med_rel_err",
                                           "torque_from_profile_N_m",
                                           "density_drift"])
        w.writeheader()
        for r in rows:
            w.writerow({"gamma_1_s": r["gamma"], "omega_rad_s": r["omega"],
                        "tau_Pa": r["tau"], "eta_app_Pa_s": r["eta_app"],
                        "torque_N_m": r["T"],
                        "torque_analytic_N_m": r["T_analytic"],
                        "rel_err": r["rel_err"],
                        "R_yield_mm": r["R_y_mm"],
                        "yielded_throughout": r["yielded_throughout"],
                        "torque_sd_N_m": r["torque_sd"],
                        "steady_drift": r["drift"],
                        "profile_med_rel_err": r["prof_err"],
                        "torque_from_profile_N_m": r["prof_torque"],
                        "density_drift": r["rho_drift"]})

    gammas = [r["gamma"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(gammas, [r["tau"] for r in rows], "o-", label="OpenLB (TC)")
    ax[0].loglog(gammas,
                 [r["T_analytic"] / (2 * math.pi * Ri * Ri * L) for r in rows],
                 "k--", label="analytic TC")
    ax[0].set_xlabel("shear rate label (1/s)")
    ax[0].set_ylabel("shear stress at inner wall (Pa)")
    ax[0].set_title("Couette flow curve — TC fluid (OpenLB)")
    ax[0].legend()
    ax[1].loglog(gammas, [r["eta_app"] for r in rows], "o-")
    ax[1].set_xlabel("shear rate label (1/s)")
    ax[1].set_ylabel("apparent viscosity (Pa.s)")
    ax[1].set_title("Apparent viscosity (label rate)")
    fig.tight_layout()
    fig.savefig(HERE / "flow_curve.png", dpi=150)
    print("wrote flow_curve.csv and flow_curve.png")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
