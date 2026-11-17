#!/usr/bin/env python3
"""Physical -> lattice parameters for the OpenLB Couette-TC case.

Reads examples/couette-tc/case.yaml. The TC viscosity spans orders of
magnitude across the gap (eta ~ 1 Pa*s at the inner wall at high rate,
~7e4 Pa*s capped in the plug), so dt is solved PER SWEEP POINT from
the viscosity at the label rate, eta_label = eta_TC(gamma_label), at
omega = 1.0 — exactly what the C++ case does. Writes lattice.json
next to case.yaml. Run this first; the C++ case must agree with it.

Formulas (D2Q9, cs^2 = 1/3):
    dx     = gap / N_gap
    nu_lb  = nu * dt / dx^2
    omega  = 1 / (3*nu_lb + 0.5)      ->  nu_lb = (1/omega - 0.5)/3
    dt     = nu_lb * dx^2 / nu_label  (nu_label = eta_label / rho)
    u_lb   = u_phys * dt / dx         ->  must stay < max_lattice_velocity
    run_time = diffusion_times * gap^2 / nu_label
    torque conversion: lattice force -> N uses rho * dx^4 / dt^2

With run_time = 40 gap^2/nu_label and omega = 1.0, the step count is
40*6*(gap/dx)^2 = 24,000 at every point.

The local omega at the inner wall is estimated from the analytic TC
inner-wall shear rate (tc-analytic.py math, copied here so the repo is
self-contained); the plug omega from eta_max. The yielded-region omega
must stay in ~[0.8, 1.9]; the plug omega is ~0.02 by construction of
the eta cap — see case.yaml, it is reported, not hidden.
"""
import json
import math
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # keep the file usable without PyYAML for a quick check
    yaml = None

CASE = Path(__file__).with_name("case.yaml")


def load_case():
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    return yaml.safe_load(CASE.read_text())


# --- TC analytic reference (copy of ~/workspace/tc-analytic.py) ---
def tc_eta(gamma, sigma_y, gamma_c, eta_bg, eta_max):
    if gamma <= 0:
        return eta_max
    return min(sigma_y / gamma + sigma_y / math.sqrt(gamma * gamma_c)
               + eta_bg, eta_max)


def tc_gamma_of_tau(tau, sigma_y, gamma_c, eta_bg):
    if tau <= sigma_y:
        return 0.0
    b = sigma_y / math.sqrt(gamma_c)
    disc = b * b - 4.0 * eta_bg * (sigma_y - tau)
    return ((-b + math.sqrt(disc)) / (2.0 * eta_bg)) ** 2


def tc_omega_of_T(T, Ri, Ro, L, sigma_y, gamma_c, eta_bg, n=2000):
    if T <= 0:
        return 0.0
    r_yield = math.sqrt(T / (2.0 * math.pi * L * sigma_y))
    r_hi = min(Ro, r_yield)
    if r_hi <= Ri:
        return 0.0
    ln_lo, ln_hi = math.log(Ri), math.log(r_hi)
    dln = (ln_hi - ln_lo) / n
    return sum(tc_gamma_of_tau(
        T / (2.0 * math.pi * L * math.exp(ln_lo + (i + 0.5) * dln) ** 2),
        sigma_y, gamma_c, eta_bg) * dln for i in range(n))


def tc_torque_of_omega(omega, Ri, Ro, L, sigma_y, gamma_c, eta_bg):
    T_lo = 2.0 * math.pi * L * sigma_y * Ri * Ri
    T_hi = T_lo * 2.0
    while tc_omega_of_T(T_hi, Ri, Ro, L, sigma_y, gamma_c, eta_bg) < omega:
        T_hi *= 2.0
    for _ in range(200):
        T_mid = 0.5 * (T_lo + T_hi)
        if tc_omega_of_T(T_mid, Ri, Ro, L, sigma_y, gamma_c, eta_bg) < omega:
            T_lo = T_mid
        else:
            T_hi = T_mid
        if (T_hi - T_lo) / T_mid < 1e-12:
            break
    return 0.5 * (T_lo + T_hi)


def lattice_params(case):
    g = case["geometry"]
    f = case["fluid"]
    s = case["sweep"]
    lat = case["lattice"]
    Ri, Ro, L = g["R_i"], g["R_o"], g["L"]
    gap = Ro - Ri
    rho = f["rho_kg_m3"]
    sy, gc, ebg = f["sigma_y_Pa"], f["gamma_dot_c_1_s"], f["eta_bg_Pa_s"]
    eta_max = f["eta_max_Pa_s"]
    N = lat["cells_across_gap"]
    dx = gap / N
    omega = lat["omega"]
    if not (0.5 < omega < 2.0):
        raise ValueError(f"omega={omega} outside (0.5, 2.0); unstable or overdamped")
    nu_lb = (1.0 / omega - 0.5) / 3.0
    points = {}
    for gamma in s["shear_rate_1_s"]:
        eta_label = tc_eta(gamma, sy, gc, ebg, eta_max)
        nu_label = eta_label / rho
        dt = nu_lb * dx * dx / nu_label
        run_time = s["diffusion_times"] * gap * gap / nu_label
        ramp_time = s["ramp_fraction"] * run_time
        n_steps = int(math.ceil(run_time / dt))
        Omega = gamma * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro)
        u_wall = Omega * Ri
        u_lb = u_wall * dt / dx
        T_an = tc_torque_of_omega(Omega, Ri, Ro, L, sy, gc, ebg)
        tau_i = T_an / (2.0 * math.pi * Ri * Ri * L)
        gamma_inner = tc_gamma_of_tau(tau_i, sy, gc, ebg)
        eta_inner = tc_eta(gamma_inner, sy, gc, ebg, eta_max)
        omega_inner = 1.0 / (3.0 * nu_lb * eta_inner / eta_label + 0.5)
        omega_plug = 1.0 / (3.0 * nu_lb * eta_max / eta_label + 0.5)
        r_yield = math.sqrt(T_an / (2.0 * math.pi * L * sy))
        points[str(gamma)] = {
            "gamma_1_s": gamma,
            "omega_rad_s": Omega,
            "eta_label_Pa_s": eta_label,
            "eta_inner_Pa_s": eta_inner,
            "gamma_inner_1_s": gamma_inner,
            "dt_s": dt,
            "run_time_s": run_time,
            "ramp_time_s": ramp_time,
            "n_steps": n_steps,
            "u_lb": u_lb,
            "omega_inner": omega_inner,
            "omega_plug": omega_plug,
            "T_analytic_N_m": T_an,
            "tau_inner_Pa": tau_i,
            "R_yield_mm": min(r_yield, Ro) * 1000,
            "yielded_throughout": r_yield >= Ro,
        }
    gamma_max = max(s["shear_rate_1_s"])
    return {
        "dx_m": dx,
        "gap_m": gap,
        "cells_across_gap": N,
        "omega_target": omega,
        "nu_lb": nu_lb,
        "eta_max_Pa_s": eta_max,
        "gamma_dot_reg_1_s": f["gamma_dot_reg_1_s"],
        "u_lb_max": max(p["u_lb"] for p in points.values()),
        "max_lattice_velocity": lat["max_lattice_velocity"],
        "mach_ok": max(p["u_lb"] for p in points.values()) < lat["max_lattice_velocity"],
        "yielded_omega_ok": all(0.8 <= p["omega_inner"] <= 1.9
                                for p in points.values()),
        "points": points,
    }


def main():
    case = load_case()
    p = lattice_params(case)
    out = CASE.with_name("lattice.json")
    out.write_text(json.dumps(p, indent=2) + "\n")
    print(f"dx = {p['dx_m']*1e3:.3f} mm   omega_target = {p['omega_target']:.3f}   "
          f"nu_lb = {p['nu_lb']:.4f}   eta_max = {p['eta_max_Pa_s']:.1f} Pa.s")
    print(f"{'gamma':>6} {'eta_lab':>9} {'dt':>10} {'steps':>6} "
          f"{'runTime':>9} {'u_lb':>7} {'om_in':>6} {'om_plug':>8} {'R_y_mm':>7}")
    for k, q in p["points"].items():
        print(f"{q['gamma_1_s']:>6} {q['eta_label_Pa_s']:9.3f} {q['dt_s']:10.3e} "
              f"{q['n_steps']:6d} {q['run_time_s']:9.3e} {q['u_lb']:7.4f} "
              f"{q['omega_inner']:6.2f} {q['omega_plug']:8.4f} {q['R_yield_mm']:7.2f}")
    if not p["mach_ok"]:
        print(f"FAIL: u_lb,max {p['u_lb_max']:.4f} >= "
              f"{p['max_lattice_velocity']} — reduce dt or coarsen dx")
        return 1
    if not p["yielded_omega_ok"]:
        print("FAIL: inner-wall omega outside [0.8, 1.9] at some point")
        return 1
    print("lattice check PASSED (yielded-region omega in [0.8, 1.9]; "
          "plug omega ~0.02 by cap design) — wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
