#!/usr/bin/env python3
"""Physical -> lattice parameters for the OpenLB cavity-TC case.

Reads examples/cavity-tc/case.yaml, solves dt from omega = 1.0 at the
characteristic viscosity eta_char = eta_TC(U_lid/L), and checks the
low-Mach condition. Also holds the Ghia, Ghia & Shin (1982) Re = 100
centreline reference used to validate the Newtonian run of the same
binary (copied from the values embedded in OpenLB's own
laminar/cavity2d/cavity2d.cpp, which cites Ghia et al. Table 1).

Formulas (D2Q9, cs^2 = 1/3): same as examples/couette-tc/lattice.py.
    dx     = L / N
    dt     = nu_lb * dx^2 / nu_char,  nu_lb = (1/omega - 0.5)/3
    u_lb   = U_lid * dt / dx          ->  must stay < 0.1
    run_time / diffusion time: L^2/nu_char = 1.08 s at eta_char = 9.26
"""
import json
import math
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

CASE = Path(__file__).with_name("case.yaml")

# Ghia, Ghia & Shin (1982), J. Comput. Phys. 48, 387–411, Table 1,
# u_x along the vertical centreline at Re = 100 (as embedded in
# OpenLB's laminar/cavity2d/cavity2d.cpp). y is the fraction y/L.
GHIA_RE100_Y = [128, 125, 124, 123, 122, 109, 94, 79, 64, 58, 36, 22,
                13, 9, 8, 7, 0]
GHIA_RE100_U = [1.0, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151,
                0.00332, -0.13641, -0.20581, -0.21090, -0.15662,
                -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.0]


def load_case():
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    return yaml.safe_load(CASE.read_text())


def tc_eta(gamma, sigma_y, gamma_c, eta_bg, eta_max):
    if gamma <= 0:
        return eta_max
    return min(sigma_y / gamma + sigma_y / math.sqrt(gamma * gamma_c)
               + eta_bg, eta_max)


def lattice_params(case):
    g = case["geometry"]
    f = case["fluid"]
    fl = case["flow"]
    lat = case["lattice"]
    L = g["L"]
    rho = f["rho_kg_m3"]
    sy, gc, ebg = f["sigma_y_Pa"], f["gamma_dot_c_1_s"], f["eta_bg_Pa_s"]
    eta_max = f["eta_max_Pa_s"]
    U = fl["U_lid_m_s"]
    N = lat["N"]
    dx = L / N
    omega = lat["omega"]
    nu_lb = (1.0 / omega - 0.5) / 3.0
    gamma_char = U / L
    eta_char = tc_eta(gamma_char, sy, gc, ebg, eta_max)
    nu_char = eta_char / rho
    dt = nu_lb * dx * dx / nu_char
    u_lb = U * dt / dx
    n_steps = int(math.ceil(fl["run_time_s"] / dt))
    eta_newt = case["validation"]["newtonian_eta_Pa_s"]
    # Newtonian check runs at omega=1.5, not 1.0: at omega=1.0 its dt
    # gives u_lb = 0.13 (Ma = 0.23) — too compressible for Ghia.
    omega_newt = case["validation"].get("newtonian_omega", 1.5)
    nu_lb_newt = (1.0 / omega_newt - 0.5) / 3.0
    dt_newt = nu_lb_newt * dx * dx / (eta_newt / rho)
    return {
        "L_m": L, "U_lid_m_s": U, "N": N, "dx_m": dx,
        "omega_target": omega, "nu_lb": nu_lb,
        "gamma_char_1_s": gamma_char,
        "tau_char_Pa": sy + sy * math.sqrt(gamma_char / gc) + ebg * gamma_char,
        "eta_char_Pa_s": eta_char,
        "eta_max_Pa_s": eta_max,
        "dt_s": dt, "u_lb": u_lb,
        "max_lattice_velocity": lat["max_lattice_velocity"],
        "mach_ok": u_lb < lat["max_lattice_velocity"],
        "run_time_s": fl["run_time_s"],
        "n_steps": n_steps,
        "diffusion_time_s": L * L / nu_char,
        "snap_interval_s": fl["snapshot_interval_s"],
        "newtonian": {
            "eta_Pa_s": eta_newt,
            "omega": omega_newt,
            "u_lb": U * dt_newt / dx,
            "Re": U * L / (eta_newt / rho),
            "dt_s": dt_newt,
            "n_steps": int(math.ceil(fl["run_time_s"] / dt_newt)),
        },
    }


def main():
    case = load_case()
    p = lattice_params(case)
    out = CASE.with_name("lattice.json")
    out.write_text(json.dumps(p, indent=2) + "\n")
    print(f"L = {p['L_m']} m  U_lid = {p['U_lid_m_s']} m/s  N = {p['N']}  "
          f"dx = {p['dx_m']*1e3:.3f} mm")
    print(f"char rate U/L = {p['gamma_char_1_s']} 1/s  tau_char = "
          f"{p['tau_char_Pa']:.1f} Pa  (> 70 Pa: yielded at the lid)")
    print(f"eta_char = {p['eta_char_Pa_s']:.3f} Pa.s  dt = {p['dt_s']:.3e} s  "
          f"u_lb = {p['u_lb']:.4f}  steps = {p['n_steps']:,}  "
          f"t_diff = {p['diffusion_time_s']:.2f} s")
    print(f"Newtonian check: eta = {p['newtonian']['eta_Pa_s']} Pa.s  "
          f"omega = {p['newtonian']['omega']}  Re = {p['newtonian']['Re']:.0f}  "
          f"dt = {p['newtonian']['dt_s']:.3e} s  u_lb = {p['newtonian']['u_lb']:.4f}  "
          f"steps = {p['newtonian']['n_steps']:,}")
    if not p["mach_ok"]:
        print("FAIL: u_lb too high")
        return 1
    print("lattice check PASSED — wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
