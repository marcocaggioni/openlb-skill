#!/usr/bin/env python3
"""Sweep the Couette-TC case over the shear rates in case.yaml.

For each gamma: run the compiled OpenLB binary (./couette-tc, built in
examples/couette-tc/ by `make`), and collect
runs/gamma_<value>/torque.csv. Default is a dry run: print the plan
and the per-point lattice parameters, run nothing. Serial only.

Usage:
    python3 couette_sweep.py            # dry run
    python3 couette_sweep.py --run      # actually run the sweep
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from lattice import lattice_params, load_case  # noqa: E402

HERE = Path(__file__).parent
BIN = HERE / "couette-tc"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="actually run the solver")
    args = ap.parse_args()

    case = load_case()
    lat = lattice_params(case)
    if not lat["mach_ok"] or not lat["yielded_omega_ok"]:
        print("FAIL: lattice check failed — fix case.yaml first")
        return 1
    for k, q in lat["points"].items():
        print(f"gamma={q['gamma_1_s']:>6}  dt={q['dt_s']:.3e} s  "
              f"steps={q['n_steps']:,}  runTime={q['run_time_s']:.3e} s  "
              f"u_lb={q['u_lb']:.4f}  omega_inner={q['omega_inner']:.2f}  "
              f"R_y={q['R_yield_mm']:.2f} mm")

    gammas = case["sweep"]["shear_rate_1_s"]
    if not args.run:
        print(f"DRY RUN: would run {len(gammas)} points: {gammas} 1/s")
        print("Pass --run to execute. Binary expected at", BIN)
        return 0
    if not BIN.exists():
        print(f"Binary not found: {BIN} — build it first (cd examples/couette-tc && make)")
        return 1

    import time
    t0 = time.time()
    for gamma in gammas:
        out = HERE / "runs" / f"gamma_{gamma}"
        out.mkdir(parents=True, exist_ok=True)
        cmd = [str(BIN), "--gamma", str(gamma), "--out", str(out)]
        print("running:", " ".join(cmd), flush=True)
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"run failed for gamma={gamma}")
            return r.returncode
        print(f"  gamma={gamma} done in {time.time()-t0:.1f} s total", flush=True)
    print(f"sweep done in {time.time()-t0:.1f} s — post-process with: python3 couette_post.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
