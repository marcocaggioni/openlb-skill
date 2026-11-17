#!/usr/bin/env python3
"""Run the cavity case: Newtonian validation first, then the TC run.

    python3 cavity_run.py              # dry run: print the plan
    python3 cavity_run.py --run        # run Newtonian (Re=100) + TC
    python3 cavity_run.py --run --tc-only
    python3 cavity_run.py --run --newtonian-only

Each run writes runs/<mode>/{energy.csv, centerline.csv, field_*.csv,
meta.json}. Serial only. Then post-process with cavity_post.py.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from cavity import lattice_params, load_case  # noqa: E402

HERE = Path(__file__).parent
BIN = HERE / "cavity"


def run_one(cmd, label):
    t0 = time.time()
    print("running:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    dt = time.time() - t0
    print(f"  {label} done in {dt:.1f} s (exit {r.returncode})", flush=True)
    return r.returncode, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--tc-only", action="store_true")
    ap.add_argument("--newtonian-only", action="store_true")
    args = ap.parse_args()

    case = load_case()
    lat = lattice_params(case)
    if not lat["mach_ok"]:
        print("FAIL: lattice check failed")
        return 1
    print(f"TC: dt={lat['dt_s']:.3e} s steps={lat['n_steps']:,} "
          f"runTime={lat['run_time_s']} s u_lb={lat['u_lb']:.4f}")
    print(f"Newtonian: dt={lat['newtonian']['dt_s']:.3e} s "
          f"steps={lat['newtonian']['n_steps']:,} Re={lat['newtonian']['Re']:.0f}")
    if not args.run:
        print("DRY RUN: pass --run to execute. Binary expected at", BIN)
        return 0
    if not BIN.exists():
        print(f"Binary not found: {BIN} — build it first (cd examples/cavity-tc && make)")
        return 1

    t00 = time.time()
    if not args.tc_only:
        out = HERE / "runs" / "newtonian"
        out.mkdir(parents=True, exist_ok=True)
        rc, _ = run_one([str(BIN), "--newtonian", "--eta", "1.0",
                         "--omega", "1.5",
                         "--out", str(out)], "newtonian Re=100")
        if rc:
            return rc
    if not args.newtonian_only:
        out = HERE / "runs" / "tc"
        out.mkdir(parents=True, exist_ok=True)
        rc, _ = run_one([str(BIN), "--out", str(out)], "TC")
        if rc:
            return rc
    print(f"all runs done in {time.time()-t00:.1f} s — post-process with: "
          f"python3 cavity_post.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
