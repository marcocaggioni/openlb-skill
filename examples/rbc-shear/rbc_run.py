#!/usr/bin/env python3
"""Run the RBC single-cell shear test.

    python3 rbc_run.py            # dry run: print the plan
    python3 rbc_run.py --run      # build (if needed) + run the case

The shipped example binary reads parametersMeshes.xml,
parametersPositions.xml and sph_ico_1280.msh from its working
directory and writes ./tmp/ there, so the run happens in
runs/single/ (inputs symlinked in). Serial only — the example
undefines PARALLEL_MODE_MPI. Then post-process with rbc_post.py.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
BIN = HERE / "rbc-shear"
RUN_DIR = HERE / "runs" / "single"
INPUTS = ["parametersMeshes.xml", "parametersPositions.xml", "sph_ico_1280.msh"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="actually run the case")
    args = ap.parse_args()

    print("single deformable capsule in shear flow (Krüger et al. 2011):")
    print("  domain 50^3, sphere r=5, G=0.01, Re=0.02, 9000 steps to kappa*t=120")
    print("  outputs: 16 snapshots (every kappa*t = 8) in runs/single/tmp/")
    if not args.run:
        print("DRY RUN: pass --run to execute. Binary expected at", BIN)
        return 0
    if not BIN.exists():
        print(f"Binary not found: {BIN} — build it first (cd examples/rbc-shear && make)")
        return 1

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    for name in INPUTS:
        link = RUN_DIR / name
        if not link.exists():
            link.symlink_to(HERE / name)
    t0 = time.time()
    print("running:", BIN, "(cwd runs/single)", flush=True)
    r = subprocess.run([str(BIN)], cwd=RUN_DIR)
    dt = time.time() - t0
    print(f"run done in {dt:.1f} s (exit {r.returncode}) — "
          f"post-process with: python3 rbc_post.py")
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
