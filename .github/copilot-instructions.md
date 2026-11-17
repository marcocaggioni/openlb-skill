# Copilot instructions — openlb-skill

This repo teaches an agent to run open-source OpenLB lattice
Boltzmann simulations from an English discussion. The skill is
`skills/openlb-simulation/SKILL.md`, mirrored at
`.github/skills/openlb-simulation/SKILL.md`.

Rules for every simulation in this repo:

1. Discuss in English first, then write `case.yaml`, show it, and
   get a yes before building. Never invent a number.
2. Convert physical → lattice units with `lattice.py` and check
   the low-Mach condition before running.
3. Default to a dry run (lattice check + build only). Never run the
   solver without saying so first.
4. Validate a Newtonian case against the analytic Couette torque
   (within ~10%) before any non-Newtonian flow curve.
5. OpenLB has no rheology dropdown and no GUI: non-Newtonian
   means a local relaxation frequency, and the C++ case is a
   program you compile with GNU Make, not a file you open.
6. The Couette C++ case follows OpenLB's own examples
   (`forBeginners/cylinder2d`, `laminar/cavity2d`,
   `fsi/taylorCouette3d`) — check `examples/couette/RUN_REPORT.md`
   for its run status before trusting it.
