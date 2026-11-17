# AGENTS.md — openlb-skill

Open-source counterpart of the private `mstar-skill` repo: same
discussion → case → run → report workflow, same Couette example,
OpenLB solver instead of licensed M-Star. No license, no GUI.

- Skill: `skills/openlb-simulation/SKILL.md`, mirrored byte-for-byte
  at `.github/skills/openlb-simulation/SKILL.md` (the copy Copilot
  CLI discovers). Keep the two copies in sync; edit both per commit.
- No `allowed-tools` in the skill frontmatter, deliberately — Copilot
  must ask before shell/solver runs.
- `case.yaml` is the contract. Defaults are labelled; confirm before
  running. Never invent a number.
- Lattice units are the agent's job: `lattice.py` converts physical
  → lattice and enforces u_lb < 0.1. Never run with omega outside
  roughly [0.8, 1.9] without saying so.
- Validation gate: Newtonian torque vs analytic Couette torque within
  ~10% before any non-Newtonian curve. `couette_post.py` enforces it.
- Build is GNU Make (`make` in examples/couette/), not CMake.
  OpenLB lives at $OPENLB_ROOT (default ~/openlb). Use plain `make`
  (one job) — parallel core builds can exhaust memory on small machines.
- v0.1 status: `lattice.py` and `couette_post.py` are tested here;
  `couette.cpp` compiled and run against OpenLB 1.9r0 (commit
  4b911c3) on 2026-11-17 — Newtonian validation PASSED (max error
  0.94%, η_app = 1.0093 Pa·s); see `examples/couette/RUN_REPORT.md`.
- OpenLB itself is GPLv2; this repo's own code is MIT.
