# Quickstart — openlb-skill

```bash
git clone https://github.com/marcocaggioni/openlb-skill.git
cd openlb-skill
uv sync            # Python side: numpy, matplotlib, pyyaml (see note below)
bash scripts/check-setup.sh
```

OpenLB itself is C++, not Python — there is no PyPI package for
it, so `uv sync` cannot install it. Fetch and build it separately:

```bash
bash scripts/install-openlb.sh   # clones OpenLB to ~/openlb
export OPENLB_ROOT=~/openlb
```

Then build the Couette case and run it:

```bash
cd examples/couette && make
./../../run-couette.sh --run
```

(Use plain `make`, one job: parallel builds of OpenLB's
template-heavy core can exhaust memory on small machines. The first
build compiles the core library — a few minutes, once.)

Then start Copilot CLI in the repo root and invoke the skill:

```bash
copilot
```

```
/openlb-simulation run the Couette example in examples/couette/case.yaml
```

Or in one go:

```bash
copilot -i "Use the /openlb-simulation skill to run the Couette example in examples/couette/case.yaml"
```

Copilot must be started from the repo root for skill discovery.
If the skill was added during a session, run `/skills reload` first.
The skill has no `allowed-tools` frontmatter on purpose — Copilot
must ask before running shell commands or a solver run.
