#!/usr/bin/env bash
# Fetch the OpenLB C++ library (there is no PyPI package for it) and
# prove the toolchain by building one of OpenLB's own examples.
# Usage:
#   bash scripts/install-openlb.sh              # clone to ~/openlb
#   OPENLB_ROOT=/somewhere bash scripts/install-openlb.sh
#   bash scripts/install-openlb.sh --build-example   # also build cavity2d
# After this: export OPENLB_ROOT=~/openlb
# Then build this repo's Couette case:
#   cd examples/couette && make
# Use plain `make` (one job): parallel builds of the template-heavy
# core can exhaust memory on small machines.
set -e
ROOT="${OPENLB_ROOT:-$HOME/openlb}"

if [ -d "$ROOT/.git" ]; then
  echo "OpenLB already cloned at $ROOT — pulling latest."
  git -C "$ROOT" pull --ff-only
else
  echo "Cloning OpenLB into $ROOT ..."
  git clone --depth 1 https://gitlab.com/openlb/release.git "$ROOT"
fi

if [ "$1" = "--build-example" ]; then
  echo "Building OpenLB's own cavity2d example to prove the toolchain ..."
  cd "$ROOT/examples/laminar/cavity2d"
  make
  echo "cavity2d built OK."
fi

echo "OpenLB is at $ROOT"
echo "Add to your shell: export OPENLB_ROOT=$ROOT"
