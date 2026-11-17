# Couette example — OpenLB case (GNU Make build, like OpenLB's own examples).
#
# Build (serial, no MPI):
#   cd examples/couette && make
# The first build also compiles OpenLB's core library (libolbcore.a) and
# its embedded dependencies — a few minutes, once. Use plain `make`
# (one job): parallel builds of the template-heavy core can exhaust
# memory on small machines.
#
# OPENLB_ROOT: where OpenLB lives (default ~/openlb). Override with:
#   make OPENLB_ROOT=/path/to/openlb
EXAMPLE = couette
OPENLB_ROOT ?= $(HOME)/openlb
OLB_ROOT := $(OPENLB_ROOT)
include $(OLB_ROOT)/default.mk
