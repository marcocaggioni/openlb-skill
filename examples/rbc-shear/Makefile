# RBC single-cell shear test — OpenLB case (GNU Make build, like OpenLB's own
# examples). This is OpenLB's shipped fsi/deformableParticleInShearFlow3d
# benchmark (Krüger et al. 2011), copied here with a run/post wrapper.
#
# Build (serial, no MPI — the example undefines PARALLEL_MODE_MPI):
#   cd examples/rbc-shear && make
# Use plain `make` (one job): parallel builds of the template-heavy core
# can exhaust memory on small machines.
#
# OPENLB_ROOT: where OpenLB lives (default ~/openlb). Override with:
#   make OPENLB_ROOT=/path/to/openlb
EXAMPLE = rbc-shear
OPENLB_ROOT ?= $(HOME)/openlb
OLB_ROOT := $(OPENLB_ROOT)
include $(OLB_ROOT)/default.mk
