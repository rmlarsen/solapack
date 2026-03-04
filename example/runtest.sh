#!/bin/sh
#
# Run the SOLApack example: kernel setup followed by 2D SOLA inversion.
#
# Usage: cd example && ./runtest.sh [path-to-build-dir]
#
# The build directory defaults to ../build.
#
# Note: the full example requires a complete eigenfunction file (amde.1)
# covering all modes in the splitting data. By default we run a small
# test using the subset of modes available in the distributed amde.1
# (l=0,1,2 only). Set SOLAPACK_FULL=1 to attempt the full run.

BUILDDIR="${1:-../build}"

mkdir -p Kernels Inversion

if [ "${SOLAPACK_FULL:-0}" = "1" ]; then
  SETUP_CFG=set-2drls.cfg
  INVERSION_CFG=2dsola.cfg
else
  SETUP_CFG=set-2drls.small.cfg
  INVERSION_CFG=2dsola.small.cfg
fi

echo "********* Setting up mode kernels ($SETUP_CFG) **********"
../get-config "$SETUP_CFG" | "$BUILDDIR/set-2drls"

echo "********* Performing inversion ($INVERSION_CFG) **********"
../get-config "$INVERSION_CFG" | "$BUILDDIR/2dsola"
