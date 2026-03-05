#!/bin/sh
#
# CTest driver for the SOLApack small test.
#
# Usage: run_small_test.sh <build-dir> <source-dir>
#
# Runs both steps (kernel setup + inversion) and checks:
#   1. Both executables exit successfully
#   2. Expected output files are produced and non-empty
#   3. stdout contains expected diagnostic strings
#   4. The inversion reports "Normal termination"
#
set -e

BUILDDIR="$1"
SRCDIR="$2"
EXDIR="$SRCDIR/example"
GET_CONFIG="$SRCDIR/get-config"

# Work in a temporary directory so tests don't clobber the source tree.
# Config files reference ../modeldata relative to the working directory,
# so we create a parent dir with a modeldata symlink and run inside a subdir.
TMPROOT=$(mktemp -d)
trap 'rm -rf "$TMPROOT"' EXIT

WORKDIR="$TMPROOT/example"
mkdir -p "$WORKDIR/Kernels" "$WORKDIR/Inversion"

# Symlink modeldata so ../modeldata resolves from WORKDIR
ln -s "$SRCDIR/modeldata" "$TMPROOT/modeldata"

# Copy input data needed by both steps
cp -r "$EXDIR/Splittings" "$WORKDIR/"
cp "$EXDIR/set2d.cfg" "$EXDIR/efunc2d.cfg" "$WORKDIR/"
cp "$EXDIR/set-2drls.small.cfg" "$EXDIR/2dsola.small.cfg" "$WORKDIR/"

cd "$WORKDIR"

FAIL=0

# ── Step 1: kernel setup ────────────────────────────────────────────
echo "=== Step 1: set-2drls (kernel setup) ==="
SETUP_OUT=$("$GET_CONFIG" set-2drls.small.cfg | "$BUILDDIR/set-2drls" 2>&1) || {
    echo "FAIL: set-2drls exited with non-zero status"
    echo "$SETUP_OUT"
    FAIL=1
}

if [ $FAIL -eq 0 ]; then
    # Check expected kernel output files exist and are non-empty
    for f in Kernels/As.small.10151 Kernels/rthet.small.10151 Kernels/Sum.small.10151; do
        if [ ! -s "$f" ]; then
            echo "FAIL: expected output file $f missing or empty"
            FAIL=1
        fi
    done

    # Check diagnostic output
    if echo "$SETUP_OUT" | grep -q "No of modes read in:"; then
        NMODES=$(echo "$SETUP_OUT" | sed -n 's/.*No of modes read in: *\([0-9]*\).*/\1/p')
        if [ "$NMODES" != "25" ]; then
            echo "FAIL: expected 25 modes, got $NMODES"
            FAIL=1
        else
            echo "OK: 25 modes read"
        fi
    else
        echo "FAIL: missing 'No of modes read in' in set-2drls output"
        FAIL=1
    fi

    if echo "$SETUP_OUT" | grep -q "f matrix has been set up"; then
        echo "OK: f matrix set up"
    else
        echo "FAIL: f matrix not set up"
        FAIL=1
    fi

    if echo "$SETUP_OUT" | grep -q "g matrix has been set up"; then
        echo "OK: g matrix set up"
    else
        echo "FAIL: g matrix not set up"
        FAIL=1
    fi
fi

# ── Step 2: 2D SOLA inversion ──────────────────────────────────────
echo ""
echo "=== Step 2: 2dsola (inversion) ==="
INV_OUT=$("$GET_CONFIG" 2dsola.small.cfg | "$BUILDDIR/2dsola" 2>&1) || {
    echo "FAIL: 2dsola exited with non-zero status"
    echo "$INV_OUT"
    FAIL=1
}

if [ $FAIL -eq 0 ]; then
    # Check expected inversion output files exist and are non-empty
    for f in Inversion/sol.small Inversion/avker.small; do
        if [ ! -s "$f" ]; then
            echo "FAIL: expected output file $f missing or empty"
            FAIL=1
        fi
    done

    # Must report normal termination
    if echo "$INV_OUT" | grep -q "Normal termination"; then
        echo "OK: normal termination"
    else
        echo "FAIL: missing 'Normal termination' in 2dsola output"
        FAIL=1
    fi

    # Check that Lanczos iterations ran
    if echo "$INV_OUT" | grep -q "Total time in dlanc_b"; then
        echo "OK: Lanczos bidiagonalization completed"
    else
        echo "FAIL: Lanczos bidiagonalization did not complete"
        FAIL=1
    fi

    # Verify number of target parameters
    NTARGETS=$(echo "$INV_OUT" | sed -n 's/.*[^0-9]\([0-9][0-9]*\) *sets of target parameters.*/\1/p')
    if [ "$NTARGETS" = "5151" ]; then
        echo "OK: 5151 target parameters"
    else
        echo "FAIL: expected 5151 target parameters, got '$NTARGETS'"
        FAIL=1
    fi

    # Verify M_kers = 25
    MKERS=$(echo "$INV_OUT" | grep 'M_kers' | head -1 | sed -n 's/.*M_kers[^=]*= *\([0-9]*\).*/\1/p')
    if [ "$MKERS" = "25" ]; then
        echo "OK: M_kers = 25"
    else
        echo "FAIL: expected M_kers = 25, got '$MKERS'"
        FAIL=1
    fi

    # Check first rotation rate value (should be ~435.126)
    ROT1=$(echo "$INV_OUT" | sed -n 's/.*rot= *\([0-9.]*\).*/\1/p')
    if [ -n "$ROT1" ]; then
        # Check that it starts with 435.1 (rough sanity check)
        case "$ROT1" in
            435.1*)
                echo "OK: first rotation rate = $ROT1"
                ;;
            *)
                echo "FAIL: unexpected first rotation rate = $ROT1 (expected ~435.12)"
                FAIL=1
                ;;
        esac
    else
        echo "FAIL: could not extract rotation rate from output"
        FAIL=1
    fi
fi

echo ""
if [ $FAIL -ne 0 ]; then
    echo "FAILED"
    exit 1
fi
echo "ALL CHECKS PASSED"
exit 0
