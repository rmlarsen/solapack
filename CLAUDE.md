# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SOLApack is a Fortran 77 / C scientific computing package that implements the SOLA (Subtractive Optimally Localized Averages) inversion algorithm for helioseismology. It infers the 2-D solar rotation rate from measured rotational frequency splittings of solar oscillation eigenmodes.

## Build Commands

```bash
# Configure and build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Specify BLAS vendor (optional)
cmake -B build -DBLA_VENDOR=OpenBLAS
```

Dependencies: Fortran compiler, C compiler, CMake >= 3.15, LAPACK/BLAS.

## Running Tests

```bash
ctest --test-dir build        # run small test after building
ctest --test-dir build -V     # verbose output
```

## Running the Example

```bash
cd example
./runtest.sh                   # small test (25 modes, l=0,1,2)
SOLAPACK_FULL=1 ./runtest.sh   # full run (29617 modes)
./runtest.sh /path/to/build    # alternate build directory
```

Programs read parameters from stdin. The `get-config` script strips `#`-comment lines from `.cfg` files before piping to executables.

## Architecture

### Two executables

- **`set-2drls`** — Kernel setup: reads mode parameters + eigenfunction binary files, computes and writes mode kernels. Entry point: `src/set-2drls.nn.f` → calls `src/set2dker.nn.f`, `src/kersub.n.f`.
- **`2dsola`** — SOLA inversion: reads splittings + kernels, solves linear inverse problem via SOLA with Lanczos bidiagonalization. Entry point: `src/2dsola.main.lanczos.f` with C driver `src/2dsola.lanczos.c`.

### Shared code (`solapack_common` OBJECT library)

I/O (`2dolaIO.f`, `byteswapIO.c`), target functions (`target.f`), averaging kernels (`avker.f`), covariance (`covariance.f`), BLAS helpers (`parblas.f`, `dgemm_ovwr.f`), disk allocation (`diskmalloc.c`), timing (`timer.c`), and error handling (`xerbla.f`).

### `linpack/` — Static library

Bundled LINPACK routines (LU and Cholesky factorization/solve) used by both executables.

### Build option

`SOLAPACK_USE_RISC_MGS` selects between `mgs.pentium.f` (default) and `mgs.risc.f` for the Modified Gram-Schmidt implementation in the Lanczos solver.

### Data files

`modeldata/` contains little-endian binary eigenfunction and solar model files. Byte-swapping for big-endian is controlled by `ibyteswap` in config files.

### Configuration

`.cfg` files in `example/` control all runtime parameters. `*.small.cfg` variants use a 25-mode subset for quick testing.
