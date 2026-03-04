# SOLApack

SOLApack is a Fortran 77 software package that implements the Subtractive
Optimally Localized Averages (SOLA) inversion algorithm. It infers the
symmetric component of the 2-D solar rotation rate from measured rotational
frequency splittings of the eigenmodes of global solar oscillation.

**Authors:** Rasmus Munk Larsen (Stanford University), Jesper Schou
(Stanford University) & Jorgen Christensen-Dalsgaard (Aarhus University).

**Version:** 1.0, August 2003

## Programs

SOLApack builds two executables:

1. **`set-2drls`** - Sets up mode kernels given a list of mode parameters
   and a file with the mode eigenfunctions calculated from a solar model.

2. **`2dsola`** - Reads rotational splittings and mode kernels, then solves
   a linear inverse problem using the SOLA algorithm with Lanczos
   bidiagonalization to infer the symmetric component of the solar rotation
   rate as a function of radius and latitude.

## Dependencies

- A Fortran compiler (gfortran, ifort, etc.)
- A C compiler (gcc, clang, etc.)
- CMake >= 3.15
- LAPACK and BLAS libraries

## Building with CMake

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `SOLAPACK_USE_RISC_MGS` | `OFF` | Use RISC-optimized Modified Gram-Schmidt |

### Specifying LAPACK/BLAS

CMake's `find_package(LAPACK)` will locate your system LAPACK/BLAS. To use
a specific installation (e.g. Intel MKL, OpenBLAS), set `BLA_VENDOR`:

```bash
cmake -B build -DBLA_VENDOR=OpenBLAS
cmake -B build -DBLA_VENDOR=Intel10_64lp  # MKL
```

## Running the example

The programs read parameters from stdin. Configuration files use `#`
for comment lines; the `get-config` script strips them before piping
to the executables.

```bash
cd example
./runtest.sh                   # small test (25 modes, l=0,1,2 only)
SOLAPACK_FULL=1 ./runtest.sh   # full run (29617 modes, ~22 seconds)
```

You can pass an alternate build directory: `./runtest.sh /path/to/build`.

## Data files

The `modeldata/` directory contains little-endian binary files:

- **`amde.1`** — Mode eigenfunctions from the l5bi.d.15 solar model
  (4,644 modes, l=0–300, n=0–35).
- **`amdl.l5bi.d.15`** — Solar model for sound-speed scaling.

## Data formats

Most data files used by SOLApack are stored in raw binary format and
differ between little-endian (Intel, Alpha) and big-endian (Sun, SGI, IBM)
machines. The kernel setup program reads and writes files in the native
format. The inversion program can read and write files with the opposite
endianness by setting `ibyteswap = 1` in the configuration file.

## License

Copyright Rasmus Munk Larsen, Stanford University, 2003.
