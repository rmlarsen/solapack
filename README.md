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

See the `example/` directory and the original [README](README.orig) for
details on the configuration files and data formats. The basic workflow is:

```bash
# 1. Set up kernel matrix from mode eigenfunctions
./build/set-2drls < example/set-2drls.cfg

# 2. Perform 2D SOLA inversion
./build/2dsola < example/2dsola.cfg
```

The programs read their parameters from stdin via configuration files.
Edit the paths in the `.cfg` files to point to your eigenfunction and
model data files before running.

## Data files

The eigenfunction (`amde.1`) and model (`amdl.l5bi.d.15`) data files
are binary and are available in both big-endian (Sun/SPARC, SGI/MIPS)
and little-endian (Intel, Alpha) formats. These files were originally
hosted at `http://sun.stanford.edu/~rmunk/SOLApack/` but are no longer
available from that location.

## Data formats

Most data files used by SOLApack are stored in raw binary format and
differ between little-endian (Intel, Alpha) and big-endian (Sun, SGI, IBM)
machines. The kernel setup program reads and writes files in the native
format. The inversion program can read and write files with the opposite
endianness by setting `ibyteswap = 1` in the configuration file.

## License

Copyright Rasmus Munk Larsen, Stanford University, 2003.
