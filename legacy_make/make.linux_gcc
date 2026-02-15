###  Compiler flags for fortran and C modules.

#  Platform postfix to use in naming executables and library files
PLAT = LINUX_GCC

# Architecture specific optimization flags
# Pentium 4 flags
P4_OPT = -O3 -funroll-loops -march=pentium4 -mfpmath=sse
# Pentium 2 flags
P2_OPT =  -O3 -funroll-loops -march=pentium2


# Select which set of optimization flags to use
OPT = $(P4_OPT)
#OPT = -O0 
NO_OPT = -O0

# Profiling flags
PROFILE = -Wall -pedantic

# Compiler flags
FFLAGS = -g -fno-second-underscore -fno-f2c $(OPT) $(PROFILE)
FFLAGS_NOOPT = -g  -fno-second-underscore -fno-f2c $(NOOPT) $(PROFILE)
CFLAGS = -g $(OPT) $(PROFILE)
LINKFLAGS = -g  $(OPT) $(PROFILE) 

# Which compiler and linker to use
CC = gcc
F77 = g77 
LINKER = g77 

# Linker paths and flags
LIBPATH =  -L$(HOME)/lib/$(ARCH)
LINPACK = -llinpack
LAPACK = -llapack
#LAPACK = -lmkl_lapack -lguide # Intel MKL library
#BLAS = -lblas
#BLAS = -lmkl -lguide  # Intel MKL library
BLAS = -lgoto_p4_512-r0.6
#BLAS = -lf77blas -latlas # Atlas BLAS
RANLIB = ranlib
MGS = mgs.pentium


# MPI scripts
MPIF77 = mpif77
MPICC = mpicc
MPILINKER = mpif77
MPILIB =
