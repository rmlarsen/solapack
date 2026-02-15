cmake_minimum_required(VERSION 3.15)
project(solapack
  VERSION 1.0.0
  DESCRIPTION "SOLA inversion of helioseismic rotational splitting data"
  LANGUAGES C Fortran
)

# Use standard install directories (lib, bin, etc.)
include(GNUInstallDirs)

# Options
option(SOLAPACK_USE_RISC_MGS
  "Use RISC-optimized Modified Gram-Schmidt (instead of Pentium version)" OFF)

# Find LAPACK (also provides BLAS)
find_package(LAPACK REQUIRED)

# Fortran include files are not compiled directly; mark them header-only
# so they appear in IDE project views and trigger rebuilds when changed.
set(FORTRAN_INCLUDES
  src/2dsola.n.include.f
  src/include1.f
  src/include1.n.f
  src/include2.f
)
set_source_files_properties(${FORTRAN_INCLUDES} PROPERTIES HEADER_FILE_ONLY TRUE)

# --- Linpack static library ---
add_library(linpack STATIC
  linpack/dgeco.f
  linpack/dgedi.f
  linpack/dgefa.f
  linpack/dpoco.f
  linpack/dpofa.f
  linpack/dposl.f
)

# --- Common object library (shared by both executables) ---
add_library(solapack_common OBJECT
  src/elden.f
  src/timer.c
  src/dgemm_ovwr.f
  src/byteswapIO.c
  src/parblas.f
  src/diskmalloc.c
  src/xerbla.f
  src/2dolaIO.f
  src/target.f
  src/avker.f
  src/covariance.f
  src/compat.c
  src/2dsola.n.include.f
)
target_include_directories(solapack_common PRIVATE src)

# --- Kernel setup executable: set-2drls ---
add_executable(set-2drls
  src/set-2drls.nn.f
  src/set2dker.nn.f
  src/kersub.n.f
  src/splitsub.f
  src/miscsub.f
  src/ran.f
  src/plm.f
  src/include1.f
  src/include1.n.f
  src/include2.f
  $<TARGET_OBJECTS:solapack_common>
)
target_include_directories(set-2drls PRIVATE src)
target_link_libraries(set-2drls PRIVATE linpack LAPACK::LAPACK)

# --- 2D SOLA inversion executable: 2dsola ---
if(SOLAPACK_USE_RISC_MGS)
  set(MGS_SOURCE src/mgs.risc.f)
else()
  set(MGS_SOURCE src/mgs.pentium.f)
endif()

add_executable(2dsola
  src/2dsola.lanczos.c
  src/2dsola.main.lanczos.f
  src/2dsola.sub.lanczos.f
  src/dlanc_b.f
  src/dkronecker2.f
  ${MGS_SOURCE}
  $<TARGET_OBJECTS:solapack_common>
)
target_include_directories(2dsola PRIVATE src)
target_link_libraries(2dsola PRIVATE linpack LAPACK::LAPACK)

# --- Install rules ---
install(TARGETS set-2drls 2dsola
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY example/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/solapack/example
)
