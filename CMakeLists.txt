cmake_minimum_required(VERSION 3.15)
project(solapack
  VERSION 1.0.0
  DESCRIPTION "SOLA inversion of helioseismic rotational splitting data"
  LANGUAGES C Fortran
)

# Use standard install directories (lib, bin, etc.)
include(GNUInstallDirs)

# Options
# Compiler warnings
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
  add_compile_options($<$<COMPILE_LANGUAGE:C>:-Wall>)
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Fortran 77 codebase: suppress warnings for intentional idioms
  add_compile_options(
    $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wno-conversion>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wno-surprising>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>
    $<$<COMPILE_LANGUAGE:Fortran>:-fno-automatic>
  )
endif()

# Find LAPACK (also provides BLAS).
# Prefer OpenBLAS if available, otherwise fall back to any LAPACK.
if(NOT BLA_VENDOR)
  set(BLA_VENDOR OpenBLAS)
  find_package(LAPACK QUIET)
  if(NOT LAPACK_FOUND)
    unset(BLA_VENDOR)
    find_package(LAPACK REQUIRED)
  endif()
else()
  find_package(LAPACK REQUIRED)
endif()

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
  src/diskmalloc.c
  src/xerbla.f
  src/2dolaIO.f
  src/target.f
  src/avker.f
  src/covariance.f
  src/compat.c
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
  $<TARGET_OBJECTS:solapack_common>
)
target_include_directories(set-2drls PRIVATE src)
target_link_libraries(set-2drls PRIVATE linpack ${LAPACK_LIBRARIES})

# --- 2D SOLA inversion executable: 2dsola ---
add_executable(2dsola
  src/2dsola.lanczos.c
  src/2dsola.main.lanczos.f
  src/2dsola.sub.lanczos.f
  src/dlanc_b.f
  src/dkronecker2.f
  src/mgs.f
  $<TARGET_OBJECTS:solapack_common>
)
target_include_directories(2dsola PRIVATE src)
target_link_libraries(2dsola PRIVATE linpack ${LAPACK_LIBRARIES})

# --- Tests ---
enable_testing()
add_test(
  NAME small_test
  COMMAND ${CMAKE_SOURCE_DIR}/test/run_small_test.sh
          ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
)
set_tests_properties(small_test PROPERTIES
  TIMEOUT 120
)

# --- Install rules ---
install(TARGETS set-2drls 2dsola
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY example/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/solapack/example
)
