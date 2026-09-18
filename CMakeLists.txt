cmake_minimum_required(VERSION 3.20)

project(
  sqcalc
  VERSION 0.1.0
  DESCRIPTION "Total structure factors from LAMMPS dump trajectories"
  LANGUAGES Fortran C CXX
)

option(SQC_BUILD_TESTS "Build the sqcalc tests" ON)
option(SQC_USE_SYSTEM_FINUFFT "Use an installed FINUFFT instead of external/finufft" OFF)
option(SQC_ENABLE_CUDA "Build the cuFFT based GPU backend (cufinufft)" OFF)
set(SQC_FINUFFT_ROOT "" CACHE PATH "Prefix of an installed FINUFFT (SQC_USE_SYSTEM_FINUFFT=ON)")
set(SQC_CUDA_ARCHITECTURES "" CACHE STRING "CUDA compute capability to build for (default: detected)")

# FINUFFT's own build system obtains its helper packages (CPM.cmake, findFFTW,
# xsimd, poet) with CPM.  Keep that download cache next to the sources instead
# of inside the build tree, so that reconfiguring, or building a second build
# directory, does not need network access again.
if(NOT DEFINED CPM_SOURCE_CACHE)
  set(CPM_SOURCE_CACHE "${CMAKE_CURRENT_SOURCE_DIR}/.cpm-cache" CACHE PATH "CPM download cache")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")

find_package(OpenMP REQUIRED COMPONENTS Fortran CXX)

# ---------------------------------------------------------------------------
# FFT backend for FINUFFT.
#
# FINUFFT links against FFTW unless it is built with DUCC0.  Point CMake at a
# local FFTW installation here (for example a conda/pixi environment) so that
# FINUFFT's find_package(FFTW) succeeds without downloading anything.
# ---------------------------------------------------------------------------
if(NOT FFTW_ROOT AND NOT DEFINED ENV{FFTW_ROOT})
  file(GLOB _sqc_pixi_envs "$ENV{HOME}/.pixi/envs/*")
  foreach(_sqc_candidate IN ITEMS "$ENV{CONDA_PREFIX}" ${_sqc_pixi_envs} "/usr/local" "/usr")
    if(_sqc_candidate AND EXISTS "${_sqc_candidate}/include/fftw3.h")
      list(APPEND CMAKE_PREFIX_PATH "${_sqc_candidate}")
      set(FFTW_ROOT "${_sqc_candidate}" CACHE PATH "FFTW prefix used by FINUFFT" FORCE)
      message(STATUS "sqcalc: FFTW found in ${_sqc_candidate}")
      break()
    endif()
  endforeach()
endif()

# ---------------------------------------------------------------------------
# FINUFFT: vendored under external/ (git subtree) or an external installation.
# ---------------------------------------------------------------------------
if(SQC_USE_SYSTEM_FINUFFT)
  if(SQC_FINUFFT_ROOT)
    list(APPEND CMAKE_PREFIX_PATH "${SQC_FINUFFT_ROOT}")
  endif()
  find_package(finufft REQUIRED CONFIG)
  set(_sqc_finufft_target finufft::finufft)
else()
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/finufft/CMakeLists.txt")
    message(
      FATAL_ERROR
        "external/finufft is missing. Add it with\n"
        "  git subtree add --prefix=external/finufft <finufft-url> <tag> --squash\n"
        "or configure with -DSQC_USE_SYSTEM_FINUFFT=ON."
    )
  endif()
  set(FINUFFT_USE_CPU ON CACHE BOOL "" FORCE)
  set(FINUFFT_USE_CUDA ${SQC_ENABLE_CUDA} CACHE BOOL "" FORCE)
  set(FINUFFT_USE_OPENMP ON CACHE BOOL "" FORCE)
  if(SQC_ENABLE_CUDA)
    # cufinufft is linked as a shared library: its device code is then linked
    # inside libcufinufft and a plain Fortran link against it is enough.
    set(FINUFFT_STATIC_LINKING OFF CACHE BOOL "" FORCE)
    set(FINUFFT_SHARED_LINKING ON CACHE BOOL "" FORCE)
  else()
    set(FINUFFT_STATIC_LINKING ON CACHE BOOL "" FORCE)
  endif()
  set(FINUFFT_BUILD_FORTRAN OFF CACHE BOOL "" FORCE)
  set(FINUFFT_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(FINUFFT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(FINUFFT_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
  set(FINUFFT_BUILD_MATLAB OFF CACHE BOOL "" FORCE)
  set(FINUFFT_BUILD_DOCS OFF CACHE BOOL "" FORCE)
  set(FINUFFT_BUILD_DEVEL OFF CACHE BOOL "" FORCE)
  set(FINUFFT_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  add_subdirectory(external/finufft)
  set(_sqc_finufft_target finufft)
endif()

if(SQC_ENABLE_CUDA)
  # FINUFFT defaults to "native", which needs a visible GPU while configuring;
  # fall back to an explicit compute capability so that CPU-only build hosts
  # (or containers without device access) can still cross compile.
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    if(SQC_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES "${SQC_CUDA_ARCHITECTURES}" CACHE STRING "CUDA SMs" FORCE)
    else()
      find_program(_sqc_nvidia_smi nvidia-smi)
      if(_sqc_nvidia_smi)
        execute_process(
          COMMAND ${_sqc_nvidia_smi} --query-gpu=compute_cap --format=csv,noheader
          OUTPUT_VARIABLE _sqc_compute_caps OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
        )
      endif()
      if(_sqc_compute_caps MATCHES "^[0-9]+\\.[0-9]+")
        string(REPLACE "." "" _sqc_arch "${_sqc_compute_caps}")
        set(CMAKE_CUDA_ARCHITECTURES "${_sqc_arch}" CACHE STRING "CUDA SMs" FORCE)
      else()
        message(FATAL_ERROR "Cannot detect a CUDA GPU; pass -DSQC_CUDA_ARCHITECTURES=<cc>, e.g. 75.")
      endif()
    endif()
  endif()
  find_package(CUDAToolkit REQUIRED)
  message(STATUS "sqcalc: CUDA backend enabled (arch ${CMAKE_CUDA_ARCHITECTURES})")
endif()

# ---------------------------------------------------------------------------
# sqcalc
# ---------------------------------------------------------------------------
add_library(
  sqcalc_lib
  STATIC
  src/sqc_kinds.f90
  src/sqc_element_data.f90
  src/sqc_elements.f90
  src/sqc_cell.f90
  src/sqc_dump.f90
  src/sqc_weights.f90
  src/sqc_finufft.f90
  src/sqc_structure.f90
  src/sqc_options.f90
)
if(SQC_ENABLE_CUDA)
  target_sources(sqcalc_lib PRIVATE src/sqc_cufinufft.f90 src/sqc_gpu.f90)
  target_link_libraries(sqcalc_lib PUBLIC cufinufft CUDA::cudart)
  target_compile_definitions(sqcalc_lib PUBLIC SQC_ENABLE_CUDA)
endif()
target_link_libraries(sqcalc_lib PUBLIC OpenMP::OpenMP_Fortran ${_sqc_finufft_target})
target_include_directories(sqcalc_lib PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(sqcalc_lib PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")

target_compile_options(
  sqcalc_lib
  # -Wno-unused-dummy-argument: the deferred type bound procedures share one
  # signature, so a method may legitimately ignore e.g. the weight scheme.
  # -cpp: the driver selects the GPU backend at compile time.
  PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Wall;-Wextra;-Wno-unused-dummy-argument;-cpp>
          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>
)

add_executable(sqcalc app/main.f90)
target_link_libraries(sqcalc PRIVATE sqcalc_lib)
target_compile_options(sqcalc PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)

include(GNUInstallDirs)
install(TARGETS sqcalc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(SQC_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
