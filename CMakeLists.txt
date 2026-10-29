cmake_minimum_required(VERSION 3.24)
project(CPPWRF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(CPPWRF_DOUBLE "Use double precision (real = double)" OFF)
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  set(_cuda_default ON)
else()
  set(_cuda_default OFF)
endif()
option(CPPWRF_ENABLE_CUDA "Build the CUDA backend" ${_cuda_default})

if(CPPWRF_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)  # the GPU of the build machine; pass e.g. 80;90 for servers
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  find_package(CUDAToolkit REQUIRED)
else()
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Strict IEEE: no fast-math anywhere (determinism and honest stability tests).
# -ffp-contract=off in every build type, after the user's flags: no FMA contraction on the CPU, as WRF's
# gfortran build on x86-64 (no -march: no FMA instructions); with -march=native or skylake g++ would
# otherwise fuse a*b+c (its default is -ffp-contract=fast). CUDA builds: the same for nvcc's host code.
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off> $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-ffp-contract=off>)
option(CPPWRF_NATIVE "Tune for the build machine's CPU (-march=native); OFF for portable binaries" ON)
if(CPPWRF_NATIVE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()
set(CMAKE_CUDA_FLAGS_RELEASE "-O3")
# flags that change results: reassociation, reciprocals, flush-to-zero, approximate division / sqrt
string(TOUPPER "${CMAKE_BUILD_TYPE}" _cppwrf_bt)
if("${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_cppwrf_bt}}" MATCHES
   "(-Ofast|-f(fast-math|unsafe-math-optimizations|associative-math|reciprocal-math|finite-math-only)|-ffp-contract=(fast|on))")
  message(FATAL_ERROR "CMAKE_CXX_FLAGS: '${CMAKE_MATCH_1}' changes the floating-point results (CPPWRF needs strict IEEE)")
endif()
if("${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${_cppwrf_bt}}" MATCHES
   "(use_fast_math|ftz=true|prec-div=false|prec-sqrt=false|-Xcompiler=?[^ ]*(-Ofast|-ffast-math))")
  message(FATAL_ERROR "CMAKE_CUDA_FLAGS: '${CMAKE_MATCH_1}' changes the floating-point results (CPPWRF needs strict IEEE)")
endif()

enable_testing()
add_subdirectory(src/io)

file(GLOB CPPWRF_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/src/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/dyn/*.cpp
  ${CMAKE_SOURCE_DIR}/src/cases/*.cpp)
# Column physics used by the dycore (P1/P2 suite: sfclayrev + Noah + YSU + KF and its driver glue)
list(APPEND CPPWRF_SOURCES
  ${CMAKE_SOURCE_DIR}/src/phys/driver/suite.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/surface_driver.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/pbl_driver.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/cumulus_driver.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/phy_prep.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/tendency_coupling.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/landuse_tables.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/radiation.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw_batch.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw_tables.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw_batch.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw_tables.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_common/cam_gases.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_common/cam_ozone.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/sf_sfclayrev/sfclayrev.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/sf_myjsfc/myjsfc.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/bl_ysu/ysu.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/bl_myj/myjpbl.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/lsm_noah/noah.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/lsm_noah/noah_tables.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/lsm_noahmp/noahmp.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/lsm_noahmp/noahmp_tables.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/cu_kf/kf.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/cu_kf/kf_tables.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/mp_thompson/thompson.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/mp_thompson/thompson_tables.cpp)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/phys/mp_thompson/thompson.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/suite.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/surface_driver.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/pbl_driver.cpp
  ${CMAKE_SOURCE_DIR}/src/phys/driver/cumulus_driver.cpp
  PROPERTIES COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>")

add_library(cppwrf_core STATIC ${CPPWRF_SOURCES})
# tables shipped with the sources (data/rrtmg)
target_compile_definitions(cppwrf_core PRIVATE CPPWRF_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_include_directories(cppwrf_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(cppwrf_core PUBLIC cppwrf_io)
# glibc libmvec: Grell-Freitas-Li uses the 2-lane vector pow where WRF's gfortran build vectorises (cu_gfl)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
  target_link_libraries(cppwrf_core PUBLIC mvec)
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/src/dyn/small_step.cpp)
  target_compile_definitions(cppwrf_core PUBLIC CPPWRF_HAVE_SMALL_STEP)
endif()
if(CPPWRF_DOUBLE)
  target_compile_definitions(cppwrf_core PUBLIC CPPWRF_DOUBLE)
endif()

# NetCDF-C library
find_package(netCDF CONFIG QUIET)
if(NOT TARGET netCDF::netcdf)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(NETCDF netcdf)
  endif()
endif()
find_path(NETCDF_INCDIR netcdf.h HINTS /usr/include /usr/local/include)
if(TARGET netCDF::netcdf OR NETCDF_FOUND OR NETCDF_INCDIR)
  if(EXISTS ${CMAKE_SOURCE_DIR}/src/io_netcdf)
    add_subdirectory(src/io_netcdf)
    target_link_libraries(cppwrf_core PUBLIC cppwrf_io_netcdf)
    target_compile_definitions(cppwrf_core PUBLIC CPPWRF_HAVE_NETCDF)
  endif()
  if(EXISTS ${CMAKE_SOURCE_DIR}/tests/io_netcdf)
    add_subdirectory(tests/io_netcdf)
  endif()
  # C++ real.exe (host only)
  find_library(NETCDF_LIBRARY NAMES netcdf HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)
  set(NETCDF_INCLUDE_DIR ${NETCDF_INCDIR})
  add_subdirectory(src/wrf_io)
  target_link_libraries(cppwrf_core PUBLIC cppwrf_wrfio)
  add_subdirectory(src/init_real)
  # cppwrf_real against WRF real.exe output (skipped, exit 77, without ${CPPWRF_WRF_RUNS}/real)
  add_test(NAME real_vs_wrf COMMAND bash ${CMAKE_SOURCE_DIR}/tests/init_real/run_real_test.sh ${CMAKE_BINARY_DIR}
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties(real_vs_wrf PROPERTIES SKIP_RETURN_CODE 77)
  # wrfout/wrfrst files identical to WRF v4.8.0's (CPU FP64 builds; skipped, exit 77, without ${CPPWRF_WRF_RUNS})
  if(CPPWRF_DOUBLE AND NOT CPPWRF_ENABLE_CUDA)
    add_test(NAME wrf_files_vs_wrf COMMAND bash ${CMAKE_SOURCE_DIR}/tests/verify/run_wrf_files_test.sh ${CMAKE_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    set_tests_properties(wrf_files_vs_wrf PROPERTIES SKIP_RETURN_CODE 77)
  endif()
endif()

# Single-source kernels: every source of the core is compiled by nvcc in CUDA builds.
function(cppwrf_device_sources)
  if(CPPWRF_ENABLE_CUDA)
    set_source_files_properties(${ARGN} PROPERTIES LANGUAGE CUDA)
  endif()
endfunction()

if(CPPWRF_ENABLE_CUDA)
  cppwrf_device_sources(${CPPWRF_SOURCES})
  # RRTMG host code (tables, constants, the column drivers of the checks): C++; the device path is
  # radiation.cpp + rrtmg_lw_batch.cpp + rrtmg_sw_batch.cpp
  set_source_files_properties(
    ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw.cpp ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw_tables.cpp
    ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw.cpp ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw_tables.cpp
    ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_common/cam_gases.cpp ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_common/cam_ozone.cpp
    PROPERTIES LANGUAGE CXX)
  # RRTMG on the device: no FMA contraction (WRF's x86 arithmetic; table indices int(tblint*x + 0.5)
  # are discontinuous, so contracted products move them)
  set_source_files_properties(
    ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw_batch.cpp ${CMAKE_SOURCE_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw_batch.cpp
    PROPERTIES COMPILE_OPTIONS "--fmad=false")
  target_compile_definitions(cppwrf_core PUBLIC CPPWRF_USE_CUDA)
  target_compile_options(cppwrf_core PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --extended-lambda>)
  target_link_libraries(cppwrf_core PUBLIC CUDA::cudart)
  # host OpenMP for the radiation columns (C++ sources only)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(cppwrf_core PUBLIC OpenMP::OpenMP_CXX)
else()
  target_link_libraries(cppwrf_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/src/main.cpp)
  add_executable(cppwrf src/main.cpp)
  cppwrf_device_sources(src/main.cpp)
  target_link_libraries(cppwrf PRIVATE cppwrf_core)
endif()

# Unit tests (tests/unit/*.cpp; tests/unit/io belongs to cppwrf_io)
file(GLOB CPPWRF_UNIT_TESTS CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/unit/*.cpp)
foreach(src ${CPPWRF_UNIT_TESTS})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  cppwrf_device_sources(${src})
  target_link_libraries(${name} PRIVATE cppwrf_core)
  add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  # OpenMP threads that spin while waiting stall the tests for hours on a loaded machine
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;OMP_WAIT_POLICY=passive")
endforeach()

# Namelist translator (tools/wrf2cppwrf.py; python3, no WRF data needed)
add_subdirectory(tests/tools)

# Differential tests against the original WRF Fortran (gfortran; slow): physics schemes, the
# physics driver glue, diffusion kernels, lateral boundaries. ctest -L diff runs them.
option(CPPWRF_DIFF_TESTS "Register the WRF-Fortran differential tests with ctest" OFF)
if(CPPWRF_DIFF_TESTS)
  find_package(OpenMP REQUIRED COMPONENTS CXX)  # the CPU reference drivers (also in CUDA builds)
  add_subdirectory(tests/phys)
  foreach(t dyn/diffusion dyn/halo dyn/lateral_bc dyn/scalars phys/suite nest)  # need WRF_DIR (a WRF v4.8.0 -r8 build); nest also NETCDF
    string(REPLACE "/" "_" tname "diff_${t}")
    add_test(NAME ${tname} COMMAND bash ${CMAKE_SOURCE_DIR}/tests/${t}/run_test.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    set_tests_properties(${tname} PROPERTIES LABELS diff SKIP_RETURN_CODE 77)  # 77: no WRF build / runs
  endforeach()
  # column replays of the model's cumulus calls through WRF's objects (CPU FP64 cppwrf; need
  # CPPWRF_WRF_RUNS and WRF_R8D8_DIR, a WRF v4.8.0 -r8 -double-8 build)
  if(CPPWRF_DOUBLE AND NOT CPPWRF_ENABLE_CUDA)
    foreach(t cu_gfl cu_ntiedtke)
      add_test(NAME diff_phys_${t} COMMAND bash ${CMAKE_SOURCE_DIR}/tests/phys/${t}/run_test.sh ${CMAKE_BINARY_DIR}
               WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      set_tests_properties(diff_phys_${t} PROPERTIES LABELS diff SKIP_RETURN_CODE 77)
    endforeach()
  endif()
endif()
