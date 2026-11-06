# SPDX-License-Identifier: GPL-3.0-or-later
cmake_minimum_required(VERSION 3.24)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(cppwrf_phys_tests LANGUAGES C CXX)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  enable_testing()
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  set(CPPWRF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
else()
  set(CPPWRF_ROOT "${CMAKE_SOURCE_DIR}")
endif()

# Column tests of the physics schemes against WRF v4.8.0's own compiled code. The C++ drivers are built here
# (CPU FP64); each run_test.sh compiles its Fortran driver and links it to the objects of the WRF build that made
# the model's references (tests/phys/wrf_fortran.sh: WRF_DIR = -r8, WRF_R8D8_DIR = -r8 -fdefault-double-8). Without
# that build a test exits 77 (skipped). No WRF source is compiled by a test.

# C++ drivers of Kessler, WSM6 (and its radar reflectivity), the revised MM5 surface layer and YSU (CPU FP64)
foreach(t "mp_kessler;driver_kessler_cpu;driver_kessler.cpp;mp_kessler/kessler.cpp"
          "mp_wsm6;driver_wsm6_cpu;driver_wsm6.cpp;mp_wsm6/wsm6.cpp"
          "mp_radar;driver_radar_cpu;driver_radar.cpp;mp_wsm6/wsm6.cpp"
          "sf_sfclayrev;driver_sfclayrev_cpu;driver_sfclayrev.cpp;sf_sfclayrev/sfclayrev.cpp"
          "bl_ysu;driver_ysu_cpu;driver_ysu.cpp;bl_ysu/ysu.cpp")
  list(GET t 0 dir)
  list(GET t 1 exe)
  list(GET t 2 drv)
  list(GET t 3 src)
  add_executable(${exe} ${dir}/${drv} ${CPPWRF_ROOT}/src/phys/${src} ${CPPWRF_ROOT}/src/core/backend.cpp)
  target_include_directories(${exe} PRIVATE ${CPPWRF_ROOT}/src ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions(${exe} PRIVATE CPPWRF_DOUBLE)
  target_compile_options(${exe} PRIVATE -O3 -march=native -ffp-contract=off)
  target_link_libraries(${exe} PRIVATE OpenMP::OpenMP_CXX)
endforeach()

# Kain-Fritsch, Noah LSM, Thompson (own CMake files: tables and column drivers)
foreach(scheme cu_kf lsm_noah mp_thompson)
  add_subdirectory(${scheme})
endforeach()

# ctest entries (label diff). Schemes with CMake drivers get this build directory; the others (MYJ, Noah ice,
# Noah-MP, MYNN) build their drivers in their own work directory.
set(PHYS_BIN ${CMAKE_CURRENT_BINARY_DIR})
foreach(t "kessler;mp_kessler;${PHYS_BIN}"
          "wsm6;mp_wsm6;${PHYS_BIN}"
          "radar;mp_radar;${PHYS_BIN}"
          "sfclayrev;sf_sfclayrev;${PHYS_BIN}"
          "ysu;bl_ysu;${PHYS_BIN}"
          "bl_myj;bl_myj;${PHYS_BIN}/bl_myj_test"
          "lsm_noah_ice;lsm_noah_ice;${PHYS_BIN}/lsm_noah_ice_test"
          "lsm_noahmp;lsm_noahmp;${PHYS_BIN}/lsm_noahmp_test"
          "bl_mynn;bl_mynn;${PHYS_BIN}/bl_mynn_test"
          "sf_mynn;sf_mynn;${PHYS_BIN}/sf_mynn_test")
  list(GET t 0 name)
  list(GET t 1 dir)
  list(GET t 2 work)
  add_test(NAME test_phys_${name} COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/run_test.sh ${work}
           WORKING_DIRECTORY ${PHYS_BIN})
  set_tests_properties(test_phys_${name} PROPERTIES LABELS diff SKIP_RETURN_CODE 77)
endforeach()
