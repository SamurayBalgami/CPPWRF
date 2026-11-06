# SPDX-License-Identifier: GPL-3.0-or-later
cmake_minimum_required(VERSION 3.24)

set(CPPWRF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

# The Fortran side (WRF's module_sf_noahlsm.o, module_sf_noahdrv.o) is compiled and linked by run_test.sh
# (tests/phys/wrf_fortran.sh)

# 1. C++ Driver (CPU FP64)
add_executable(driver_noah_cpu
  driver_noah.cpp
  ${CPPWRF_ROOT}/src/phys/lsm_noah/noah.cpp
  ${CPPWRF_ROOT}/src/phys/lsm_noah/noah_tables.cpp
  ${CPPWRF_ROOT}/src/core/backend.cpp
)
target_include_directories(driver_noah_cpu PRIVATE
  ${CPPWRF_ROOT}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(driver_noah_cpu PRIVATE CPPWRF_DOUBLE)
target_compile_options(driver_noah_cpu PRIVATE -O3 -march=native -ffp-contract=off)
target_link_libraries(driver_noah_cpu PRIVATE
  OpenMP::OpenMP_CXX
)

# 3. CUDA Drivers (if CUDA enabled)
if(CMAKE_CUDA_COMPILER AND (NOT DEFINED CPPWRF_ENABLE_CUDA OR CPPWRF_ENABLE_CUDA))
  add_executable(driver_noah_cuda_dp
    driver_noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah_tables.cpp
    ${CPPWRF_ROOT}/src/core/backend.cpp
  )
  set_source_files_properties(
    driver_noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah_tables.cpp
    ${CPPWRF_ROOT}/src/core/backend.cpp
    PROPERTIES LANGUAGE CUDA
  )
  target_include_directories(driver_noah_cuda_dp PRIVATE
    ${CPPWRF_ROOT}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_compile_definitions(driver_noah_cuda_dp PRIVATE CPPWRF_USE_CUDA CPPWRF_DOUBLE)
  target_compile_options(driver_noah_cuda_dp PRIVATE --expt-relaxed-constexpr)

  add_executable(driver_noah_cuda_sp
    driver_noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah_tables.cpp
    ${CPPWRF_ROOT}/src/core/backend.cpp
  )
  set_source_files_properties(
    driver_noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah.cpp
    ${CPPWRF_ROOT}/src/phys/lsm_noah/noah_tables.cpp
    ${CPPWRF_ROOT}/src/core/backend.cpp
    PROPERTIES LANGUAGE CUDA
  )
  target_include_directories(driver_noah_cuda_sp PRIVATE
    ${CPPWRF_ROOT}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_compile_definitions(driver_noah_cuda_sp PRIVATE CPPWRF_USE_CUDA)
  target_compile_options(driver_noah_cuda_sp PRIVATE --expt-relaxed-constexpr)
endif()

# 4. CTest entry
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(
  NAME test_phys_noah
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_phys_noah PROPERTIES LABELS diff SKIP_RETURN_CODE 77)
