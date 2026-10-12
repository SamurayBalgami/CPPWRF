cmake_minimum_required(VERSION 3.24)

set(CPPWRF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

# Fortran flags for bit-for-bit WRF verification
set(PHYS_FORTRAN_FLAGS "-cpp -ffree-form -fdefault-real-8 -fdefault-double-8 -O3 -march=native -ffp-contract=off -fno-tree-vectorize")

# 1. Fortran Driver
add_executable(driver_noah_fortran
  fortran/driver_noah.f90
  fortran/module_sf_noahlsm.F
  fortran/stubs.f90
)
set_target_properties(driver_noah_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules"
  COMPILE_FLAGS "${PHYS_FORTRAN_FLAGS}"
)

# 2. C++ Driver (CPU FP64)
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
if(CMAKE_CUDA_COMPILER)
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
