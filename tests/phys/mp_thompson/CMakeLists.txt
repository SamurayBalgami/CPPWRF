cmake_minimum_required(VERSION 3.24)
project(test_phys_mp_thompson LANGUAGES C CXX Fortran)

# CUDA drivers: standalone builds with a CUDA compiler, or CPPWRF builds with CPPWRF_ENABLE_CUDA
set(THOMPSON_TEST_CUDA OFF)
if(NOT DEFINED CPPWRF_ENABLE_CUDA OR CPPWRF_ENABLE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(THOMPSON_TEST_CUDA ON)
  endif()
endif()

find_package(OpenMP REQUIRED)

set(CPPWRF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

# Fortran flags for bit-for-bit WRF verification
# WRF's own r8d8 configure flags (FCOPTIM): gfortran vectorises through libmvec as in wrf.exe
set(PHYS_FORTRAN_FLAGS "-cpp -ffree-form -fdefault-real-8 -fdefault-double-8 -O2 -ftree-vectorize -funroll-loops -fopenmp -DTHOMPSON_OMP")

# 1. Fortran Driver
add_executable(fortran_driver_thompson
  fortran/driver_thompson.f90
  fortran/module_mp_thompson.F
  fortran/module_mp_radar.F
  fortran/stubs.f90
)
set_target_properties(fortran_driver_thompson PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules"
  COMPILE_FLAGS "${PHYS_FORTRAN_FLAGS}"
)
target_link_libraries(fortran_driver_thompson PRIVATE OpenMP::OpenMP_Fortran)

# 2. Table Test: WRF's thompson_init tables (dump_tables_thompson) against thompson_tables.cpp
add_executable(dump_tables_thompson
  dump_tables_f.f90
  fortran/module_mp_thompson.F
  fortran/module_mp_radar.F
  fortran/stubs.f90
)
set_target_properties(dump_tables_thompson PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules_dump"
  COMPILE_FLAGS "${PHYS_FORTRAN_FLAGS}"
)
target_link_libraries(dump_tables_thompson PRIVATE OpenMP::OpenMP_Fortran)

add_executable(test_tables
  test_tables.cpp
  ${CPPWRF_ROOT}/src/phys/mp_thompson/thompson_tables.cpp
)
target_include_directories(test_tables PRIVATE
  ${CPPWRF_ROOT}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(test_tables PRIVATE CPPWRF_DOUBLE)
target_compile_options(test_tables PRIVATE -O3 -march=native -fopenmp -ffp-contract=off -fno-tree-vectorize)
target_link_libraries(test_tables PRIVATE OpenMP::OpenMP_CXX)

# 3. C++ Driver (CPU FP64)
add_executable(driver_thompson_cpu
  driver_thompson.cpp
  ${CPPWRF_ROOT}/src/phys/mp_thompson/thompson.cpp
  ${CPPWRF_ROOT}/src/phys/mp_thompson/thompson_tables.cpp
  ${CPPWRF_ROOT}/src/core/backend.cpp
)
target_include_directories(driver_thompson_cpu PRIVATE
  ${CPPWRF_ROOT}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(driver_thompson_cpu PRIVATE CPPWRF_DOUBLE)
target_compile_options(driver_thompson_cpu PRIVATE -O3 -march=native -fopenmp -ffp-contract=off -fno-tree-vectorize)
target_link_libraries(driver_thompson_cpu PRIVATE
  OpenMP::OpenMP_CXX
)

# 4. CUDA Drivers (if CUDA enabled)
if(THOMPSON_TEST_CUDA)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/driver_thompson_cuda.cu "#include \"${CMAKE_CURRENT_SOURCE_DIR}/driver_thompson.cpp\"\n")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/thompson_cuda.cu "#include \"${CPPWRF_ROOT}/src/phys/mp_thompson/thompson.cpp\"\n")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/thompson_tables_cuda.cu "#include \"${CPPWRF_ROOT}/src/phys/mp_thompson/thompson_tables.cpp\"\n")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/backend_cuda.cu "#include \"${CPPWRF_ROOT}/src/core/backend.cpp\"\n")

  add_executable(driver_thompson_cuda_dp
    ${CMAKE_CURRENT_BINARY_DIR}/driver_thompson_cuda.cu
    ${CMAKE_CURRENT_BINARY_DIR}/thompson_cuda.cu
    ${CMAKE_CURRENT_BINARY_DIR}/thompson_tables_cuda.cu
    ${CMAKE_CURRENT_BINARY_DIR}/backend_cuda.cu
  )
  target_include_directories(driver_thompson_cuda_dp PRIVATE
    ${CPPWRF_ROOT}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_compile_definitions(driver_thompson_cuda_dp PRIVATE CPPWRF_USE_CUDA CPPWRF_DOUBLE)
  target_compile_options(driver_thompson_cuda_dp PRIVATE --expt-relaxed-constexpr --fmad=false)

  add_executable(driver_thompson_cuda_sp
    ${CMAKE_CURRENT_BINARY_DIR}/driver_thompson_cuda.cu
    ${CMAKE_CURRENT_BINARY_DIR}/thompson_cuda.cu
    ${CMAKE_CURRENT_BINARY_DIR}/thompson_tables_cuda.cu
    ${CMAKE_CURRENT_BINARY_DIR}/backend_cuda.cu
  )
  target_include_directories(driver_thompson_cuda_sp PRIVATE
    ${CPPWRF_ROOT}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_compile_definitions(driver_thompson_cuda_sp PRIVATE CPPWRF_USE_CUDA)
  target_compile_options(driver_thompson_cuda_sp PRIVATE --expt-relaxed-constexpr)
endif()

# 5. CTest entry
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(
  NAME test_phys_mp_thompson
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_phys_mp_thompson PROPERTIES LABELS diff)
