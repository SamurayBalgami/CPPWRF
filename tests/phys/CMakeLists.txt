cmake_minimum_required(VERSION 3.24)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(cppwrf_phys_tests LANGUAGES CXX Fortran)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  enable_testing()
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  set(CPPWRF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
else()
  enable_language(Fortran)
  set(CPPWRF_ROOT "${CMAKE_SOURCE_DIR}")
endif()

# Fortran flags for bit-for-bit WRF verification
set(PHYS_FORTRAN_FLAGS "-ffree-form -fdefault-real-8 -O3 -march=native -ffp-contract=off -fno-tree-vectorize")

set(FORTRAN_KESSLER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/fortran/module_mp_kessler.F")

# 1. Fortran Driver
add_executable(driver_kessler_fortran
  mp_kessler/driver_kessler.f90
  ${FORTRAN_KESSLER_SRC}
)
set_target_properties(driver_kessler_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules"
  COMPILE_FLAGS "${PHYS_FORTRAN_FLAGS}"
)

# 2. C++ Driver (CPU FP64)
add_executable(driver_kessler_cpu
  mp_kessler/driver_kessler.cpp
  ${CPPWRF_ROOT}/src/phys/mp_kessler/kessler.cpp
  ${CPPWRF_ROOT}/src/core/backend.cpp
)
target_include_directories(driver_kessler_cpu PRIVATE
  ${CPPWRF_ROOT}/src
)
target_compile_definitions(driver_kessler_cpu PRIVATE CPPWRF_DOUBLE)
target_compile_options(driver_kessler_cpu PRIVATE -O3 -march=native -ffp-contract=off)
target_link_libraries(driver_kessler_cpu PRIVATE
  OpenMP::OpenMP_CXX
)

# 3. CTest entry
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(
  NAME test_phys_kessler
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mp_kessler/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# ==============================================================================
# WSM6 Scheme
# ==============================================================================
set(FORTRAN_WSM6_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/ccpp_kind_types.F"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/module_libmassv.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/mp_radar.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/mp_wsm6_effectRad.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/mp_wsm6.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/module_mp_wsm6.F"
)

add_executable(driver_wsm6_fortran
  mp_wsm6/driver_wsm6.f90
  ${FORTRAN_WSM6_SRCS}
)
set_target_properties(driver_wsm6_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules_wsm6"
  COMPILE_FLAGS "${PHYS_FORTRAN_FLAGS}"
)

add_executable(driver_wsm6_cpu
  mp_wsm6/driver_wsm6.cpp
  ${CPPWRF_ROOT}/src/phys/mp_wsm6/wsm6.cpp
  ${CPPWRF_ROOT}/src/core/backend.cpp
)
target_include_directories(driver_wsm6_cpu PRIVATE
  ${CPPWRF_ROOT}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(driver_wsm6_cpu PRIVATE CPPWRF_DOUBLE)
target_compile_options(driver_wsm6_cpu PRIVATE -O3 -march=native -ffp-contract=off)
target_link_libraries(driver_wsm6_cpu PRIVATE
  OpenMP::OpenMP_CXX
)

add_test(
  NAME test_phys_wsm6
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mp_wsm6/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

