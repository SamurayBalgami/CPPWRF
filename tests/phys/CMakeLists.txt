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

# WRF's own -r8 configure flags (FCOPTIM): gfortran vectorises through libmvec as wrf.exe does
set(WRF_FORTRAN_FLAGS "-ffree-form -fdefault-real-8 -DDOUBLE_PRECISION -O2 -ftree-vectorize -funroll-loops")

set(FORTRAN_KESSLER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/fortran/module_mp_kessler.F")

# 1. Fortran Driver
add_executable(driver_kessler_fortran
  mp_kessler/driver_kessler.f90
  ${FORTRAN_KESSLER_SRC}
)
set_target_properties(driver_kessler_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules"
  COMPILE_FLAGS "${WRF_FORTRAN_FLAGS}"
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
set_tests_properties(test_phys_kessler PROPERTIES LABELS diff)

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
  COMPILE_FLAGS "${WRF_FORTRAN_FLAGS}"
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
set_tests_properties(test_phys_wsm6 PROPERTIES LABELS diff)

# ==============================================================================
# Revised MM5 Surface Layer Scheme (sf_sfclayrev)
# ==============================================================================
set(FORTRAN_SFCLAYREV_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/ccpp_kind_types.F"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/sf_sfclayrev.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/module_sf_sfclayrev.F"
)

add_executable(driver_sfclayrev_fortran
  sf_sfclayrev/driver_sfclayrev.f90
  ${FORTRAN_SFCLAYREV_SRCS}
)
# WRF's own -r8 configure flags: gfortran vectorises the pow/exp loops through libmvec as in wrf.exe
set_target_properties(driver_sfclayrev_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules_sfclayrev"
  COMPILE_FLAGS "-ffree-form -fdefault-real-8 -DDOUBLE_PRECISION -O2 -ftree-vectorize -funroll-loops -DEM_CORE=1"
)

add_executable(driver_sfclayrev_cpu
  sf_sfclayrev/driver_sfclayrev.cpp
  ${CPPWRF_ROOT}/src/phys/sf_sfclayrev/sfclayrev.cpp
  ${CPPWRF_ROOT}/src/core/backend.cpp
)
target_include_directories(driver_sfclayrev_cpu PRIVATE
  ${CPPWRF_ROOT}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(driver_sfclayrev_cpu PRIVATE CPPWRF_DOUBLE)
target_compile_options(driver_sfclayrev_cpu PRIVATE -O3 -march=native -ffp-contract=off)
target_link_libraries(driver_sfclayrev_cpu PRIVATE
  OpenMP::OpenMP_CXX
)

add_test(
  NAME test_phys_sfclayrev
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/sf_sfclayrev/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_phys_sfclayrev PROPERTIES LABELS diff)

# ==============================================================================
# YSU PBL Scheme (bl_ysu)
# ==============================================================================
set(FORTRAN_YSU_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/ccpp_kind_types.F"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/bl_ysu.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/fortran/module_bl_ysu.F"
)

add_executable(driver_ysu_fortran
  bl_ysu/driver_ysu.f90
  ${FORTRAN_YSU_SRCS}
)
set_target_properties(driver_ysu_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules_ysu"
  COMPILE_FLAGS "${WRF_FORTRAN_FLAGS} -ffree-line-length-none"
)

add_executable(driver_ysu_cpu
  bl_ysu/driver_ysu.cpp
  ${CPPWRF_ROOT}/src/phys/bl_ysu/ysu.cpp
  ${CPPWRF_ROOT}/src/core/backend.cpp
)
target_include_directories(driver_ysu_cpu PRIVATE
  ${CPPWRF_ROOT}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(driver_ysu_cpu PRIVATE CPPWRF_DOUBLE)
target_compile_options(driver_ysu_cpu PRIVATE -O3 -march=native -ffp-contract=off)
target_link_libraries(driver_ysu_cpu PRIVATE
  OpenMP::OpenMP_CXX
)

add_test(
  NAME test_phys_ysu
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bl_ysu/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_phys_ysu PROPERTIES LABELS diff)

# ==============================================================================
# MYJ PBL + Eta Surface Layer Scheme (bl_myj + sf_myjsfc)
# ==============================================================================
add_test(
  NAME test_phys_bl_myj
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/bl_myj/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}/bl_myj_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_phys_bl_myj PROPERTIES LABELS diff)

# ==============================================================================
# Kain-Fritsch, Noah LSM, Thompson (own CMake projects: tables and column drivers)
# ==============================================================================
foreach(scheme cu_kf lsm_noah mp_thompson)
  add_subdirectory(${scheme})
endforeach()

# ==============================================================================
# Noah-MP v5.1.2 Land Surface Model (lsm_noahmp, sf_surface_physics=4)
# ==============================================================================
add_test(
  NAME test_phys_lsm_noahmp
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/lsm_noahmp/run_test.sh ${CMAKE_CURRENT_BINARY_DIR}/lsm_noahmp_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_phys_lsm_noahmp PROPERTIES LABELS diff)
