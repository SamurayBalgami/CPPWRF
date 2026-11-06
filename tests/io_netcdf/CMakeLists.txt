# SPDX-License-Identifier: GPL-3.0-or-later
find_package(netCDF CONFIG QUIET)
if(NOT TARGET netCDF::netcdf)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(NETCDF netcdf)
  endif()
endif()

find_path(NETCDF_INCDIR netcdf.h HINTS /usr/include /usr/local/include)
find_library(NETCDF_LIB NAMES netcdf HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)

add_executable(test_restart test_restart.cpp)
target_link_libraries(test_restart PRIVATE cppwrf_io_netcdf)
target_include_directories(test_restart PRIVATE ${NETCDF_INCDIR} ${CMAKE_SOURCE_DIR}/src)
target_compile_features(test_restart PRIVATE cxx_std_20)
add_test(NAME test_restart COMMAND test_restart WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_wrf_readers test_wrf_readers.cpp)
target_link_libraries(test_wrf_readers PRIVATE cppwrf_io_netcdf)
target_include_directories(test_wrf_readers PRIVATE ${NETCDF_INCDIR} ${CMAKE_SOURCE_DIR}/src)
target_compile_features(test_wrf_readers PRIVATE cxx_std_20)
add_test(NAME test_wrf_readers COMMAND test_wrf_readers WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(test_wrf_readers PROPERTIES SKIP_RETURN_CODE 77)  # without ${CPPWRF_WRF_RUNS}/wrf_run

# NetCDF-4 output compressed on host threads: the same bytes as HDF5's own deflate filter (src/wrf_io/nc_deflate.hpp)
add_executable(test_nc_deflate test_nc_deflate.cpp)
target_link_libraries(test_nc_deflate PRIVATE cppwrf_wrfio)
target_include_directories(test_nc_deflate PRIVATE ${NETCDF_INCDIR} ${CMAKE_SOURCE_DIR}/src)
target_compile_features(test_nc_deflate PRIVATE cxx_std_20)
add_test(NAME test_nc_deflate COMMAND bash ${CMAKE_SOURCE_DIR}/tests/io_netcdf/nc_deflate_test.sh
         $<TARGET_FILE:test_nc_deflate> ${CMAKE_BINARY_DIR}/test_nc_deflate_work)
set_tests_properties(test_nc_deflate PROPERTIES SKIP_RETURN_CODE 77)
