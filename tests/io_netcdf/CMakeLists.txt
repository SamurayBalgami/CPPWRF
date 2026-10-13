find_package(netCDF CONFIG QUIET)
if(NOT TARGET netCDF::netcdf)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(NETCDF netcdf)
  endif()
endif()

find_path(NETCDF_INCDIR netcdf.h HINTS /usr/include /usr/local/include)
find_library(NETCDF_LIB NAMES netcdf HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)

add_executable(test_wrfout_writer test_wrfout_writer.cpp)
target_link_libraries(test_wrfout_writer PRIVATE cppwrf_io_netcdf)
target_include_directories(test_wrfout_writer PRIVATE ${NETCDF_INCDIR} ${CMAKE_SOURCE_DIR}/src)
target_compile_features(test_wrfout_writer PRIVATE cxx_std_20)
add_test(NAME test_wrfout_writer COMMAND test_wrfout_writer WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

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

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  add_test(NAME test_wrfout_python
    COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test_wrfout_python.py"
            "${CMAKE_SOURCE_DIR}/tests/io_netcdf/test_out/wrfout_d01_test.nc"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  set_tests_properties(test_wrfout_python PROPERTIES DEPENDS test_wrfout_writer)

  # real R0 (10 steps) written as wrfout, then compared field by field with WRF's wrfout
  # (skips inside the script if ~/wrf_real/run_R0_r8 is absent)
  set(_r0dir "${CMAKE_CURRENT_BINARY_DIR}/real_r0_wrfout")
  add_test(NAME test_real_r0_wrfout_run
    COMMAND $<TARGET_FILE:cppwrf> "${CMAKE_SOURCE_DIR}/configs/real_r0_ref.ini" --set "output.dir=${_r0dir}"
            --set output.format=both --set run.run_time=600
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties(test_real_r0_wrfout_run PROPERTIES FIXTURES_SETUP real_r0_wrfout)
  add_test(NAME test_compare_wrfout_real
    COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/compare_wrfout_real.py"
            "${_r0dir}/wrfout_d01_2026-10-12_00:00:00"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  set_tests_properties(test_compare_wrfout_real PROPERTIES FIXTURES_REQUIRED real_r0_wrfout)
endif()
