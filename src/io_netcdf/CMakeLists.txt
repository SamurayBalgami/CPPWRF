find_package(netCDF CONFIG QUIET)
if(NOT TARGET netCDF::netcdf)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(NETCDF netcdf)
  endif()
endif()

find_path(NETCDF_INCDIR netcdf.h HINTS /usr/include /usr/local/include)
find_library(NETCDF_LIB NAMES netcdf HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)

add_library(cppwrf_io_netcdf STATIC
  wrfout_writer.cpp
  restart.cpp
)

target_include_directories(cppwrf_io_netcdf PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${NETCDF_INCDIR}
)

target_compile_features(cppwrf_io_netcdf PUBLIC cxx_std_20)

if(TARGET netCDF::netcdf)
  target_link_libraries(cppwrf_io_netcdf PUBLIC netCDF::netcdf)
elseif(NETCDF_LIB)
  target_link_libraries(cppwrf_io_netcdf PUBLIC ${NETCDF_LIB})
else()
  target_link_libraries(cppwrf_io_netcdf PUBLIC ${NETCDF_LIBRARIES})
endif()
