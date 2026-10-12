find_package(PkgConfig QUIET)
find_path(NETCDF_INCLUDE_DIR netcdf.h HINTS /usr/include /usr/local/include)
find_library(NETCDF_C_LIBRARY NAMES netcdf HINTS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)

if(NETCDF_INCLUDE_DIR AND NETCDF_C_LIBRARY)
  message(STATUS "Found NetCDF: ${NETCDF_C_LIBRARY} (include: ${NETCDF_INCLUDE_DIR})")
  
  add_library(cppwrf_init_real STATIC
    vcoord.cpp
    met_em_reader.cpp
    base_state.cpp
    vert_interp.cpp
    hydrostatic.cpp
    bdy.cpp
    init_real.cpp
  )
  target_include_directories(cppwrf_init_real PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${NETCDF_INCLUDE_DIR}
  )
  target_link_libraries(cppwrf_init_real PUBLIC
    ${NETCDF_C_LIBRARY}
  )
  target_compile_features(cppwrf_init_real PUBLIC cxx_std_20)
else()
  message(STATUS "NetCDF C library not found; skipping cppwrf_init_real")
endif()
