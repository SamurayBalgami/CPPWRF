# WRF I/O layer shared by cppwrf_real and the model: namelist.input reader, the WRF I/O schema
# (tools/wrf_schema), the NetCDF writer with WRF's conventions and output_wrf's global attributes.
add_library(cppwrf_wrfio STATIC namelist.cpp nc_out.cpp wrf_schema.cpp global_atts.cpp)
target_include_directories(cppwrf_wrfio PUBLIC ${CMAKE_SOURCE_DIR}/src ${NETCDF_INCLUDE_DIR})
target_link_libraries(cppwrf_wrfio PUBLIC ${NETCDF_LIBRARY})
target_compile_features(cppwrf_wrfio PUBLIC cxx_std_20)
