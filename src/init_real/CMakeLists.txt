# C++ port of WRF's real.exe (met_em -> wrfinput_d01, wrfbdy_d01); host code, both precisions.
add_library(cppwrf_init_real STATIC namelist.cpp met_em.cpp nc_out.cpp real_domain.cpp wrf_schema.cpp)
target_include_directories(cppwrf_init_real PUBLIC ${CMAKE_SOURCE_DIR}/src ${NETCDF_INCLUDE_DIR})
target_link_libraries(cppwrf_init_real PUBLIC ${NETCDF_LIBRARY} mvec m)
target_compile_features(cppwrf_init_real PUBLIC cxx_std_20)
target_compile_options(cppwrf_init_real PRIVATE -ffp-contract=off)

add_executable(cppwrf_real main_real.cpp)
target_link_libraries(cppwrf_real PRIVATE cppwrf_init_real)
target_compile_options(cppwrf_real PRIVATE -ffp-contract=off)
target_compile_definitions(cppwrf_real PRIVATE CPPWRF_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
