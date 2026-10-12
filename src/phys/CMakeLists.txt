add_library(cppwrf_phys STATIC
  mp_kessler/kessler.cpp
  mp_wsm6/wsm6.cpp
)

target_include_directories(cppwrf_phys PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_compile_features(cppwrf_phys PUBLIC cxx_std_20)

if(TARGET cppwrf_core)
  target_link_libraries(cppwrf_phys PUBLIC cppwrf_core)
endif()

if(COMMAND cppwrf_device_sources)
  cppwrf_device_sources(
    ${CMAKE_CURRENT_SOURCE_DIR}/mp_kessler/kessler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mp_wsm6/wsm6.cpp
  )
endif()

