add_library(cppwrf_phys_noah STATIC
  noah.cpp
  noah_tables.cpp
)

target_include_directories(cppwrf_phys_noah PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_compile_features(cppwrf_phys_noah PUBLIC cxx_std_20)

if(TARGET cppwrf_core)
  target_link_libraries(cppwrf_phys_noah PUBLIC cppwrf_core)
endif()

if(COMMAND cppwrf_device_sources)
  cppwrf_device_sources(
    ${CMAKE_CURRENT_SOURCE_DIR}/noah.cpp
  )
endif()
