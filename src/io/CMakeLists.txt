add_library(cppwrf_io STATIC
  config.cpp
  output.cpp
)

target_include_directories(cppwrf_io PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_compile_features(cppwrf_io PUBLIC cxx_std_20)

# Resolve test and script directories
if (EXISTS "${CMAKE_SOURCE_DIR}/tests/unit/io")
  set(CPPWRF_IO_TEST_DIR "${CMAKE_SOURCE_DIR}/tests/unit/io")
else()
  set(CPPWRF_IO_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../tests/unit/io")
endif()

if (EXISTS "${CMAKE_SOURCE_DIR}/scripts")
  set(CPPWRF_SCRIPTS_DIR "${CMAKE_SOURCE_DIR}/scripts")
else()
  set(CPPWRF_SCRIPTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../scripts")
endif()

# Unit tests
add_executable(test_io_config "${CPPWRF_IO_TEST_DIR}/test_config.cpp")
target_link_libraries(test_io_config PRIVATE cppwrf_io)
target_compile_features(test_io_config PRIVATE cxx_std_20)
add_test(NAME test_io_config COMMAND test_io_config)

add_executable(test_io_output "${CPPWRF_IO_TEST_DIR}/test_output.cpp")
target_link_libraries(test_io_output PRIVATE cppwrf_io)
target_compile_features(test_io_output PRIVATE cxx_std_20)
add_test(NAME test_io_output COMMAND test_io_output)

add_executable(test_io_cross_writer "${CPPWRF_IO_TEST_DIR}/test_cross_writer.cpp")
target_link_libraries(test_io_cross_writer PRIVATE cppwrf_io)
target_compile_features(test_io_cross_writer PRIVATE cxx_std_20)

find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
  add_test(NAME test_python_io
    COMMAND ${Python3_EXECUTABLE} "${CPPWRF_SCRIPTS_DIR}/test_cppwrf_io.py" --cross-writer "$<TARGET_FILE:test_io_cross_writer>"
  )
endif()
