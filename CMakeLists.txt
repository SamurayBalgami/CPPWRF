cmake_minimum_required(VERSION 3.24)
project(CPPWRF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(CPPWRF_DOUBLE "Use double precision (real = double)" OFF)
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  set(_cuda_default ON)
else()
  set(_cuda_default OFF)
endif()
option(CPPWRF_ENABLE_CUDA "Build the CUDA backend" ${_cuda_default})

if(CPPWRF_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  find_package(CUDAToolkit REQUIRED)
else()
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Strict IEEE: no fast-math anywhere (determinism and honest stability tests).
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3")

enable_testing()
add_subdirectory(src/io)

file(GLOB CPPWRF_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/src/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/dyn/*.cpp
  ${CMAKE_SOURCE_DIR}/src/cases/*.cpp)

add_library(cppwrf_core STATIC ${CPPWRF_SOURCES})
target_include_directories(cppwrf_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(cppwrf_core PUBLIC cppwrf_io)
if(EXISTS ${CMAKE_SOURCE_DIR}/src/dyn/small_step.cpp)
  target_compile_definitions(cppwrf_core PUBLIC CPPWRF_HAVE_SMALL_STEP)
endif()
if(CPPWRF_DOUBLE)
  target_compile_definitions(cppwrf_core PUBLIC CPPWRF_DOUBLE)
endif()

# Single-source kernels: every source of the core is compiled by nvcc in CUDA builds.
function(cppwrf_device_sources)
  if(CPPWRF_ENABLE_CUDA)
    set_source_files_properties(${ARGN} PROPERTIES LANGUAGE CUDA)
  endif()
endfunction()

if(CPPWRF_ENABLE_CUDA)
  cppwrf_device_sources(${CPPWRF_SOURCES})
  target_compile_definitions(cppwrf_core PUBLIC CPPWRF_USE_CUDA)
  target_compile_options(cppwrf_core PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --extended-lambda>)
  target_link_libraries(cppwrf_core PUBLIC CUDA::cudart)
else()
  target_link_libraries(cppwrf_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/src/main.cpp)
  add_executable(cppwrf src/main.cpp)
  cppwrf_device_sources(src/main.cpp)
  target_link_libraries(cppwrf PRIVATE cppwrf_core)
endif()

# Unit tests (tests/unit/*.cpp; tests/unit/io belongs to cppwrf_io)
file(GLOB CPPWRF_UNIT_TESTS CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/unit/*.cpp)
foreach(src ${CPPWRF_UNIT_TESTS})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  cppwrf_device_sources(${src})
  target_link_libraries(${name} PRIVATE cppwrf_core)
  add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
