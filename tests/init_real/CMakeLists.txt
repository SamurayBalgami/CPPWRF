if(TARGET cppwrf_init_real)
  add_executable(compare_wrfinput compare_wrfinput.cpp)
  target_link_libraries(compare_wrfinput PRIVATE cppwrf_init_real)
  target_compile_features(compare_wrfinput PRIVATE cxx_std_20)

  add_test(NAME test_compare_wrfinput COMMAND compare_wrfinput)
endif()
