# tools/wrf2cppwrf.py (WRF namelist -> CPPWRF configuration): ctest -L unit; the 2-step runs use this build's cppwrf
add_test(NAME wrf2cppwrf_translator COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(wrf2cppwrf_translator PROPERTIES LABELS unit SKIP_RETURN_CODE 77)
if(TARGET cppwrf)
  set_tests_properties(wrf2cppwrf_translator PROPERTIES ENVIRONMENT "CPPWRF_BIN=$<TARGET_FILE:cppwrf>")
endif()
