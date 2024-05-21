include (CTest)
if ("${CMAKE_GENERATOR}" STREQUAL "Unix Makefiles")
  set (_cmd ${CMAKE_CTEST_COMMAND} --output-on-failure $(ARGS))
else ()
  set (_cmd ${CMAKE_CTEST_COMMAND} --output-on-failure)
endif ()
add_custom_target (check COMMAND ${_cmd} COMMENT "Running the tests")

### List of tests
set (TESTS_LIST test_gmp_extras
                test_qfi
                test_CL_HSMqk
                test_CL_HSM2k
                test_Paillier
                test_Joye_Libert
                test_ec
                test_threshold_ECDSA
                )

### Automatically create target and test from TESTS_LIST
foreach (test ${TESTS_LIST})
  # /!\ The test and the target building the executable have the same name
  add_executable (${test} ${test}.cpp)
  target_link_libraries (${test} PUBLIC bicycl)
  add_dependencies (check ${test})
  add_test (NAME ${test} COMMAND ${test}
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach ()
