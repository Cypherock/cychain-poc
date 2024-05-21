add_custom_target (benchs_build COMMAND ${CMAKE_COMMAND} -E sleep 0)
add_custom_target (benchs COMMAND ${CMAKE_COMMAND} -E sleep 0)
add_dependencies (benchs benchs_build)

### List of benchmarks
set (BENCHS_LIST benchmark_CL_HSMqk
                 benchmark_CL_HSM2k
                 benchmark_Paillier
                 benchmark_JoyeLibert
                 benchmark_threshold_ECDSA
                 )

### Automatically create target and test from TESTS_LIST
foreach (bench ${BENCHS_LIST})
  add_executable (${bench}_build ${bench}.cpp)
  set_target_properties (${bench}_build PROPERTIES OUTPUT_NAME ${bench})
  target_link_libraries (${bench}_build bicycl)
  add_custom_target (${bench} COMMAND ${bench}_build COMMENT "Running ${bench}")
  add_dependencies (benchs_build ${bench}_build)
  add_dependencies (benchs ${bench})
endforeach ()
