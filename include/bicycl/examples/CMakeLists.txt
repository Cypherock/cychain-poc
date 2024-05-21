add_custom_target (examples COMMAND ${CMAKE_COMMAND} -E sleep 0)

add_executable (chia-vdf chia-vdf.cpp)
target_link_libraries (chia-vdf PUBLIC bicycl)
add_dependencies (examples chia-vdf)

add_executable (CL_HSMqk_setup CL_HSMqk_setup.cpp)
target_link_libraries (CL_HSMqk_setup PUBLIC bicycl)
add_dependencies (examples CL_HSMqk_setup)

add_executable (CL_HSM2k_setup CL_HSM2k_setup.cpp)
target_link_libraries (CL_HSM2k_setup PUBLIC bicycl)
add_dependencies (examples CL_HSM2k_setup)
