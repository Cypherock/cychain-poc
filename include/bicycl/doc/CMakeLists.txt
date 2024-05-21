configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

find_program (DOXYGEN doxygen)
if (NOT DOXYGEN)
  message (STATUS "Warning, will not be able to build doc: doxygen not found")
  set (MSG "The documentation cannot be build as CMake failed to find doxygen \
during the configuration step. If you want to build the documentation, you \
need to install doxygen and re-run the configuration step of CMake.")
  add_custom_target (doc COMMAND ${CMAKE_COMMAND} -E echo ${MSG})
else ()
  add_custom_target (doc COMMAND ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
endif ()
