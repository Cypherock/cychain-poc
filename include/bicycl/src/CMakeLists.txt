set (BICYCL_MAIN_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/bicycl.hpp)
set (BICYCL_INSTALLABLE_HEADERS_AND_INL
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/seclevel.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/seclevel.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/gmp_extras.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/gmp_extras.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/openssl_wrapper.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/openssl_wrapper.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/ec.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/ec.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/qfi.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/qfi.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/CL_HSM_utils.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/CL_HSM_utils.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/CL_HSMqk.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/CL_HSMqk.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/CL_HSM2k.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/CL_HSM2k.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/Paillier.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/Paillier.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/Joye_Libert.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/Joye_Libert.inl
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/threshold_ECDSA.hpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/bicycl/threshold_ECDSA.inl
                      )
add_library (bicycl INTERFACE)
target_sources (bicycl INTERFACE ${BICYCL_MAIN_HEADER}
                                    ${BICYCL_INSTALLABLE_HEADERS_AND_INL})
target_include_directories (bicycl INTERFACE
                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries (bicycl INTERFACE gmp::gmp OpenSSL::Crypto)
if (THREADS_FOUND AND CMAKE_USE_PTHREADS_INIT)
  target_link_libraries (bicycl INTERFACE Threads::Threads)
  target_compile_definitions (bicycl INTERFACE BICYCL_WITH_PTHREADS)
endif ()

install (TARGETS bicycl EXPORT bicyclExport
                          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (FILES ${BICYCL_MAIN_HEADER}
                          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install (FILES ${BICYCL_INSTALLABLE_HEADERS_AND_INL}
                          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bicycl)
