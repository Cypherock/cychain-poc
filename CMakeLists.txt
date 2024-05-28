cmake_minimum_required(VERSION 3.5.1)
project(cychain-poc)

# Include the BICYCL library
add_subdirectory("include/bicycl" EXCLUDE_FROM_ALL)

# Add custom target for poc
add_custom_target(poc COMMAND ${CMAKE_COMMAND} -E sleep 0)

# Add the executable for cychain_sig_gen
add_executable(cychain_sig_gen src/cychain_sig_gen.cpp src/ThreadPool.cpp src/CryptoUtils.cpp)

# Include directories
include_directories(include/trezor-firmware/crypto)

# Ensure BICYCL is linked
target_link_libraries(cychain_sig_gen PUBLIC bicycl)

# Set the directory for trezor-firmware
set(TREZOR_FIRMWARE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/trezor-firmware")

# Add custom target to build the trezor-crypto library
add_custom_target(
    trezor_crypto_build
    COMMAND make -C ${TREZOR_FIRMWARE_DIR}/crypto all
    COMMENT "Building trezor-firmware/crypto library"
)

# Add dependency to ensure trezor_crypto_build is built first
add_dependencies(cychain_sig_gen trezor_crypto_build)

# Link the trezor-crypto library
target_link_libraries(cychain_sig_gen PRIVATE "${TREZOR_FIRMWARE_DIR}/crypto/tests/libtrezor-crypto.so")

# Add dependencies for the poc target
add_dependencies(poc cychain_sig_gen)
