# Installation {#install}

## Download

The code source is available at https://gite.lirmm.fr/crypto/bicycl.

## Build from source

To compile the code, a C++ compiler and CMake 3.5.1 or later are necessary, and
the following libraries are required:
  - GMP
  - openSSL

\note On Debian and Ubuntu, the necessary files can be installed with

    apt install g++ libgmp-dev libssl-dev cmake


Compilation can be done by running the following commands:

    mkdir build
    cd build
    cmake ..

\note The directory name `build` can be replaced by any name that do not already
exist. This documentation assumes that the name `build` was used.

The following arguments can be added to the `cmake ..` command-line to customize
the build:

  - `CMAKE_INSTALL_PREFIX`: a path to the desired install directory (default
    install directory is chosen by CMake). For example,

        cmake -DCMAKE_INSTALL_PREFIX=/alternative/install/path ..

  - `WITH_THREADS`: a boolean to specify whether or not to use multithread code
    (default is ON). For example, 

        cmake -DWITH_THREADS=OFF ..

  - `GMP_DIR`: a path to help CMake find the GMP library.

        cmake -DGMP_DIR=/path/to/gmp ..

  - `OPENSSL_ROOT_DIR`: a path to help CMake find the OpenSSL library. For,
    example,

        cmake -DOPENSSL_ROOT_DIR=/path/to/openssl ..


## Installation

To install the library, run the following command in the `build` directory:

    make install


Depending on where the install directory is located, you may need to prefix the
command with sudo for the installation to succeed:

    sudo make install

## Examples

Some examples can be found in the `examples` directory of the source tree. To
compile all the examples, run the following command in the `build` directory:

    make examples

To compile only a specific example (`chia-vdf` in this example), run the
following command in the `build` directory:

    make chia-vdf

## Tests

To test the code, run the following command in the `build` directory:

    make check

## Benchmarks

To run the benchmarks, run the following command in the `build` directory:

    make benchs


To build the benchmarks (but not run them), run the following command in the
`build` directory:

    make benchs_build

