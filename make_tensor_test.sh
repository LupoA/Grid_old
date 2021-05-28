#!/bin/bash

# This should ideally integrate with the make system rather than being a shell script
# But since I don't understand GNU autotools, this is a quick and dirty solution
# This should be removed before pull requesting (either by tidying up the unit tests and
# integrating them with the Grid build system, or by stripping them out entirely)


# Beyond Grid's usual dependencies, the tests require Boost to be installed.
# Configure Grid in the usual way (e.g. mkdir build && cd build && ../configure).
# Set where you have OpenSSL and Boost installed, if they are not in the system path.
# Call this script from within the Grid subdirectory of your build directory.


CXX=g++
BOOST_PREFIX=/opt/homebrew
OPENSSL_PREFIX=/opt/homebrew/Cellar/openssl@1.1/1.1.1k
WARNINGS="-Wall -Wno-unused-private-filed -Wno-unused-local-typedef -Wno-unused-variable"
INCLUDE="-I${BOOST_PREFIX}/include -I${OPENSSL_PREFIX}/include -I$(git rev-parse --show-toplevel) -I."
LIBDIRS="-L${BOOST_PREFIX}/lib -L${OPENSSL_PREFIX}/lib"
LIBS="-lboost_unit_test_framework"

SOURCE=../../Grid/tensors/tests/test_imatrix.cc
OBJ=tensors/tests/test_imatrix

# Stop if anything breaks
set -eux

# Ensure directory to contain binary exists
if [ ! -d tensors/tests ]; then
    mkdir tensors/tests
fi

# Do compilation
${CXX} -std=c++11 ../../Grid/tensors/tests/test_imatrix.cc ${INCLUDE} ${LIBDIRS} ${LIBS} -o tensors/tests/test_imatrix

# Run tests
./tensors/tests/test_imatrix
