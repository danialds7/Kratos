#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list with all the compilation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC="/globalfs/opt/gcc/10.2.0/bin/gcc"
export CXX="/globalfs/opt/gcc/10.2.0/bin/g++"

# Set variables
export CMAKE_PATH="/globalfs/opt/cmake/3.13.0/bin:${PATH}"
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}
export BOOST_ROOT=${BOOST_ROOT:-"/globalfs/opt/boost/1.78.0"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/globalfs/opt/python/3.12.1/bin/python3"}
export TRILINOS_INCLUDE_DIR=${TRILINOS_INCLUDE_DIR:-"/globalfs/opt/trilinos/12.10.1/include"}
export TRILINOS_LIBRARY_DIR=${TRILINOS_INCLUDE_DIR:-"/globalfs/opt/trilinos/12.10.1/lib"}
# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication;
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
add_app ${KRATOS_APP_DIR}/FreeSurfaceApplication;
add_app ${KRATOS_APP_DIR}/FluidDynamicsHydraulicsApplication;


# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-DUSE_MPI=ON                                                       \
-DUSE_EIGEN_MKL=OFF                                                 \
-DKRATOS_BUILD_TESTING=OFF                                          \
-DKRATOS_GENERATE_PYTHON_STUBS=ON                                    \
-DTRILINOS_ROOT=/globalfs/opt/trilinos/12.10.1/                      \
#-DCMAKE_INSTALL_PREFIX="/home/ddehghan/Kratos"                        \
-DCMAKE_BUILD_TYPE=Release
#-DCMAKE_CXX_FLAGS="-g3" 

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j8