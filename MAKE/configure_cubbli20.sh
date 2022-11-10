#!/bin/bash
# This script is assumed to be loacted in MAKE folder

MAKE_FOLDER=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Create external_libs folder

NAME_OF_EXTERNAL_LIBS_FOLDER="external_libs"

mkdir -p $MAKE_FOLDER/../$NAME_OF_EXTERNAL_LIBS_FOLDER

EXTERNAL_LIBS_FOLDER=$( cd $(dirname ${MAKE_FOLDER}); cd $MAKE_FOLDER/../$NAME_OF_EXTERNAL_LIBS_FOLDER; pwd -P )

# Libs are installed in different subshells

    # Zoltan
(
    cd $EXTERNAL_LIBS_FOLDER
    git clone https://github.com/sandialabs/Zoltan
    cd Zoltan
    mkdir build
    cd build
    ../configure --prefix="$EXTERNAL_LIBS_FOLDER/Zoltan/build" --enable-mpi --with-mpi-compilers=yes --with-gnumake --with-id-type=ullong CC=mpicc CXX=mpicxx
    make -j 4
    make install
)
    # VLSV 
(
    cd $EXTERNAL_LIBS_FOLDER
    git clone https://github.com/fmihpc/vlsv
    cd vlsv
    ARCH=arch make

)
    # DCCRG
(
    cd $EXTERNAL_LIBS_FOLDER
    git clone https://github.com/fmihpc/dccrg
    cd dccrg
    git checkout vlasiator-version
)
    
    # FSGRID
(
    cd $EXTERNAL_LIBS_FOLDER
    git clone https://github.com/fmihpc/fsgrid
)
    # PHIPROF
(
    cd $EXTERNAL_LIBS_FOLDER
    git clone https://github.com/fmihpc/phiprof/
    cd phiprof/src
    make

)
    # Eigen
(
    cd $EXTERNAL_LIBS_FOLDER
    EIGEN_VERSION="3.4.0"
    wget https://gitlab.com/libeigen/eigen/-/archive/$EIGEN_VERSION/eigen-$EIGEN_VERSION.tar.bz2
    tar -xvf eigen-$EIGEN_VERSION.tar.bz2
    rm -f eigen-$EIGEN_VERSION.tar.bz2
)
                                                                                                                      
    # Vectorclass/version1
(
    cd $EXTERNAL_LIBS_FOLDER
    git clone https://github.com/vectorclass/version1
    git clone https://github.com/vectorclass/add-on
    cp add-on/vector3d/vector3d.h version1
    rm -rf add-on
)

