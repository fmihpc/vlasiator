#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
else 
   PLATFORM=""
fi

# Clean up old library build dirs and libraries for this platform
rm -rf library-build libraries${PLATFORM}

# Create new ones
mkdir -p libraries${PLATFORM}/include
mkdir -p libraries${PLATFORM}/lib
mkdir library-build
cd library-build

# Build phiprof
git clone https://github.com/fmihpc/phiprof/ 
cd phiprof/src
make -j 4 CCC=mpic++
cp ../include/* $WORKSPACE/libraries${PLATFORM}/include
cp ../lib/* $WORKSPACE/libraries${PLATFORM}/lib
cd ../..

# Build VLSV
git clone https://github.com/fmihpc/vlsv.git  
cd vlsv
make
cp libvlsv.a $WORKSPACE/libraries${PLATFORM}/lib
cp *.h $WORKSPACE/libraries${PLATFORM}/include
cd ..

# Build papi
if [[ $PLATFORM != "-arriesgado" ]]; then  # This fails on RISCV
   git clone https://github.com/icl-utk-edu/papi
   cd papi/src
   ./configure --prefix=$WORKSPACE/libraries${PLATFORM} && make -j 4 CC=gcc && make install
   cd ../..
fi

# Build jemalloc
wget https://github.com/jemalloc/jemalloc/releases/download/4.0.4/jemalloc-4.0.4.tar.bz2
tar xjf jemalloc-4.0.4.tar.bz2
cd jemalloc-4.0.4
./configure --prefix=$WORKSPACE/libraries${PLATFORM} --with-jemalloc-prefix=je_ && make -j 4 && make install
cd ..

# Build Zoltan
git clone https://github.com/sandialabs/Zoltan.git
mkdir zoltan-build
cd zoltan-build
if [[ $PLATFORM != "-arriesgado" ]]; then
   ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong && make -j 4 && make install
else
   ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong --build=arm-linux-gnu && make -j 4 && make install
cd ..
fi

git clone https://gitlab.com/libeigen/eigen.git
cp -ua eigen/Eigen $WORKSPACE/libraries${PLATFORM}/include
