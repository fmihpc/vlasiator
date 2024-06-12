#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
else 
   PLATFORM=""
fi
echo "Using platform $PLATFORM"

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

if [[ $PLATFORM == "-arriesgado" ]]; then
   # Special workaround for missing include paths on arriesgado
   make -j 4 CCC=mpic++ CCFLAGS="-I /usr/lib/gcc/riscv64-linux-gnu/11/include -fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -fopenmp -W -Wall -Wextra -pedantic"
elif [[ $PLATFORM == "-appleM1" ]]; then
   make -j 4 CCC=mpic++ CC=appleLLVM CCFLAGS="-fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -fopenmp" LDFLAGS="-fopenmp"
else
   make -j 4 CCC=mpic++
fi
cp ../include/* $WORKSPACE/libraries${PLATFORM}/include
cp ../lib/* $WORKSPACE/libraries${PLATFORM}/lib
cd ../..

# Build VLSV
if [[ $PLATFORM != "-appleM1" ]]; then
   git clone https://github.com/fmihpc/vlsv.git
else
   git clone -b appleM1Build https://github.com/ursg/vlsv.git
fi
cd vlsv
make
cp libvlsv.a $WORKSPACE/libraries${PLATFORM}/lib
cp *.h $WORKSPACE/libraries${PLATFORM}/include
cd ..

# Build papi
if [[ $PLATFORM != "-arriesgado" && $PLATFORM != "-appleM1" ]]; then  # This fails on RISCV and MacOS
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
if [[ $PLATFORM == "-arriesgado" ]]; then
   ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong --host=riscv64-unknown-linux-gnu --build=arm-linux-gnu && make -j 4 && make install
elif [[ $PLATFORM == "-appleM1" ]]; then
   ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong CC=mpicc CXX=mpic++ && make -j 4 && make install
else
   ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong && make -j 4 && make install
cd ..
fi
