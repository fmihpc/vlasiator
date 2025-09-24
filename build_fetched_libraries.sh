#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
else
   PLATFORM=""
fi
echo "Using platform $PLATFORM"

# Clean up old libraries for this platform
rm -rf libraries${PLATFORM}

# Create new ones
mkdir -p libraries${PLATFORM}/include
mkdir -p libraries${PLATFORM}/lib

# Assumes required files are available in this directory
cd library-build

# Some platforms allow for nice parallel builds.
if [[ $PLATFORM == "-pioneer" ]]; then
   PARALLEL=64
elif [[ $PLATFORM == "-hile_cpu" || $PLATFORM == "-hile_gpu" || $PLATFORM == "-lumi_2403" ]]; then
   PARALLEL=128
else
   # Otherwise we are friendly to other users and limit our parallelism
   PARALLEL=4
fi



# Build phiprof
#git clone https://github.com/fmihpc/phiprof/
cd phiprof/src
make clean
if [[ $PLATFORM == "-pioneer" ]]; then
   # Special workaround for missing include paths on pioneer
   make -j $PARALLEL CCC=mpic++ CCFLAGS="-I /usr/lib/gcc/riscv64-linux-gnu/11/include -fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -fopenmp -W -Wall -Wextra -pedantic"
elif [[ $PLATFORM == "-appleM1" ]]; then
   make -j $PARALLEL CCC=mpic++ CC=appleLLVM CCFLAGS="-fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -fopenmp" LDFLAGS="-fopenmp"
elif [[ $PLATFORM == "-leonardo_dcgp_intel" ]]; then
   make -j $PARALLEL CCC="mpiicpc -cxx=icpx" CC="mpiicc -cc=icx" CCFLAGS="-fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -qopenmp" LDFLAGS="-qopenmp"
elif [[ $PLATFORM == "-hile_cpu" || $PLATFORM == "-hile_gpu" || $PLATFORM == "-lumi_2403" ]]; then
   make -j $PARALLEL CCC=CC CC=cc CCFLAGS="-fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -fopenmp" LDFLAGS="-fopenmp"
elif [[ $PLATFORM == "-lumi_hipcc" ]]; then
    make -j $PARALLEL CC=hipcc
    CCFLAGS="-fpic -O2 -std=c++17 -DCLOCK_ID=CLOCK_MONOTONIC -fopenmp" LDFLAGS="-lstdc++ -fopenmp -lgomp -lroctx64 -lroctracer64"
else
   make -j $PARALLEL CCC=mpic++
fi
cp ../include/* $WORKSPACE/libraries${PLATFORM}/include
cp ../lib/* $WORKSPACE/libraries${PLATFORM}/lib
cd ../..

# Build VLSV
# if [[ $PLATFORM != "-appleM1" ]]; then
#    git clone https://github.com/fmihpc/vlsv.git
# else
#    git clone -b appleM1Build https://github.com/ursg/vlsv.git
# fi
cd vlsv
make clean ARCH=arch
if [[ $PLATFORM == "-leonardo_dcgp_intel" ]]; then
   make ARCH=arch CMP="mpiicpc -cxx=icpx"
elif [[ $PLATFORM == "-hile_cpu" || $PLATFORM == "-hile_gpu" || $PLATFORM == "-lumi_hipcc" || $PLATFORM == "-lumi_2403" ]]; then
   make ARCH=arch CMP=CC
else
   make ARCH=arch
fi
cp libvlsv.a $WORKSPACE/libraries${PLATFORM}/lib
cp *.h $WORKSPACE/libraries${PLATFORM}/include
cd ..

# Build papi
if [[ $PLATFORM != "-pioneer" && $PLATFORM != "-appleM1" && $PLATFORM != "-ukkogpu" && $PLATFORM != "-hile_cpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc" && $PLATFORM != "-lumi_2403" && $PLATFORM != "-mahti_gcc_build" && $PLATFORM != "-mahti_cuda" ]]; then
    # This fails on RISCV and MacOS
    # LUMI, UkkoGPU and HILE use system module
    # git clone https://github.com/icl-utk-edu/papi
    cd papi/src
    if [[ $PLATFORM == "-leonardo_dcgp_intel" ]]; then
        # OneAPI compilers should use CC="mpiicc -cc=iccx" but this fails in configure. Needed to modify few files to pass configure and compilation phases
        sed -i 's/(MAKE) CC=$(CC)/(MAKE) CC="$(CC)"/g' Makefile.inc
        sed -i 's/DBG?=-g -Wall -Werror -Wextra -Wno-unused-parameter/DBG?=-g -Wall -Wextra/g' libpfm4/config.mk
        ./configure --prefix=$WORKSPACE/libraries${PLATFORM} CC="mpiicc -cc=icx" CXX="mpiicpc -cxx=icpx" MPICC="mpiicc -cc=icx"
    else
        ./configure --prefix=$WORKSPACE/libraries${PLATFORM} CC=mpicc CXX=mpic++
    fi
    make clean
    make -j $PARALLEL && make install
    cd ../..
fi

# Build jemalloc (not for GPU versions or Mahti)
if [[ $PLATFORM != "-leonardo_booster" && $PLATFORM != "-karolina_cuda" && $PLATFORM != "-ukkogpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" ]]; then
    # curl -O -L https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2
    # tar xjf jemalloc-5.3.0.tar.bz2
    cd jemalloc-5.3.0
    if [[ $PLATFORM == "-pioneer" ]]; then
        ./configure --prefix=$WORKSPACE/libraries${PLATFORM} --with-jemalloc-prefix=je_
    elif [[ $PLATFORM == "-leonardo_dcgp_intel" ]]; then
        ./configure --prefix=$WORKSPACE/libraries${PLATFORM} --with-jemalloc-prefix=je_ CC="mpiicc -cc=icx" CXX="mpiicpc -cxx=icpx"
    elif [[ $PLATFORM == "-hile_cpu" || $PLATFORM == "-lumi_2403" ]]; then
        ./configure --prefix=$WORKSPACE/libraries${PLATFORM} --with-jemalloc-prefix=je_ CC="cc" CXX="CC"
    else
        ./configure --prefix=$WORKSPACE/libraries${PLATFORM} --with-jemalloc-prefix=je_ CC=mpicc CXX=mpic++
    fi
    make clean
    make -j $PARALLEL
    make install
    cd ..
fi

# Build Zoltan
# git clone https://github.com/sandialabs/Zoltan.git
rm -rf zoltan-build
mkdir zoltan-build
cd zoltan-build
if [[ $PLATFORM == "-pioneer" ]]; then
    ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong --host=riscv64-unknown-linux-gnu --build=arm-linux-gnu
elif [[ $PLATFORM == "-appleM1" || $PLATFORM == "-meluxina" ]]; then
    ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong CC=mpicc CXX=mpic++
elif [[ $PLATFORM == "-leonardo_dcgp_intel" ]]; then
    ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong CC="mpiicc -cc=icx" CXX="mpiicpc -cxx=icpx"
    # Although configured with new compilers, the compilations ignores the -cc=icx and -cxx=icpx flags. Need to add them manually.
    sed -i 's/mpiicc/mpiicc -cc=icx/g' Makefile
    sed -i 's/mpiicc/mpiicc -cc=icx/g' src/Makefile
    sed -i 's/mpiicc/mpiicc -cc=icx/g' src/driver/Makefile
    sed -i 's/mpiicpc/mpiicpc -cxx=icpx/g' Makefile
    sed -i 's/mpiicpc/mpiicpc -cxx=icpx/g' src/Makefile
    sed -i 's/mpiicpc/mpiicpc -cxx=icpx/g' src/driver/Makefile
elif [[ $PLATFORM == "-hile_cpu" ||  $PLATFORM == "-hile_gpu" || $PLATFORM == "-lumi_2403" ]]; then
   ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong CC=cc CXX=CC
else
    ../Zoltan/configure --prefix=$WORKSPACE/libraries${PLATFORM} --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ullong CC=mpicc CXX=mpic++
fi
make clean
make -j $PARALLEL && make install
cd ..

# Build boost (only if system module is not available)
if [[ $PLATFORM == "-leonardo_booster" || $PLATFORM == "-leonardo_dcgp" || $PLATFORM == "-karolina_cuda" || $PLATFORM == "-karolina_gcc" || $PLATFORM == "-ukkogpu" ||  $PLATFORM == "-mahti_gcc_build" ]]; then
    # echo "### Downloading boost. ###"
    # wget -q https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
    # echo "### Extracting boost. ###"
    # tar -xzf boost_1_86_0.tar.gz
    echo "### Building boost. ###"
    cd boost_1_86_0
    ./bootstrap.sh --with-libraries=program_options --prefix=$WORKSPACE/libraries${PLATFORM} stage
    echo "using mpi ;" >> ./tools/build/src/user-config.jam
    ./b2
    echo "### Installing boost. ###"
    ./b2 --prefix=$WORKSPACE/libraries${PLATFORM} install > /dev/null
    cd ..
fi

# Clean up build directory
#rm -rf $BUILDDIR
#
