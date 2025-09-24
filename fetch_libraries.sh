#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
else
   PLATFORM=""
fi
echo "Fetching library files for platform $PLATFORM"

mkdir -p library-build
cd library-build

# Phiprof
git clone https://github.com/fmihpc/phiprof/

# VLSV
if [[ $PLATFORM != "-appleM1" ]]; then
   git clone https://github.com/fmihpc/vlsv.git
else
   git clone -b appleM1Build https://github.com/ursg/vlsv.git
fi

# PAPI
if [[ $PLATFORM != "-arriesgado" && $PLATFORM != "-appleM1" && $PLATFORM != "-ukkogpu" && $PLATFORM != "-hile_cpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc"  && $PLATFORM != "-lumi_2403" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" ]]; then
    # This fails on RISCV and MacOS
    # Mahti, LUMI, UkkoGPU and HILE use system module
    git clone https://github.com/icl-utk-edu/papi
fi

# jemalloc (not for GPU versions, on Mahti use system module)
if [[ $PLATFORM != "-leonardo_booster" && $PLATFORM != "-karolina_cuda" && $PLATFORM != "-ukkogpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" ]]; then
    curl -O -L https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2
    tar xjf jemalloc-5.3.0.tar.bz2
fi

# Zoltan
git clone https://github.com/sandialabs/Zoltan.git

# Boost (only if system module not available)
if [[ $PLATFORM == "-leonardo_booster" || $PLATFORM == "-leonardo_dcgp" || $PLATFORM == "-karolina_cuda" || $PLATFORM == "-karolina_gcc" || $PLATFORM == "-ukkogpu" || $PLATFORM == "-mahti_gcc_build" ]]; then
    echo "### Downloading boost. ###"
    wget -q https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
    echo "### Extracting boost. ###"
    tar -xzf boost_1_86_0.tar.gz
fi
