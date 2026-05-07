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
cd phiprof
git checkout 605a7247c85d967fe22fe079c96c817b461c92b1
cd ..

# VLSV
if [[ $PLATFORM != "-appleM1" ]]; then
   git clone https://github.com/fmihpc/vlsv.git
   cd vlsv
   git checkout 0d06db7078ee7066f69180b559c506c4cb0d7f1b
   cd ..
else
   git clone -b appleM1Build https://github.com/ursg/vlsv.git
   cd vlsv
   git checkout 0d06db7078ee7066f69180b559c506c4cb0d7f1b
   cd ..
fi

# PAPI
if [[ $PLATFORM != "-arriesgado" && $PLATFORM != "-appleM1" && $PLATFORM != "-ukkogpu" && $PLATFORM != "-hile_cpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc"  && $PLATFORM != "-lumi_2403" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" && $PLATFORM != "-frankenstein_hopper2_cuda" && $PLATFORM != "-roihu_cpu" ]]; then
    # This fails on RISCV and MacOS
    # Mahti, LUMI, UkkoGPU and HILE use system module
    git clone https://github.com/icl-utk-edu/papi
    cd papi
    git checkout 25a278ee5f4ccc9a2263e90ff8c15a1a58b2b7ed
    cd ..
fi

# jemalloc (not for GPU versions, on Mahti use system module)
if [[ $PLATFORM != "-leonardo_booster" && $PLATFORM != "-karolina_cuda" && $PLATFORM != "-ukkogpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" && $PLATFORM != "-frankenstein_hopper2_cuda" ]]; then
    curl -O -L https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2
    tar xjf jemalloc-5.3.0.tar.bz2
fi

# Zoltan
git clone https://github.com/sandialabs/Zoltan.git
cd Zoltan
git checkout f6361719dd66cac62db8dbed120704e436a5ee81
cd ..

# Boost (only if system module not available)
if [[ $PLATFORM == "-leonardo_booster" || $PLATFORM == "-leonardo_dcgp" || $PLATFORM == "-karolina_cuda" || $PLATFORM == "-karolina_gcc" || $PLATFORM == "-ukkogpu" || $PLATFORM == "-mahti_gcc_build" || $PLATFORM == "-frankenstein_hopper2_cuda" && $PLATFORM != "-roihu_cpu" ]]; then
    echo "### Downloading boost. ###"
    wget -q https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
    echo "### Extracting boost. ###"
    tar -xzf boost_1_86_0.tar.gz
fi
