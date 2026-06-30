#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
else
   PLATFORM=""
fi
echo "Fetching library files for platform $PLATFORM"

PHIPROF_COMMIT="fb9b37b3d7ce592069143de465c53caf5fb80e1d"
VLSV_COMMIT="0d06db7078ee7066f69180b559c506c4cb0d7f1b"
PAPI_COMMIT="25a278ee5f4ccc9a2263e90ff8c15a1a58b2b7ed"

TRILINOS_BRANCH="zoltanLBSafeAllreduce-issue15235"
TRILINOS_COMMIT="16ceeebdfbe0809a549e0543f2824a79ebc2aa2d"

ZFP_COMMIT="f2046180a8fea296646236d7d612d89b52841d46"
TUCKER_OCTREE_COMMIT="3bc42470a5b486d947005bebf03d8846a4af9aa4"

git_use_commit() {
	git fetch origin "$1"
	git checkout "$1"
}

mkdir -p library-build
cd library-build

# Phiprof
git clone --depth=1 https://github.com/fmihpc/phiprof/
cd phiprof
git_use_commit $PHIPROF_COMMIT
cd ..

# VLSV
if [[ $PLATFORM != "-appleM1" ]]; then
   git clone --depth=1 https://github.com/fmihpc/vlsv.git
   cd vlsv
   git_use_commit $VLSV_COMMIT
   cd ..
else
   git clone --depth=1 -b appleM1Build https://github.com/ursg/vlsv.git
   cd vlsv
   git_use_commit $VLSV_COMMIT
   cd ..
fi

# PAPI
if [[ $PLATFORM != "-arriesgado" && $PLATFORM != "-appleM1" && $PLATFORM != "-ukko_dgx" && $PLATFORM != "-hile_cpu" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc"  && $PLATFORM != "-lumi_2403" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" && $PLATFORM != "-frankenstein_hopper2_cuda" && $PLATFORM != "-roihu_cpu" && $PLATFORM != "-roihu_cpu_aocc" && $PLATFORM != "-roihu_gpu" ]]; then
    # This fails on RISCV and MacOS
    # Mahti, LUMI, UkkoGPU and HILE use system module
    git clone --depth=1 https://github.com/icl-utk-edu/papi
    cd papi
    git_use_commit "$PAPI_COMMIT"
    cd ..
fi

# jemalloc (not for GPU versions, on Mahti use system module)
if [[ $PLATFORM != "-leonardo_booster" && $PLATFORM != "-karolina_cuda" && $PLATFORM != "-ukko_dgx" && $PLATFORM != "-hile_gpu" && $PLATFORM != "-lumi_hipcc" && $PLATFORM != "-mahti_cuda" && $PLATFORM != "-mahti_gcc_build" && $PLATFORM != "-frankenstein_hopper2_cuda" && $PLATFORM != "-roihu_gpu" ]]; then
    curl -O -L https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2
    tar xjf jemalloc-5.3.0.tar.bz2
fi

# Zoltan
git clone --depth=1 --branch="$TRILINOS_BRANCH" https://github.com/ykempf/Trilinos.git
cd Trilinos
git_use_commit $TRILINOS_COMMIT
cd ..

#ZFP and OCTREE
git clone --depth=1 https://github.com/LLNL/zfp.git
cd zfp
git_use_commit $ZFP_COMMIT
cd ..

git clone --depth=1 https://github.com/cschpc/tucker-octree.git
cd tucker-octree
git_use_commit $TUCKER_OCTREE_COMMIT
cd ..
