#!/bin/bash

export HTTP_PROXY=http://www-cache.cs.helsinki.fi:3128
export HTTPS_PROXY=http://www-cache.cs.helsinki.fi:3128
export https_proxy=http://www-cache.cs.helsinki.fi:3128
export http_proxy=http://www-cache.cs.helsinki.fi:3128

export MPICH_GPU_SUPPORT_ENABLED=1
module load papi
module load rocm
module load cray-pmi
module load craype-accel-amd-gfx90a
module load libfabric

module use /appl/hile/modules
module load CMake/3.27.7

