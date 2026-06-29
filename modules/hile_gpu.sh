#!/bin/bash

export MPICH_GPU_SUPPORT_ENABLED=1
module load papi
module load rocm/6.3.4
module load cray-pmi
module load craype-accel-amd-gfx90a
module load libfabric/1.22.0

module use /appl/hile/modules
module load CMake/3.27.7

