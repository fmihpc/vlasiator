#!/bin/bash

# GCC 14.3.0 + CUDA 12.9.1: the newest pair that supports sm_70 without -allow-unsupported-compiler.
module purge
module load PrgEnv-gnu-v100/14.3.0
module load CUDA/12.9.1
module load MPICH/5.0.1-CUDA-12.9.1
module load bzip2
module load CMake
module load HY_proxy
module load Automake
# no module load git because that only exists for gcc 15. but the system git is new enough
