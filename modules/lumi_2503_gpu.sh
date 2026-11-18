#!/bin/bash

module load LUMI/25.03
module load partition/G
# HIP needs a Clang-based compiler behind CC; pin it instead of relying on the login default.
module load cpeCray/25.03
module load rocm/6.3.4
module load papi/7.2.0.1

export MPICH_GPU_SUPPORT_ENABLED=1
# Binaries are built for gfx90a:xnack+.
export HSA_XNACK=1
