#!/bin/bash
# nvc++ is only available via nvhpc; it bundles its own CUDA 13.1, gcc
# toolchain, and HPCX MPI (do not also load a separate openmpi module - that
# would mix two MPI/OpenMP runtime implementations in the same binary).
# No --gcc-toolchain= override needed - confirmed working against the plain
# module on the real system (an earlier concern that nvc++ couldn't find a
# usable C++ standard library without one turned out to be specific to this
# agent's own containerized test environment, not the real Roihu environment).
module load nvhpc/26.3
