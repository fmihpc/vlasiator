#!/bin/bash

module purge
module load GCC/13.2.0
module load OpenMPI/4.1.6-GCC-13.2.0
module load PMIx/4.2.6-GCCcore-13.2.0
module load PAPI/7.1.0-GCCcore-13.2.0
module load Boost/1.83.0-GCC-13.2.0
module load CMake/3.27.6-GCCcore-13.2.0
export VLASIATOR_ARCH=carrington_gcc_openmpi

