#!/bin/bash

# Performs a clean install of Vlasiator, including libraries and tools.
# Made to be able to be thrown at e.g. UH compute nodes with srun, for example
# export VLASIATOR_ARCH=turso-amd_GNU_MPICH
# srun --constraint=amd --mem=256G -c 32 -t 00:30:00 clean_install.sh
# or
# srun --constraint=amd --mem=256G -c 32 -t 00:30:00 clean_install.sh turso-amd_GNU_MPICH

set -e

if [[ z$VLASIATOR_ARCH == "z" ]]
then
   if [[ z$1 != "z" ]]
   then
     export VLASIATOR_ARCH=$1 
     echo "using VLASIATOR_ARCH=$VLASIATOR_ARCH"
   else
     echo "No VLASIATOR_ARCH set or given"
     echo "usage: ./clean_install.sh VLASIATOR_ARCH"
     echo "With e.g. turso-amd_GNU_MPICH as VLASIATOR_ARCH, either as parameter or as an environment variable"
     exit 1
   fi
else
   echo "using VLASIATOR_ARCH=$VLASIATOR_ARCH"
fi

rm -rf submodules library-build "libraries-$VLASIATOR_ARCH"

. ./modules/modules.sh

./fetch_and_build_libraries.sh

make clean
make -j 32 testpackage


make -j 32 tools
cp vlsvdiff_DP ../

echo "Compiled on host: $(hostname)"
module list

echo "LD_LIBRARY_PATH was $LD_LIBRARY_PATH"
echo "Inspecting ldd vlasiator"
ldd vlasiator

