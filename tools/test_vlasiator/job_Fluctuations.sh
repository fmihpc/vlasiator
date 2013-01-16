#!/bin/bash -l
#PBS -l mppwidth=2
#PBS -l mppdepth=12
#PBS -l walltime=00:15:00
#PBS -V  
#PBS -N short_test
##PBS -q short

cd $PBS_O_WORKDIR 

export OMP_NUM_THREADS=12
export MPICH_MAX_THREAD_SAFETY=funneled

aprun -n 2 -d 12 vlasiator --run_config=Fluctuations.cfg
