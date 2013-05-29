#!/bin/bash
#PBS -N c16384_4h
#PBS -j oe
# Request the number of cores that you need in total
#PBS -l mppwidth=16384
#PBS -l mppnppn=32
# Request the time you need for computation
#PBS -l walltime=4:00:00
#PBS -W umask=002

umask 002

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Set job list to choose from
JOBLIST=${PBS_O_WORKDIR}/joblist_16384.txt
# Set the number of OpenMP threads per node
export OMP_NUM_THREADS=8
# Determine the number of processes
NUM_PROCESSES=$(( $(qstat -f $BATCH_JOBID | grep "Resource_List.mppwidth" | cut -d "=" -f 2) / $OMP_NUM_THREADS ))
# Set the number of restarts to write
NUMBER_OF_RESTARTS=2
# Allow for this buffer for the restart IO (in s)
RESTART_IO_EXTRA_TIME=900.0
# Determine walltime
WALLTIME=$(qstat -f $BATCH_JOBID | grep "Resource_List.walltime" | cut -d "=" -f 2 | cut -d ":" -f 1)

source vlasiator_run_next.sh
