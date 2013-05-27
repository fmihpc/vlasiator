#!/bin/bash
#PBS -N c2048_24h
# Request the number of cores that you need in total
#PBS -l mppwidth=2048
#PBS -l mppnppn=32
# Request the time you need for computation
#PBS -l walltime=24:00:00
#PBS -W umask=002

umask 002

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Set job list to choose from
JOBLIST=joblist_2048_24h.txt
# Set the number of OpenMP threads per node
export OMP_NUM_THREADS=8
# Determine the number of processes
NUM_PROCESSES=$(qstat -f $BATCH_JOBID | grep mppwidth | cut -d " " -f 7) / $OMP_NUM_THREADS
# Set the number of restarts to write
NUMBER_OF_RESTARTS=2
# Allow for this buffer for the restart IO (in s)
RESTART_IO_EXTRA_TIME=1800.0
# Determine walltime
WALLTIME=$(qstat -f $BATCH_JOBID | grep walltime | cut -d " " -f 7 | cut -d ":" -f 1)

source vlasiator_run_next.sh
