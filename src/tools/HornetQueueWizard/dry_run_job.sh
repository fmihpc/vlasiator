#!/bin/bash
#PBS -N general_400
#PBS -l nodes=400
#PBS -l walltime=24:00:00
#PBS -W umask=007
ht=2                   #hyper threads per physical core
t=6                   #threads per process

export DRY_RUN=1
export PBS_O_WORKDIR=$(pwd)
export PBS_JOBID=dryrun

# set manually these
nodes=400

# Set job list to choose from
JOBLIST=joblist_dryrun.txt

#hornet has 2 x 12 cores
cores_per_node=24
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t
umask 007

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Determine the number of processes
NUM_PROCESSES=$(( $nodes * $units_per_node / $OMP_NUM_THREADS ))
# Set the number of restarts to write
NUMBER_OF_RESTARTS=2
# Allow for this buffer for the restart IO (in s)
RESTART_IO_EXTRA_TIME=1800.0
# Determine walltime
WALLTIME=1

source vlasiator_run_next.sh

