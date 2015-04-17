#!/bin/bash
#PBS -N general_1
#PBS -l nodes=1
#PBS -l walltime=1:00:00
#PBS -W umask=007
ht=2                   #hyper threads per physical core
t=6                   #threads per process


#Compute and set  stuff, do not change
nodes=$( qstat -f $PBS_JOBID | grep   Resource_List.nodes | gawk '{print $3}' )
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

# Set job list to choose from
JOBLIST=/zhome/academic/HLRS/pri/ipryakem/joblist_1.txt

# Determine the number of processes
NUM_PROCESSES=$(( $(qstat -f $PBS_JOBID | grep "Resource_List.nodect" | cut -d "=" -f 2) * $units_per_node / $OMP_NUM_THREADS ))
# Set the number of restarts to write
NUMBER_OF_RESTARTS=1
# Allow for this buffer for the restart IO (in s)
RESTART_IO_EXTRA_TIME=3400.0
# Determine walltime
WALLTIME=$(qstat -f $PBS_JOBID | grep "Resource_List.walltime" | cut -d "=" -f 2 | cut -d ":" -f 1)

source vlasiator_run_next.sh

