#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -J small_test
#SBATCH -p test
#SBATCH -N 2
#SBATCH --constraint=hsw
#SBATCH --no-requeue          

ht=1                   #hyper threads per physical core, can only be 1
t=3                   #threads per process

#Compute and set  stuff, do not change
if [ -z $SLURM_NNODES ]
then
    #if running interactively we use 2 nodes
    nodes=2
else
    nodes=$SLURM_NNODES
fi

#sisu has 2 x 12 cores
cores_per_node=24
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t


module purge
module load gcc
module load openmpi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/proj/vlasiato/libraries/taito/openmpi/1.10.2/gcc/4.9.3/boost/1.61.0/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/proj/vlasiato/libraries/taito/openmpi/1.10.2/gcc/4.9.3/jemalloc/4.0.4/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/proj/vlasiato/libraries/taito/openmpi/1.10.2/gcc/4.9.3/papi/5.5.0/lib/


umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"
#command for running stuff
run_command="srun"
small_run_command="srun -n 1"
run_command_tools="srun -n 1"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0


#folder for all reference data 
reference_dir="/proj/vlasiato/testpackage/"
#compare agains which revision
reference_revision="c36241b84ce8179f7491ebf2a94c377d7279e8c9__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"



# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
