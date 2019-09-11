#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --job-name=testpackage
#SBATCH --account=project_2000203
#SBATCH --partition=small
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000

ht=1                   #hyper threads per physical core, can only be 1
t=5                   #threads per process

#Compute and set  stuff, do not change
if [ -z $SLURM_NNODES ]
then
    #if running interactively we use 2 nodes
    nodes=2
else
    nodes=$SLURM_NNODES
fi

#sisu has 2 x 12 cores
cores_per_node=40
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t


module load gcc
module load boost
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/proj/vlasiato/libraries/taito/openmpi/1.10.2/gcc/4.9.3/papi/5.5.0/lib/


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
create_verification_files=1


#folder for all reference data 
reference_dir="/scratch/project_2000203/testpackage_data"
#compare agains which revision
reference_revision="current"



# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
