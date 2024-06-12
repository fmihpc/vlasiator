#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --job-name=Vlasiator_tp
#SBATCH --account=project_2000203
#SBATCH --partition=test
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --hint=multithread

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores

module purge
module load StdEnv gcc openmpi boost papi jemalloc

ulimit -c unlimited

ht=2    #hyper threads per physical core
t=$OMP_NUM_THREADS     #threads per process
nodes=$SLURM_NNODES

#mahti has 2 x 64 cores
cores_per_node=128
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')

# hint from George and https://github.com/openucx/ucx/issues/5504
export UCX_TLS=ud,ud_v
export OMPI_MCA_coll=^hcoll

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

bin="../vlasiator"
diffbin="/projappl/project_2000203/vlsvdiff_DP"


#folder for all reference data 
reference_dir="/scratch/project_2000203/testpackage_data"
#compare agains which revision
reference_revision="current"



# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
