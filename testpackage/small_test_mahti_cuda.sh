#!/bin/bash -l
##SBATCH --time=0:15:00
##SBATCH --partition=gputest
#SBATCH --time=12:00:00
#SBATCH --partition=gpusmall
##SBATCH --partition=gpumedium
# gpusmall allows up to 2 gpus, medium up to 4
#SBATCH --cpus-per-task=32
#SBATCH --job-name=vlasigputp
#SBATCH --account=project_2004522
##SBATCH --gres=gpu:a100:4
##SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks-per-node=1
##SBATCH --gres=gpu:a100:2
##SBATCH --ntasks-per-node=2
#SBATCH --nodes=1

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores

module purge
module load gcc/10.4.0 openmpi/4.1.5-cuda cuda/12.1.1 boost/1.82.0-mpi papi/7.1.0

ulimit -c unlimited

ht=1    #hyper threads per physical core
t=$OMP_NUM_THREADS     #threads per process
echo "threads ${t}"
nodes=$SLURM_NNODES

#mahti has 2 x 64 cores
cores_per_node=128
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
#override
tasks=1
tasks_per_node=1

# hint from George and https://github.com/openucx/ucx/issues/5504
export UCX_TLS=ud,ud_v
export OMPI_MCA_coll=^hcoll

umask 007

# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core) and $tasks_per_node tasks per node"
#command for running stuff (add compute-sanitizer here if necessary)

#run_command="srun compute-sanitizer "
run_command="srun "
small_run_command="srun -n 1 "
#small_run_command="srun -n 1 compute-sanitizer "
#run_command_tools="srun -n 1 "
run_command_tools=" "

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

bin="./vlasiator_gpu_wid4_10ad0ccc52_tp"

diffbin="./vlsvdiff_DP_dt"
#folder for all reference data 
reference_dir="/scratch/project_2004522/testpackage"
#compare agains which revision
reference_revision="current"

echo "Using executable $bin and comparing against $reference_revision in directory $reference_dir"
# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
