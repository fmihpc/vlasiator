#!/bin/bash -l
#SBATCH --time=0:15:00
#SBATCH --partition=gputest
##SBATCH --time=24:00:00
##SBATCH --partition=gpusmall
#SBATCH --cpus-per-task=32
#SBATCH --job-name=vlasigputest
#SBATCH --account=project_2004522
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores

module purge
#mahti_cuda
module load StdEnv
module load cuda/11.5.0
module load gcc/9.4.0
module load openmpi/4.1.2-cuda
module load jemalloc
module load papi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projappl/project_2000203/libraries_rhel8_gcccuda/boost/lib

ulimit -c unlimited

ht=1    #hyper threads per physical core
t=$OMP_NUM_THREADS     #threads per process
echo "threads ${t}"
nodes=$SLURM_NNODES

#mahti has 2 x 64 cores
cores_per_node=128
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
#tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
#tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
tasks=1
tasks_per_node=1

export OMPI_MCA_io=^ompio

# hint from George and https://github.com/openucx/ucx/issues/5504
export UCX_TLS=ud,ud_v
export OMPI_MCA_coll=^hcoll

umask 007

# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core) and $tasks_per_node tasks per node"
#command for running stuff (add compute-sanitizer here if necessary)
run_command="srun "
small_run_command="srun -n 1 "
run_command_tools="srun -n 1 "

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

bin="./vlasiator_gpu_tp"
diffbin="/projappl/project_2000203/vlsvdiff_DP"

#folder for all reference data 
reference_dir="/scratch/project_2000203/testpackage_data"
#compare agains which revision
#reference_revision="current"
reference_revision="5e656cf2e4e2dcc1fd4d521010fa655ef1b2fa09__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
