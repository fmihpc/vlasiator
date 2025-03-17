#!/bin/bash
#SBATCH -t 02:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=hile-g-tp
#SBATCH -C g
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -G 1
#SBATCH -c 16                 # CPU cores per task
#SBATCH -n 1                  # number of tasks
#SBATCH --mem-per-cpu=2G
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

# folder for all reference data
reference_dir="/wrk-kappa/group/spacephysics/vlasiator/testpackage"
#cd $SLURM_SUBMIT_DIR
#cd $reference_dir # don't run on /proj

bin="vlasiator"
diffbin="vlsvdiff_DP_hile_gpu"

#compare agains which revision
reference_revision="current"

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
# Turn off forced managed memory paging
export HSA_XNACK=0
# use extra threads for MPI in background
#export MPICH_ASYNC_PROGRESS=1
# allow more in-parallel queues (should be 2x threads)
export GPU_MAX_HW_QUEUES=32
# Load custom memory pool
export LD_PRELOAD=/wrk-kappa/users/markusb/vlasiator-mempool/libpreload-me.so

# This is important for multi-node performance on Carrington, but not required with Hile.
# export UCX_NET_DEVICES=eth0

# Would allow oversubscription of cores with hyperthreading, do not use.
# export OMP_WAIT_POLICY=PASSIVE

module load papi
module load cray-pmi
module load craype-accel-amd-gfx90a
module load rocm/6.2.0
module load libfabric/1.22.0
module list

# threads per job (equal to -c )
t=$SLURM_CPUS_PER_TASK
tasks=$SLURM_NTASKS
export OMP_NUM_THREADS=$t

#command for running stuff
run_command="srun --mpi=pmi2 -n $tasks "
small_run_command="srun --mpi=pmi2 -n 1 -N 1 "
run_command_tools="srun --mpi=pmi2 -n 1 "

# # Placement debugging commands
# lscpu | grep NUMA
# echo
# rocm-smi --showtoponuma
# echo
# srun -n 1 --mpi=pmi2 /appl/bin/hostinfo
# echo
# srun --mpi=pmi2 bash -c ' \
#           echo -n "task $SLURM_PROCID (node $SLURM_NODEID): "; \
#           taskset -cp $$' | sort
# echo
# srun --mpi=pmi2 -c $t -n 1 /appl/cray/experimental/xthi/xthi_mpi_mp
# echo
# export CRAY_OMP_CHECK_AFFINITY=TRUE

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
wait 

