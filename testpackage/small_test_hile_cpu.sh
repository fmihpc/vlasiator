#!/bin/bash
#SBATCH -t 02:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=hile_c_tp
#SBATCH -C c
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -c 8                 # CPU cores per task
#SBATCH -n 16                  # number of tasks
#SBATCH --mem=0
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
diffbin="vlsvdiff_DP_hile_cpu"

#compare agains which revision
reference_revision="current"

# This is important for multi-node performance on Carrington, but not required with Hile.
# export UCX_NET_DEVICES=eth0

# Would allow oversubscription of cores with hyperthreading, do not use.
# export OMP_WAIT_POLICY=PASSIVE

module load papi
module load cray-pmi
module load libfabric/1.22.0
#module load gdb4hpc
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
# srun -n 1 --mpi=pmi2 /appl/bin/hostinfo
# echo
# srun --mpi=pmi2 bash -c ' \
#           echo -n "task $SLURM_PROCID (node $SLURM_NODEID): "; \
#           taskset -cp $$' | sort
# echo
# srun --mpi=pmi2 -c $SLURM_CPUS_PER_TASK -n $tasks /appl/cray/experimental/xthi/xthi_mpi_mp
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

