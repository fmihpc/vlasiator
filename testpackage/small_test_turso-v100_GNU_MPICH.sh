#!/bin/bash
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --job-name=TP_turso-v100
#SBATCH --gres=gpu:V100:1
#SBATCH --cpus-per-task 10                 # CPU cores per task
#SBATCH --hint=nomultithread
#SBATCH --nodes=1
#SBATCH -n 1                  # number of tasks
#SBATCH --mem=40G

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

# folder for all reference data
reference_dir="/turso/group/spacephysics/vlasiator/testpackage/"
cd $SLURM_SUBMIT_DIR

#compare agains which revision
#reference_revision="CI_reference"
reference_revision="current"

bin="../vlasiator"
diffbin="../vlsvdiff_DP"

ulimit -c unlimited
source ../modules/turso-v100_GNU_MPICH.sh

nodes=$SLURM_NNODES
t=$SLURM_CPUS_PER_TASK # used by TP script
export OMP_NUM_THREADS=$t
export tasks=$SLURM_NTASKS

#command for running stuff
# MPICH here is PMIx-only (no Hydra), so no mpirun.
run_command="srun --mpi=pmix -c $t -n $tasks"
small_run_command="srun --mpi=pmix -c $t -n 1"
run_command_tools="srun --mpi=pmix -n 1"

# Informational / debugging placement outputs
#
# lscpu | grep NUMA
# echo
# nvidia-smi topo -m
# echo
# srun --mpi=pmix --ntasks-per-node=1 bash -c ' \
#           echo -n "task $SLURM_PROCID (node $SLURM_NODEID): "; \
#           taskset -cp $$' | sort
# echo

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
wait
