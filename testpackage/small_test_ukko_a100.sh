#!/bin/bash
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=TP_ukko_a100
#SBATCH -M ukko
#SBATCH -p gpu
##SBATCH -p gpu-oversub # Oversub affinities can be whatever
#SBATCH --constraint=a100
#SBATCH --gres=gpu:1
##SBATCH --cpus-per-gpu=8
#SBATCH --hint=nomultithread
#SBATCH --nodes=1
#SBATCH -c 8                 # CPU cores per task
#SBATCH -n 1                  # number of tasks
##SBATCH --mem=0 # do not request all node memory or it's equal to exclusive
#SBATCH --mem=60G
#SBATCH --exclude=ukko3-g602 # exclude oversubscription node

# Debugging: different placements
##SBATCH --distribution=block:cyclic
##SBATCH --distribution=block:block
#echo "dist block:cyclic"
#echo "dist default"
#echo "dist block:block"

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

# folder for all reference data
reference_dir="/proj/group/spacephysics/vlasiator_testpackage/"
cd $SLURM_SUBMIT_DIR

#compare agains which revision
reference_revision="CI_reference"
#reference_revision="current"

bin="$SLURM_SUBMIT_DIR/vlasiator"
diffbin="$SLURM_SUBMIT_DIR/vlsvdiff_DP"

module purge;  ml OpenMPI/4.1.6.withucx-GCC-13.2.0 PAPI/7.1.0-GCCcore-13.2.0 CUDA/12.6.0

nodes=$SLURM_NNODES
t=$SLURM_CPUS_PER_TASK # used by TP script
export OMP_NUM_THREADS=$t
export tasks=$SLURM_NTASKS

#command for running stuff
run_command="srun --mpi=pmix -c $t -n $tasks"
small_run_command="srun --mpi=pmix -c $t -n 1"
run_command_tools="mpirun -n 1 -N 1"

# Placement debugging commands

# lscpu | grep NUMA
# echo
# nvidia-smi topo -m
# echo
# srun --mpi=pmix /appl/bin/hostinfo
# echo
# srun --mpi=pmix --ntasks-per-node=1 bash -c ' \
#           echo -n "task $SLURM_PROCID (node $SLURM_NODEID): "; \
#           taskset -cp $$' | sort
# echo
# module load xthi
# srun --mpi=pmix -c $t -n 1 xthi
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

