#!/bin/bash
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=ctestpackage
##SBATCH -A spacephysics 
#SBATCH -M carrington
# test short medium 20min1d 3d
#SBATCH -p short
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -c 4                 # CPU cores per task
#SBATCH -n 16                  # number of tasks
#SBATCH --mem-per-cpu=5G
#SBATCH --hint=multithread
##SBATCH -x carrington-[801-808]

# If 1, the reference vlsv files are generated
# if 0 then we check the v1 against reference files
create_verification_files=0

# folder for all reference data 
reference_dir="/proj/group/spacephysics/vlasiator_testpackage/"
cd $SLURM_SUBMIT_DIR

bin="/proj/USERNAME/BINARYNAME"
diffbin="/proj/group/spacephysics/vlasiator_testpackage/vlsvdiff_DP_carrington"

#compare agains which revision
#reference_revision="CI_reference"
reference_revision="current"

module purge
module load GCC/13.2.0
module load OpenMPI/4.1.6-GCC-13.2.0
module load PMIx/4.2.6-GCCcore-13.2.0
module load PAPI/7.1.0-GCCcore-13.2.0
#module load xthi
export UCX_NET_DEVICES=eth0 # This is important for multi-node performance!

#Carrington has 2 x 16 cores per node, plus hyperthreading
ht=2
t=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$t

#command for running stuff
run_command="srun --mpi=pmix_v3 -n $SLURM_NTASKS "
small_run_command="srun --mpi=pmix_v3 -n 1"
run_command_tools="mpirun -np 1 "

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $SLURM_NTASKS mpi tasks, with $t threads per task on $SLURM_NNODES nodes ($ht threads per physical core)"

# Optional debug printouts
# srun -np 1 /appl/bin/hostinfo
# srun --cpu-bind=cores bash -c 'echo -n "task $SLURM_PROCID (node $SLURM_NODEID): "; taskset -cp $$' | sort
# srun --mpi=pmix --cpu-bind=cores xthi

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
wait 

