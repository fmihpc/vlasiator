#!/bin/bash
#SBATCH -t 00:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=repeat_test_acc
##SBATCH -A spacephysics
#SBATCH --constraint="amd"
# test short medium 20min1d 3d
#SBATCH -p short
#SBATCH --exclusive
##SBATCH --nodes=1
#SBATCH -c 1                 # CPU cores per task
#SBATCH -n 1                  # number of tasks
#SBATCH --array=0-9
#SBATCH --mem=10G
#SBATCH --no-requeue
##SBATCH --hint=multithread

# If 1, the reference vlsv files are generated
# if 0 then we check the v1 against reference files
create_verification_files=0

# folder for all reference data
reference_dir="/home/siclasse/"
mainfolder="/home/siclasse/vlasiator_hile/vlasiator/testpackage/"
cd $SLURM_SUBMIT_DIR

diffbin="../vlsvdiff_DP"

#compare agains which revision
#reference_revision="CI_reference"
reference_revision="rng_ref2"

module purge
module load GCC/13.2.0
module load OpenMPI/4.1.6-GCC-13.2.0
module load PMIx/4.2.6-GCCcore-13.2.0
module load PAPI/7.1.0-GCCcore-13.2.0
module load Boost/1.83.0-GCC-13.2.0
module load xthi
export UCX_TLS=dc_mlx5
export UCX_NET_DEVICES=mlx5_0:1

export OMPI_MCA_btl='^uct,ofi'
export OMPI_MCA_pml='ucx'
export OMPI_MCA_mtl='^ofi'
#Carrington has 2 x 16 cores per node, plus hyperthreading
ht=2
# t=$SLURM_CPUS_PER_TASK
# export OMP_NUM_THREADS=$t
#command for running stuff
run_command="mpirun -n 1 "
small_run_command="mpirun -n 1 "
run_command_tools="mpirun -np 1 "

umask 007
# Launch the OpenMP job to the allocated compute node
# echo "Running $exec on $SLURM_NTASKS mpi tasks, with $t threads per task on $SLURM_NNODES nodes ($ht threads per physical core)"

# Optional debug printouts
# srun -np 1 /appl/bin/hostinfo
# srun --cpu-bind=cores bash -c 'echo -n "task $SLURM_PROCID (node $SLURM_NODEID): "; taskset -cp $$' | sort
# srun --mpi=pmix --cpu-bind=cores xthi

# Define test
source test_definitions_small.sh
wait
# Run tests
jobcount=$(($SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1))
index=$(($SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN))

# echo $index $jobcount $SLURM_ARRAY_TASK_COUNT $SLURM_ARRAY_TASK_ID
# echo $index
for n in $(seq $(($index * 50)) $((($index + 1) * 50 - 1))); do
  echo "N=$n"
  bin="/home/siclasse/vlasiator_hile/vlasiator/vlasiator --map_order_shift $n"
  source run_tests.sh

  wait
  # mv logfile.txt logfile_$n.txt
  cd $mainfolder
  # cd $SLURM_SUBMIT_DIR
  wait
done
