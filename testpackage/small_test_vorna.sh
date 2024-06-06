#!/bin/bash
#SBATCH -t 00:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=testpackage
##SBATCH -A spacephysics 
#SBATCH -M vorna
# test short medium 20min1d 3d
#SBATCH -p short
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -c 8                  # CPU cores per task
#SBATCH -n 4                  # number of tasks (4xnodes)
#SBATCH --mem=55G # memory per node
#SBATCH -x vorna-[435-436]

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=1

# folder for all reference data 
reference_dir="/proj/group/spacephysics/vlasiator_testpackage/"
cd $SLURM_SUBMIT_DIR
#cd $reference_dir # don't run on /proj

bin="/proj/USERNAME/BINARYNAME"
diffbin="/proj/group/spacephysics/vlasiator_testpackage/vlsvdiff_DP_vorna"

#compare agains which revision
#reference_revision="c36241b84ce8179f7491ebf2a94c377d7279e8c9__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"
reference_revision="current"

# threads per job (equal to -c )
t=8
module purge
module load GCC/10.2.0
module load OpenMPI/4.0.5-GCC-10.2.0

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Vorna has 2 x 8 cores
cores_per_node=16
# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

#command for running stuff
run_command="srun --mpi=pmix_v2 -n $tasks -N $nodes "
small_run_command="srun --mpi=pmix_v2 -n 1 -N 1 "
run_command_tools="srun --mpi=pmix_v2 -n 1 -N 1 "

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
wait 20

