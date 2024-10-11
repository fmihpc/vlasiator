#!/bin/bash
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=ctestpackage
##SBATCH -A spacephysics 
#SBATCH -M carrington
# test short medium 20min1d 3d
#SBATCH -p short
##SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -c 4                 # CPU cores per task
#SBATCH -n 16                  # number of tasks
#SBATCH --mem-per-cpu=5G
##SBATCH -x carrington-[801-808]

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=1

# folder for all reference data 
reference_dir="/proj/group/spacephysics/vlasiator_testpackage/"
cd $SLURM_SUBMIT_DIR
#cd $reference_dir # don't run on /proj

bin="/proj/USERNAME/BINARYNAME"
diffbin="/proj/group/spacephysics/vlasiator_testpackage/vlsvdiff_DP_carrington"

#compare agains which revision
#reference_revision="c36241b84ce8179f7491ebf2a94c377d7279e8c9__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"
reference_revision="current"

# threads per job (equal to -c )
t=4
module purge
module load GCC/13.2.0
module load OpenMPI/4.1.6-GCC-13.2.0
module load PMIx/4.2.6-GCCcore-13.2.0
module load PAPI/7.1.0-GCCcore-13.2.0

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Carrington has 2 x 16 cores
cores_per_node=32
# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

#command for running stuff
run_command="mpirun --mca btl self -mca pml ^vader,tcp,openib,uct,yalla -x UCX_NET_DEVICES=mlx5_0:1 -x UCX_TLS=rc,sm -x UCX_IB_ADDR_TYPE=ib_global -np $tasks"
small_run_command="mpirun --mca btl self -mca pml ^vader,tcp,openib,uct,yalla -x UCX_NET_DEVICES=mlx5_0:1 -x UCX_TLS=rc,sm -x UCX_IB_ADDR_TYPE=ib_global -n 1 -N 1"
run_command_tools="srun --mpi=pmix_v3 -n 1 "

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
wait 

