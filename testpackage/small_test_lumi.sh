#!/bin/bash
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=ctestpackage
#SBATCH -A project_465000287
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH -c 16                 # CPU cores per task
#SBATCH -n 16                  # number of tasks
#SBATCH --mem=0
#SBATCH --hint=multithread


#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

# folder for all reference data 
reference_dir="/scratch/project_465000287/kempfyan/testpackage/dev/"
cd $SLURM_SUBMIT_DIR
#cd $reference_dir # don't run on /proj

bin="/scratch/project_465000287/kempfyan/testpackage/field_tracing/vlasiator"
diffbin="/scratch/project_465000287/kempfyan/testpackage/field_tracing/vlsvdiff_DP"

#compare agains which revision
reference_revision="d5a69b36fea0a44f62d6e89944926b092e784777__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"
#reference_revision="current"

# threads per job (equal to -c )
t=16
module --force purge
module load LUMI/22.08
module load cpeGNU
module load papi

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Carrington has 2 x 16 cores
cores_per_node=128
# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

#command for running stuff
#run_command="mpirun -n $tasks -N $nodes "
run_command="srun "
small_run_command="srun -n 1 "
run_command_tools="srun -n 1"

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
wait 20

