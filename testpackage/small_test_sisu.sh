#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -J small_test
#SBATCH -p test
#SBATCH -N 2
#SBATCH --no-requeue          

ht=2                   #hyper threads per physical core
t=6                   #threads per process

#Compute and set  stuff, do not change
if [ -z $SLURM_NNODES ]
then
    #if running interactively we use 2 nodes
    nodes=2
else
    nodes=$SLURM_NNODES
fi

#sisu has 2 x 12 cores
cores_per_node=24
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

export MPICH_GNI_MDD_SHARING=disabled 
#export MPICH_GNI_MAX_EAGER_MSG_SIZE=0
#export MPICH_GNI_LOCAL_CQ_SIZE=4096


umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"
#command for running stuff
run_command="aprun -n $tasks -N $tasks_per_node -d $OMP_NUM_THREADS -j $ht"
run_command_tools="aprun -n 1"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0


#folder for all reference data 
reference_dir="/proj/vlasiato/testpackage/"
#compare agains which revision. This can be a proper version string, or "current", which should be a symlink to the
#proper most recent one
reference_revision="current"



# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
