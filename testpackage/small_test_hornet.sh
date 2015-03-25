#!/bin/bash
#PBS -N small_test
#PBS -l nodes=2
#PBS -l walltime=0:25:00
#PBS -W umask=007
#PBS -q test
nodes=2   #Has to be identical to above! $( qstat -f $PBS_JOBID | grep   Resource_List.nodes | gawk '{print $3}' )
ht=2      #hyper threads per physical core
t=6       #threads per process
exec="./vlasiator"

#Hornet has 2 x 12 cores
cores_per_node=24
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t


umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"
#command for running stuff
run_command="aprun -n $tasks -N $tasks_per_node -d $OMP_NUM_THREADS -j $ht"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0


#folder for all reference data 
reference_dir="/zhome/academic/HLRS/pri/iprsealf/test_package_references"
#compare agains which revision
reference_revision="voima_0a8be8ac087c2da56556e4f56fcb3e5826aa6f38__DVEC4D_AGNER__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF"


# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
