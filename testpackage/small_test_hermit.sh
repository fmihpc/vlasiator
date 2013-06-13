#!/bin/bash -l
#PBS -l mppwidth=64
#PBS -l mppnppn=32
#PBS -l walltime=00:15:00
#PBS -V  
#PBS -N test_small

#threads
t=8   

#command for running stuff
run_command="aprun -n 8 -N 4 -d $t"


#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#folder for reference data 
reference_dir="/zhome/academic/HLRS/pri/iprsalft/vlasiator_reference_data/revision1956"

#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test small/medium/long
source small_test_definitions.sh
wait
#run the test
source run_tests.sh
