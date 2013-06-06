#!/bin/bash -l
#PBS -l mppwidth=192
#PBS -l mppnppn=32
#PBS -l walltime=01:00:00
#PBS -V  
#PBS -N test

#threads
t=8   

#command for running stuff
run_command="aprun -n 24 -N 4 -d $t"


#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#folder for reference data 
reference_dir="/univ_1/ws1/ws/iprshoil-test_testpackage-0/reference_data"


#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test small/medium/long
source /zhome/academic/HLRS/pri/iprshoil/vlasiator/trunk/testpackage/medium_test_definitions.sh
wait
#run the test
source /zhome/academic/HLRS/pri/iprshoil/vlasiator/trunk/testpackage/run_tests.sh
