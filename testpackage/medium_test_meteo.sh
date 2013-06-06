#!/bin/bash -l
#PBS -l mppwidth=120
#PBS -l mppnppn=12
#PBS -l walltime=01:00:00
#PBS -V  
#PBS -N test

#threads
t=12  

#command for running stuff
run_command="aprun -n 10 -N 1 -d $t"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#folder for reference data 
reference_dir="/stornext/field/users/hoilijo/Vlasiator/reference_data"


#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test
source /stornext/field/users/hoilijo/Vlasiator/vlasiator/trunk/testpackage/medium_test_definitions.sh
wait
# Run test
source /stornext/field/users/hoilijo/Vlasiator/vlasiator/trunk/testpackage/run_tests.sh
