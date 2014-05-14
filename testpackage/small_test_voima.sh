#!/bin/bash -l
#PBS -l mppwidth=40
#PBS -l mppnppn=20
#PBS -l walltime=00:15:00
#PBS -V  
#PBS -N test_small


#threads
t=10 

#command for running stuff
run_command="aprun -n 4 -d $t"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#folder for reference data 
reference_dir="/stornext/field/vlasiator/test_package_data/master_14_05_2014/"

#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
