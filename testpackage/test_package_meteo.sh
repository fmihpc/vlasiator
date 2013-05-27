#!/bin/bash -l
#PBS -l mppwidth=120
#PBS -l mppnppn=12
#PBS -l walltime=01:00:00
#PBS -V  
#PBS -N test
#command for running stuff, FIXME: should be a function or so that could easily be extended to mpirun etc
run_command="aprun"

#processes and threads
p=10    # mppwidth = p*t
t=12   # mppnppn = t

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

#folder for reference data 
reference_dir="/stornext/field/users/hoilijo/Vlasiator/reference_data"


#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test small/medium/long
source /stornext/field/users/hoilijo/Vlasiator/vlasiator/trunk/testpackage/medium_test_definitions.sh
wait
source /stornext/field/users/hoilijo/Vlasiator/vlasiator/trunk/testpackage/run_tests.sh
