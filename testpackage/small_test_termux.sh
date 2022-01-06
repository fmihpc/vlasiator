#!/bin/bash
run_command="sudo"
small_run_command="sudo"
run_command_tools="sudo"
t=4       #threads per process

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

bin="$HOME/vlasiator/vlasiator"
diffbin="$HOME/vlasiator/vlsvdiff_DP"

#folder for all reference data 
reference_dir="$HOME/vlasiator/testpackage/testpackage_data"
#compare agains which revision
reference_revision="current"

base_dir=`pwd`

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
