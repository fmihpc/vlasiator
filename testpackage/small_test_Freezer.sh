#!/bin/bash -l

run_command="mpirun -n 1"
small_run_command="mpirun -n 1"
run_command_tools="mpirun -n 1"

cd ..

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0


#folder for all reference data 
reference_dir="/scratch/testpackage/dev/"
#compare agains which revision. This can be a proper version string, or "current", which should be a symlink to the
#proper most recent one
reference_revision="current"



# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
