#!/bin/bash -l

run_command="mpirun -n 1"
small_run_command="mpirun -n 1"
run_command_tools="mpirun -n 1"

t=8
export OMP_NUM_THREADS=$t
export OMPI_MCA_io="^ompio"

bin="../vlasiator"
diffbin="../vlsvdiff_DP"

corepattern="/scratch/cores/core.%e.%p"
lines="$(cat /proc/sys/kernel/core_pattern)"
if [ "$corepattern" != "$lines" ]
then
   echo $corepattern | sudo tee /proc/sys/kernel/core_pattern
fi

ulimit -c unlimited
umask 007


#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0


#folder for all reference data 
reference_dir="/scratch/testpackage/"
#compare agains which revision. This can be a proper version string, or "current", which should be a symlink to the
#proper most recent one
reference_revision="current"



# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
