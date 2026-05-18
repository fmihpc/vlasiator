#!/bin/bash

# use this when running locally/on a frontend
mpirun_cmd=" "
# or pick the launcher you need
#mpirun_cmd="srun"
#mpirun_cmd="mpirun"
#mpirun_cmd="aprun"

vlasiator=$1
cfg=$2

# Return a nonzero exit code if anything invalid was found
retval=0

if [ $# -ne 2 ]
then
    cat    <<EOF
Prints out differences between parameters in a cfg file and the options that the vlasiator executable understands.

Usage: $0 vlasiator_executable cfg_file

EOF
    exit 127
fi


if [ ! -x $vlasiator ]
then
    echo "ERROR: Vlasiator executable $vlasiator does not exist or is not executable"
    exit 127
fi


if [ ! -e $cfg ]
then
    echo "ERROR: cfg file $cfg does not exist"
    exit 127
fi


# Extract variables from Vlasiator help output.
$mpirun_cmd $vlasiator --check_cfg --run_config="$cfg" 
exit $?
