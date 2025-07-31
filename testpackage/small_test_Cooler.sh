#!/bin/bash -l

run_command="mpirun -n 1"
small_run_command="mpirun -n 1"
run_command_tools="mpirun -n 1"

#cd ..

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0


#folder for all reference data 
reference_dir="/home/humikael/testpackage/dev/"
#compare agains which revision. This can be a proper version string, or "current", which should be a symlink to the
#proper most recent one
reference_revision="9e172e4b2937f9faca1afb0fabfc85a5408fb08a__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"
#reference_revision="b77539d572fdde7d0b49b7ac129935f544814e94__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC_FALLBACK_GENERIC__DVECL=16"

bin=../vlasiator
diffbin=../vlsvdiff_DP

t=1

export OMP_NUM_THREADS=$t

# Define test
source test_definitions_small.sh
#source additional_test_definitions.sh
wait
# Run tests
source run_tests.sh
