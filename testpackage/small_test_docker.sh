
t=1                   #threads per process

#No idea how many cores we have available on travis. On my laptop I have 4.
cores_per_node=4
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"
#command for running stuff
run_command="mpirun -n $tasks  --allow-run-as-root"
small_run_command="mpirun -n 1 --allow-run-as-root"
run_command_tools="mpirun -n 1 --allow-run-as-root"

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

#Define where binaries are
bin="../vlasiator"
diffbin="../vlsvdiff_DP"

#folder for all reference data 
reference_dir="/home/vlasiator/testpackage/"
#compare agains which revision. This can be a proper version string, or "current", which should be a symlink to the
#proper most recent one
reference_revision="current"

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
