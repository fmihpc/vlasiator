
## Define test and runs

# Assume binaries are set in job script
#bin=vlasiator
#diffbin=vlsvdiff_DP

if [ ! -f $bin ]
then
   echo Executable $bin does not exist
   exit
fi

# where the tests are run
run_dir="run"

# where the directories for different tests, including cfg and other needed data files are located 
test_dir="tests"

# choose tests to run
tests=()

# If no tests are specified, run all tests/ subdirectories.
if [ "${#tests[@]}" -eq "0" ]; then
   tests=($(ls -d tests/*))
   echo "No run_tests list specified, using all ${#tests[@]} subdirectories of tests/"
fi;


# Iterate through the specified tests and add them to the table
num_tests=0
for i in ${tests[@]}; do
   if [ -f $i/test_spec.sh ]; then
      test_name[$num_tests]=$i
      . $i/test_spec.sh
      num_tests=$((num_tests + 1))
   else 
      echo "Warning: test $i does not have a test_spec.sh."
   fi
done

run_tests=($(seq 0 $((num_tests-1))))

echo $num_tests test to run in total.

