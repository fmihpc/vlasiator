## Define test and runs
#vlasiator binary

bin=vlasiator
if [ ! -f $bin ]
then
   echo Executable $bin does not exist
   exit
fi
#where the tests are run
run_dir="run"

# choose tests to run
#run_tests=(1 2 )
run_tests=( 1 2 )

# test 1
test_name[1]="Fluctuations"
test_cfg[1]="data/Fluctuations_small.cfg"
comparison_vlsv[1]="grid.0000010.vlsv"

# test 2
test_name[2]="Magnetosphere"
test_cfg[2]="data/Magnetosphere_small.cfg"
comparison_vlsv[2]="bulk.0000005.vlsv"
# Solar wind data (named as sw1.dat) file location for magnetosphere test
sw_data="data/sw1.dat"

# define here the variables you want to be tested
variables_name=( "rho" "rho_v" "rho_v" "rho_v" "B" "B" "B" "E" "E" "E" )
# and the corresponding components to variables, 
variables_components=( 0 0 1 2 0 1 2 0 1 2 )
#arrays variables_name and variables_components should have same number of elements, e.g., 4th variable is variables_name[4] variables_components[4]=  

