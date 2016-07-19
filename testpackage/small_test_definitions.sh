
## Define test and runs
#vlasiator binary

bin=vlasiator
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
run_tests=( 1 2 3 4 5 6 7 8 9)

# acceleration test
test_name[1]="acctest_2_maxw_500k_100k_20kms_10deg"
comparison_vlsv[1]="fullf.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[1]="phiprof_full_1.txt"

# translation test
test_name[2]="transtest_2_maxw_500k_100k_20kms_20x20"
comparison_vlsv[2]="fullf.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[2]="phiprof_full_0.txt"

#very small magnetosphere, tests all at once
test_name[3]="Magnetosphere_small"
comparison_vlsv[3]="bulk.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[3]="phiprof_full_0.txt"

#Acceleration of 1 maxwellian, corresponds to SW
test_name[4]="acctest_1_maxw_500k_30kms_1deg"
comparison_vlsv[4]="fullf.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[4]="phiprof_full_1.txt"

# Field solver test
test_name[5]="test_fp_fsolver_only_3D"
comparison_vlsv[5]="fullf.0000001.vlsv"
comparison_phiprof[5]="phiprof_full_0.txt"

# Flowthrough tests
test_name[6]="Flowthrough_trans_periodic"
comparison_vlsv[6]="bulk.0000001.vlsv"
comparison_phiprof[6]="phiprof_full_0.txt"

test_name[7]="Flowthrough_x_inflow_y_outflow"
comparison_vlsv[7]="bulk.0000001.vlsv"
comparison_phiprof[7]="phiprof_full_0.txt"

test_name[8]="Flowthrough_x_inflow_y_outflow_acc"
comparison_vlsv[8]="bulk.0000001.vlsv"
comparison_phiprof[8]="phiprof_full_0.txt"

# Self-consistent wave generation test
test_name[9]="Selfgen_Waves_Periodic"
comparison_vlsv[9]="fullf.0000001.vlsv"
comparison_phiprof[9]="phiprof_full_0.txt"

# define here the variables you want to be tested
variables_name=( "rho" "rho_v" "rho_v" "rho_v" "B" "B" "B" "E" "E" "E" "proton" )
# and the corresponding components to variables, 
variables_components=( 0 0 1 2 0 1 2 0 1 2 0)
#arrays variables_name and variables_components should have same number of elements, e.g., 4th variable is variables_name[4] variables_components[4]=  

