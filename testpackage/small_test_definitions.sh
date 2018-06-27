
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
run_tests=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14)

# acceleration test
test_name[1]="acctest_2_maxw_500k_100k_20kms_10deg"
comparison_vlsv[1]="fullf.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[1]="phiprof_1.txt"

# acceleration test w/ substepping
test_name[2]="acctest_3_substeps"
comparison_vlsv[2]="fullf.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[2]="phiprof_1.txt"

# translation test
test_name[3]="transtest_2_maxw_500k_100k_20kms_20x20"
comparison_vlsv[3]="fullf.0000001.vlsv"
comparison_phiprof[3]="phiprof_0.txt"

test_name[4]="acctest_4_helium"
comparison_vlsv[4]="fullf.0000001.vlsv"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[4]="phiprof_1.txt"

# Gyration test with protons and antiprotons
test_name[5]="acctest_5_proton_antiproton"
#only one process does anything -> in _1 phiprof here
comparison_phiprof[5]="phiprof_1.txt"

# Restart tests. Writing and reading
test_name[6]="restart_write"
comparison_vlsv[6]="restart.vlsv"
comparison_phiprof[6]="phiprof_0.txt"
test_name[7]="restart_read"
comparison_vlsv[7]="restart.vlsv"
comparison_phiprof[7]="phiprof_0.txt"

#Very small ecliptic magnetosphere, no subcycling in ACC or FS
test_name[8]="Magnetosphere_small"
comparison_vlsv[8]="bulk.0000001.vlsv"
comparison_phiprof[8]="phiprof_0.txt"

#Very small polar magnetosphere, with subcycling in ACC or FS
test_name[9]="Magnetosphere_polar_small"
comparison_vlsv[9]="bulk.0000001.vlsv"
comparison_phiprof[9]="phiprof_0.txt"

# Field solver test
test_name[10]="test_fp_fsolver_only_3D"
comparison_vlsv[10]="fullf.0000001.vlsv"
comparison_phiprof[10]="phiprof_0.txt"

# Field solver test w/ subcycles
test_name[11]="test_fp_substeps"
comparison_vlsv[11]="fullf.0000001.vlsv"
comparison_phiprof[11]="phiprof_0.txt"

# Flowthrough tests
test_name[12]="Flowthrough_trans_periodic"
comparison_vlsv[12]="bulk.0000001.vlsv"
comparison_phiprof[12]="phiprof_0.txt"

test_name[13]="Flowthrough_x_inflow_y_outflow"
comparison_vlsv[13]="bulk.0000001.vlsv"
comparison_phiprof[13]="phiprof_0.txt"

test_name[14]="Flowthrough_x_inflow_y_outflow_acc"
comparison_vlsv[14]="bulk.0000001.vlsv"
comparison_phiprof[14]="phiprof_0.txt"

# Self-consistent wave generation test
test_name[15]="Selfgen_Waves_Periodic"
comparison_vlsv[15]="fullf.0000001.vlsv"
comparison_phiprof[15]="phiprof_0.txt"


# define here the variables you want to be tested
variables_name=( "proton/rho" "proton/V" "proton/V" "proton/V" "B" "B" "B" "E" "E" "E" "proton" )
# and the corresponding components to variables, 
variables_components=( 0 0 1 2 0 1 2 0 1 2 0)
#arrays variables_name and variables_components should have same number of elements, e.g., 4th variable is variables_name[4] variables_components[4]=  

