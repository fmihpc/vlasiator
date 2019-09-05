
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
comparison_phiprof[1]="phiprof_0.txt"
variable_names[1]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[1]="0 0 1 2"
single_cell[1]=1

# acceleration test w/ substepping
test_name[2]="acctest_3_substeps"
comparison_vlsv[2]="fullf.0000001.vlsv"
comparison_phiprof[2]="phiprof_0.txt"
variable_names[2]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[2]="0 0 1 2"
single_cell[2]=1

# translation test
test_name[3]="transtest_2_maxw_500k_100k_20kms_20x20"
comparison_vlsv[3]="fullf.0000001.vlsv"
comparison_phiprof[3]="phiprof_0.txt"
variable_names[3]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[3]="0 0 1 2"

test_name[4]="acctest_4_helium"
comparison_vlsv[4]="fullf.0000001.vlsv"
comparison_phiprof[4]="phiprof_0.txt"
variable_names[4]="helium/vg_rho helium/vg_v helium/vg_v helium/vg_v"
variable_components[4]="0 0 1 2"
single_cell[4]=1

# Gyration test with protons and antiprotons
test_name[5]="acctest_5_proton_antiproton"
comparison_vlsv[5]="fullf.0000001.vlsv"
comparison_phiprof[5]="phiprof_0.txt"
variable_names[5]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[5]="0 0 1 2"
single_cell[5]=1

# Restart tests. Writing and reading
test_name[6]="restart_write"
comparison_vlsv[6]="bulk.0000001.vlsv"
comparison_phiprof[6]="phiprof_0.txt"
variable_names[6]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[6]="0 0 1 2 0 1 2 0 1 2"
test_name[7]="restart_read"
comparison_vlsv[7]="initial-grid.0000000.vlsv"
comparison_phiprof[7]="phiprof_0.txt"
variable_names[7]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[7]="0 0 1 2 0 1 2 0 1 2"

#Very small ecliptic magnetosphere, no subcycling in ACC or FS
test_name[8]="Magnetosphere_small"
comparison_vlsv[8]="bulk.0000001.vlsv"
comparison_phiprof[8]="phiprof_0.txt"
variable_names[8]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton"
variable_components[8]="0 0 1 2 0 1 2 0 1 2"

#Very small polar magnetosphere, with subcycling in ACC or FS
test_name[9]="Magnetosphere_polar_small"
comparison_vlsv[9]="bulk.0000001.vlsv"
comparison_phiprof[9]="phiprof_0.txt"
variable_names[9]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton/vg_v_nonthermal proton/vg_v_nonthermal proton/vg_v_nonthermal proton/vg_ptensor_nonthermal_diagonal proton/vg_ptensor_nonthermal_diagonal proton/vg_ptensor_nonthermal_diagonal proton"
variable_components[9]="0 0 1 2 0 1 2 0 1 2 0 1 2 0 1 2"

# Field solver test
test_name[10]="test_fp_fsolver_only_3D"
comparison_vlsv[10]="fullf.0000001.vlsv"
comparison_phiprof[10]="phiprof_0.txt"
variable_names[10]="fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[10]="0 1 2 0 1 2"

# Field solver test w/ subcycles
test_name[11]="test_fp_substeps"
comparison_vlsv[11]="fullf.0000001.vlsv"
comparison_phiprof[11]="phiprof_0.txt"
variable_names[11]="fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[11]="0 1 2 0 1 2"

# Flowthrough tests
test_name[12]="Flowthrough_trans_periodic"
comparison_vlsv[12]="bulk.0000001.vlsv"
comparison_phiprof[12]="phiprof_0.txt"
variable_names[12]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[12]="0 0 1 2 0 1 2 0 1 2"

test_name[13]="Flowthrough_x_inflow_y_outflow"
comparison_vlsv[13]="bulk.0000001.vlsv"
comparison_phiprof[13]="phiprof_0.txt"
variable_names[13]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[13]="0 0 1 2 0 1 2 0 1 2"

test_name[14]="Flowthrough_x_inflow_y_outflow_acc"
comparison_vlsv[14]="bulk.0000001.vlsv"
comparison_phiprof[14]="phiprof_0.txt"
variable_names[14]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[14]="0 0 1 2 0 1 2 0 1 2"

# Self-consistent wave generation test
test_name[15]="Selfgen_Waves_Periodic"
comparison_vlsv[15]="fullf.0000001.vlsv"
comparison_phiprof[15]="phiprof_0.txt"
variable_names[15]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton"
variable_components[15]="0 0 1 2 0 1 2 0 1 2"

##AMR tests
# translation test
test_name[16]="transtest_amr"
comparison_vlsv[3]="fullf.0000001.vlsv"
comparison_phiprof[3]="phiprof_0.txt"
variable_names[3]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[3]="0 0 1 2"

# Flowthrough test
test_name[17]="Flowthrough_amr"
comparison_vlsv[12]="bulk.0000001.vlsv"
comparison_phiprof[12]="phiprof_0.txt"
variable_names[12]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[12]="0 0 1 2 0 1 2 0 1 2"
