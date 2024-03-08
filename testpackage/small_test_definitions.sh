
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

# Counter for creating tests
index=1

#######
# ACCELERATION TESTS (1..5)
#######

# basic multipeak acceleration test (fixed timestep)
test_name[${index}]="acctest_1_maxw_500k_100k_20kms_10deg"
comparison_vlsv[${index}]="fullf.0000000.vlsv fullf.0000001.vlsv fullf.0000002.vlsv fullf.0000020.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# basic multipeak acceleration test (dynamic timestep)
test_name[${index}]="acctest_2_maxw_500k_100k_20kms_10deg"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# acceleration test w/ substepping
test_name[${index}]="acctest_3_substeps"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# Helium acceleration
test_name[${index}]="acctest_4_helium"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="helium/vg_rho helium/vg_v helium/vg_v helium/vg_v"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# Gyration test with protons and antiprotons (multipop)
test_name[${index}]="acctest_5_proton_antiproton"
comparison_vlsv[${index}]="fullf.0000000.vlsv fullf.0000001.vlsv fullf.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

#######
# 1D/2D TRANSLATION TESTS (6..9)
#######

# Flowthrough tests
test_name[${index}]="Flowthrough_trans_periodic"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv bulk.0000003.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

test_name[${index}]="Flowthrough_x_inflow_y_outflow"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv bulk.0000003.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

test_name[${index}]="Flowthrough_x_inflow_y_outflow_acc"
comparison_vlsv[${index}]="bulk.0000002.vlsv bulk.0000003.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# Flowthrough with timevarying inflow
test_name[${index}]="Flowthrough_1D_timevarying"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2 0 0"
single_cell[${index}]=1
((index+=1))

#######
# 3D AMR TRANSLATION TESTS (10..13)
#######

# Three simple signal flowthroughs
test_name[${index}]="transtest_1_amr_triangle"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000010.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v"
variable_components[${index}]="0 0 1 2"
((index+=1))

test_name[${index}]="transtest_2_amr_sinewave"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000010.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v"
variable_components[${index}]="0 0 1 2"
((index+=1))

test_name[${index}]="transtest_3_amr_square"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000010.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v"
variable_components[${index}]="0 0 1 2"
((index+=1))

# Large flowthrough test
test_name[${index}]="Flowthrough_amr"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# TBA: Large flowthrough test with dAMR

#######
# RESTARTING TESTS (14..17)
#######

# Restart tests. Writing and reading
test_name[${index}]="restart_write"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

test_name[${index}]="restart_read"
comparison_vlsv[${index}]="initial-grid.0000000.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# Restart tests. Writing and reading, with acceleration and fieldsolver
test_name[${index}]="restart_write_acc"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

test_name[${index}]="restart_read_acc"
comparison_vlsv[${index}]="initial-grid.0000000.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

#######
# FIELDSOLVER TESTS (18..19)
#######

# Field solver test
test_name[${index}]="test_fp_fsolver_only_3D"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 1 2 0 1 2"
((index+=1))

# Field solver test w/ subcycles
test_name[${index}]="test_fp_substeps"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 1 2 0 1 2"
((index+=1))

#######
# GLOBAL 2D TESTS (20..21)
#######

#Very small ecliptic magnetosphere, no subcycling in ACC or FS
test_name[${index}]="Magnetosphere_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

#Very small polar magnetosphere, with subcycling in ACC or FS
test_name[${index}]="Magnetosphere_polar_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton/vg_v_nonthermal proton/vg_v_nonthermal proton/vg_v_nonthermal proton/vg_ptensor_nonthermal_diagonal proton/vg_ptensor_nonthermal_diagonal proton/vg_ptensor_nonthermal_diagonal proton"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2 0 1 2 0 1 2"
((index+=1))

#######
# GLOBAL 3D TESTS (22..23)
#######

# Magnetosphere 3D
test_name[${index}]="Magnetosphere_3D_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# Ionosphere 3D
test_name[${index}]="Ionosphere_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e ig_upmappedarea ig_fac ig_rhon ig_potential"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2 0 0 0 0"
((index+=1))


# choose tests to run (default: all tests)
run_tests=( )
for (( c=1; c<$index; c++ ))
do
    run_tests+=($c)
done

# Alternatively, set tests manually, e.g.
# run_tests=( 1 6 9 )
