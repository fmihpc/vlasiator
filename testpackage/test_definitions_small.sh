
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

# 1  basic multipeak acceleration test (fixed timestep)
test_name[${index}]="acctest_1_maxw_500k_100k_20kms_10deg"
comparison_vlsv[${index}]="fullf.0000000.vlsv fullf.0000001.vlsv fullf.0000002.vlsv fullf.0000020.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# 2  basic multipeak acceleration test (dynamic timestep)
test_name[${index}]="acctest_2_maxw_500k_100k_20kms_10deg"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# 3  acceleration test w/ substepping
test_name[${index}]="acctest_3_substeps"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v proton"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# 4  Helium acceleration
test_name[${index}]="acctest_4_helium"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="helium/vg_rho helium/vg_v helium/vg_v helium/vg_v"
variable_components[${index}]="0 0 1 2"
single_cell[${index}]=1
((index+=1))

# 5  Gyration test with protons and antiprotons (multipop)
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

# 6  Flowthrough tests (no boundaries, includes empty cells)
test_name[${index}]="Flowthrough_trans_periodic"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv bulk.0000003.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 7  Inflow and outflow boundaries
test_name[${index}]="Flowthrough_x_inflow_y_outflow"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv bulk.0000003.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 8  Inflow and outflow boundaries together with acceleration
test_name[${index}]="Flowthrough_x_inflow_y_outflow_acc"
comparison_vlsv[${index}]="bulk.0000002.vlsv bulk.0000003.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 9  Flowthrough with timevarying inflow
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

# 10 AMR translation, single v-cell, triangle signal
test_name[${index}]="transtest_1_amr_triangle"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000010.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v"
variable_components[${index}]="0 0 1 2"
((index+=1))

# 11 AMR translation, single v-cell, sinewave signal
test_name[${index}]="transtest_2_amr_sinewave"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000010.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v"
variable_components[${index}]="0 0 1 2"
((index+=1))

# 12 AMR translation, single v-cell, square signal
test_name[${index}]="transtest_3_amr_square"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000010.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v"
variable_components[${index}]="0 0 1 2"
((index+=1))

# 13 Large AMR translation flowthrough test
test_name[${index}]="Flowthrough_amr"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 14 Large (dynamic) AMR translation flowthrough test
test_name[${index}]="Flowthrough_damr"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv bulk.0000003.vlsv bulk.0000004.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v vg_amr_alpha1 vg_amr_alpha2 fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 0 0 1 2 0 1 2"
((index+=1))

# TBA: Large flowthrough test with dAMR

#######
# RESTARTING TESTS (15..18)
#######

# Restart tests. Writing and reading
# 14 Restart write with translation only
test_name[${index}]="restart_write"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 15 Restart read and propagate with translation only
test_name[${index}]="restart_read"
comparison_vlsv[${index}]="initial-grid.0000000.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 16 Restart write with translation, acceleration, and fieldsolver
test_name[${index}]="restart_write_acc"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 17 Restart read and propagate with translation, acceleration, and fieldsolver
test_name[${index}]="restart_read_acc"
comparison_vlsv[${index}]="initial-grid.0000000.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

#######
# FIELDSOLVER TESTS (19..20)
#######

# 18 3D Field solver test
test_name[${index}]="test_fp_fsolver_only_3D"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 1 2 0 1 2"
((index+=1))

# 19 3D Field solver test w/ subcycles
test_name[${index}]="test_fp_substeps"
comparison_vlsv[${index}]="fullf.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 1 2 0 1 2"
((index+=1))

#######
# GLOBAL 2D TESTS (21..22)
#######

# 20 Very small ecliptic magnetosphere, no subcycling in ACC or FS
test_name[${index}]="Magnetosphere_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 21 Very small polar magnetosphere, with subcycling in ACC or FS
test_name[${index}]="Magnetosphere_polar_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e proton/vg_v_nonthermal proton/vg_v_nonthermal proton/vg_v_nonthermal proton/vg_ptensor_nonthermal_diagonal proton/vg_ptensor_nonthermal_diagonal proton/vg_ptensor_nonthermal_diagonal proton"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2 0 1 2 0 1 2"
((index+=1))

#######
# GLOBAL 3D TESTS (23..24)
#######

# 22 Magnetosphere 3D, very small, 40 timesteps
test_name[${index}]="Magnetosphere_3D_small"
comparison_vlsv[${index}]="bulk.0000001.vlsv bulk.0000002.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"
((index+=1))

# 23 Ionosphere 3D (not a very physical or successful test at the moment but verifies some things about IG grid and outputs all datareducers)
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
