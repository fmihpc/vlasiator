
## Define test and runs

# This will overwrite small test definitions of used simulataneously, so please use separately

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
# SUBGRID TESTS (1..3)
#######

# 1  Pitch-angle diffusion test with 64 spatial cells
test_name[${index}]="subgrid_2_diffusion_spatial"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_ptensor_diagonal proton/vg_ptensor_diagonal proton/vg_ptensor_diagonal proton/vg_ptensor_offdiagonal proton/vg_ptensor_offdiagonal proton/vg_ptensor_offdiagonal proton"
variable_components[${index}]="0 0 1 2 0 1 2"
((index+=1))

# 2  Pitch-angle diffusion test with large dt
test_name[${index}]="subgrid_3_diffusion_time"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_ptensor_diagonal proton/vg_ptensor_diagonal proton/vg_ptensor_diagonal proton/vg_ptensor_offdiagonal proton/vg_ptensor_offdiagonal proton/vg_ptensor_offdiagonal proton"
variable_components[${index}]="0 0 1 2 0 1 2"
single_cell[${index}]=1
((index+=1))

# 3  Pitch-angle diffusion test with mass conservation
test_name[${index}]="subgrid_4_diffusion_mass"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_ptensor_diagonal proton/vg_ptensor_diagonal proton/vg_ptensor_diagonal proton/vg_ptensor_offdiagonal proton/vg_ptensor_offdiagonal proton/vg_ptensor_offdiagonal proton"
variable_components[${index}]="0 0 1 2 0 1 2"
single_cell[${index}]=1
((index+=1))

# choose tests to run (default: all tests)
run_tests=( )
for (( c=1; c<$index; c++ ))
do
    run_tests+=($c)
done


# choose tests to run (default: all tests)
run_tests=( )
for (( c=1; c<$index; c++ ))
do
    run_tests+=($c)
done

# Alternatively, set tests manually, e.g.
#run_tests=( 3 )