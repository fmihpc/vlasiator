#!/bin/bash -l
#PBS -l mppwidth=2
#PBS -l mppdepth=12
#PBS -l walltime=00:15:00
#PBS -V  
#PBS -N testing


## variables ##

p=2
t=12
run_command="aprun"

#IFS=$',' # separate test names by , without spaces

# choose tests to run
run_tests=(1 2)

# test 1
test_name[1]="Fluctuations"
test_cfg[1]="Fluctuations.cfg"
result_dirs[1]="reference/Fluctuations"

# test 1
test_name[2]="test_fp"
test_cfg[2]="test_fp.cfg"
result_dirs[2]="reference/test_fp"

# test 3
#test_name[3]="test_trans"
#test_cfg[3]="test_trans.cfg"
#result_dirs[3]="reference/test_trans"

# test 4
#test_name[4]="verificationLarmor"
#test_cfg[4]="verificationLarmor.cfg"
#result_dirs[4]="reference/verificationLarmor"

# give here the variables you want to be tested
variables_name=( "rho" "rho_v" "rho_v" "rho_v" "B" "B" "B" "E" "E" "E" )
# and the corresponding components to variables, 
variables_components=( 0 0 1 2 0 1 2 0 1 2)
#arrays variables_name and variables_components should have same number of elements, e.g., 4th variable is variables_name[4] variables_components[4]=  

# choose the length of the test run: short/long
short_long="short"



## code, no need to touch ##
cd $PBS_O_WORKDIR 

#print header
echo "------------------------------------------------------------------------------------------------------------" 
echo "      Test case         |     variable     |          value           |           ok/bad                    "  


# loop over different test cases
for run in ${run_tests[*]}
do

# test name
test=${test_name[$run]}
# cfg file
cfg=${test_cfg[$run]}
# directory where the good results are for comparison
result_dir=${result_dirs[$run]}
# directory for test results
vlsv_Dir=test_${test}

# Check if files for new vlsv files exists, if not create them
if [ ! -d ${vlsv_Dir} ]; then
    echo creating ${vlsv_Dir}
    mkdir ${vlsv_Dir}
else
   rm -f ${vlsv_Dir}/*.vlsv
fi

# copy cfg file and change to run directory of the test case, e.g. test_Fluctuations
cp $cfg ${vlsv_Dir}/
cd ${vlsv_Dir}

export OMP_NUM_THREADS=$t
export MPICH_MAX_THREAD_SAFETY=funneled


$run_command -n $p -d $t ../vlasiator --run_config=../$cfg


cd ..

##Compare test case with right solutions

# loop over variables
for i in ${!variables_name[*]}
do

value=$( vlsvdiff_DP ${result_dir}/grid.0000001.vlsv ${vlsv_Dir}/grid.0000001.vlsv ${variables_name[$i]} ${variables_components[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )


a=$(awk -vx="$value" 'BEGIN{print (x > 0 ? x:-x)}')
b=$(awk -vn="$a" 'BEGIN{print ( n<=10^(-14))?1:0 }')


if [ "$b" == "1" ]; then
      status="ok"
else
      status="bad"
fi

#print the results      
echo "------------------------------------------------------------------------------------------------------------" 
echo "      $test           ${variables_name[$i]} ${variables_components[$i]}            $value                $status              "
echo "------------------------------------------------------------------------------------------------------------"
done # loop over variables
 
done # loop over tests
