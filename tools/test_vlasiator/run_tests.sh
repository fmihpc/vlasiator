#!/bin/bash -l
#PBS -l mppwidth=2
#PBS -l mppdepth=12
#PBS -l walltime=00:15:00
#PBS -V  
#PBS -N testing


## variables ##
#vlasiator binary
bin=vlasiator
#processes and threads
p=2
t=12
#where the tests are run
run_dir="run"

# choose tests to run
#run_tests=(1 2
run_tests=( 1 2 )

# test 1
test_name[1]="Fluctuations_A"
test_cfg[1]="data/Fluctuations_A.cfg"
comparison_vlsv[1]="grid.0000010.vlsv"

# test 2
test_name[2]="test_fp_A"
test_cfg[2]="data/test_fp_A.cfg"
comparison_vlsv[2]="grid.0000003.vlsv"

#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0
#folder for reference data
reference_dir="/stornext/field/users/alfthan/vlasiator_reference_data/r1569"
# give here t%he variables you want to be tested
variables_name=( "rho" "rho_v" "rho_v" "rho_v" "B" "B" "B" "E" "E" "E" )
# and the corresponding components to variables, 
variables_components=( 0 0 1 2 0 1 2 0 1 2)
#arrays variables_name and variables_components should have same number of elements, e.g., 4th variable is variables_name[4] variables_components[4]=  




##--------------------------------------------------##
## code, no need to touch ##
##--------------------------------------------------##


#command for running stuff, FIXME: should be a function or so that could easily be extended to mpirun etc
run_command="aprun"

#get baseddir from PBS_O_WORKDIR if set (batch job), otherwise go to current folder
#http://stackoverflow.com/questions/307503/whats-the-best-way-to-check-that-environment-variables-are-set-in-unix-shellscr
base_dir=${PBS_O_WORKDIR:=$(pwd)}
cd  $base_dir

## add absolute paths to folder names, filenames
reference_dir=$( readlink -f $reference_dir )
run_dir=$( readlink -f $run_dir )
bin=$( readlink -f $bin )
for run in ${run_tests[*]}
  do
  test_cfg[$run]=$( readlink -f ${test_cfg[$run]} )
done  


if [ $create_verification_files == 1 ]
    then
    #if we create the references, then lets simply run in the reference dir and turn off tests below
    run_dir=$reference_dir
    echo "Computing reference results into $run_dir"
else
    #print header
    echo "Verifying $bin"
    echo "------------------------------------------------------------------------------------------------------------" 
    echo "      Test case         |     variable     |          value           |           ok/bad                    "  
fi


#create main run folder if it doesn not exist
if [ !  -d $run_dir ]
    then
    mkdir $run_dir 
fi


# loop over different test cases
for run in ${run_tests[*]}
  do
# directory for test results
  vlsv_dir=${run_dir}/${test_name[$run]}
  
# Check if folder for new run exists, if not create them, otherwise delete old results
  if [ ! -d ${vlsv_dir} ]; then
      mkdir ${vlsv_dir}
  else
      rm -f ${vlsv_dir}/*
  fi

# copy cfg file and change to run directory of the test case, e.g. test_Fluctuations
  cd ${vlsv_dir}
  
  export OMP_NUM_THREADS=$t
  export MPICH_MAX_THREAD_SAFETY=funneled
  $run_command -n $p -d $t $bin --run_config=${test_cfg[$run]}
  cd $base_dir



### TESTS #####
  if [ ! $create_verification_files == 1 ]
      then
##Compare test case with right solutions
      result_dir=${reference_dir}/${test_name[$run]}

#print header
      
      echo "------------------------------------------------------------------------------------------------------------"
      for i in ${!variables_name[*]}
        do

        relativeValue=$( vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables_name[$i]} ${variables_components[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
        absoluteValue=$( vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables_name[$i]} ${variables_components[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
        
        #if mean is zero, then relative value has no meaning, use abs instead
        value=$(awk -vrelVal="$relativeValue" -vabsVal="$absoluteValue" 'BEGIN{print (absVal==0 && relVal==-1)?absVal:relVal}')
        status=$(awk -vval="$value" 'BEGIN{print  (val<=10^(-14) && val>=-10^(-14))?"OK":"!!! ERROR !!!" }')

#print the results      
        echo "      ${test_name[$run]}          ${variables_name[$i]} ${variables_components[$i]}            $value                $status              "
      done # loop over variables
      echo "------------------------------------------------------------------------------------------------------------"
  fi
done # loop over tests

  

