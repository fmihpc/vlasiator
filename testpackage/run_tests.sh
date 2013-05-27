# source test_package_meteo.sh
# wait
# source test_definitions.sh

##--------------------------------------------------##
## code, no need to touch ##
##--------------------------------------------------##



## add absolute paths to folder names, filenames
reference_dir=$( readlink -f $reference_dir )
run_dir=$( readlink -f $run_dir )
bin=$( readlink -f $bin )


for run in ${run_tests[*]}
  do
  test_cfg[$run]=$( readlink -f ${test_cfg[$run]} )
  if [ ${test_name[$run]} = "Magnetosphere" ]; then
     sw_data=$( readlink -f $sw_data )
  fi
done  


if [ $create_verification_files == 1 ]
    then
    #if we create the references, then lets simply run in the reference dir and turn off tests below
    run_dir=$reference_dir
    echo "Computing reference results into $run_dir"
else
    #print header
    echo "Verifying $bin"
    echo "--------------------------------------------------------------------------------------------" 
    echo "      Test case            |     variable     |       value         |        ok/bad            "  
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
  
# change to run directory of the test case, e.g. test_Fluctuations
  cd ${vlsv_dir}
   if [ ${test_name[$run]} = "Magnetosphere" ]; then
      cp $sw_data .
   fi
#   
  export OMP_NUM_THREADS=$t
  export MPICH_MAX_THREAD_SAFETY=funneled
  $run_command -n $p -N 1 -d $t $bin --run_config=${test_cfg[$run]}
  cd $base_dir



### TESTS #####
  if [ ! $create_verification_files == 1 ]
      then
##Compare test case with right solutions
      result_dir=${reference_dir}/${test_name[$run]}

#print header
      
      echo "--------------------------------------------------------------------------------------------" 
      for i in ${!variables_name[*]}
        do

        relativeValue=$( vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables_name[$i]} ${variables_components[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
        absoluteValue=$( vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables_name[$i]} ${variables_components[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
        
        #if mean is zero, then relative value has no meaning, use abs instead
        value=$(awk -vrelVal="$relativeValue" -vabsVal="$absoluteValue" 'BEGIN{print (absVal==0 && relVal==-1)?absVal:relVal}')
        status=$(awk -vval="$value" 'BEGIN{print  (val<=10^(-14) && val>=-10^(-14))?"OK":"!!! ERROR !!!" }')

#print the results      
        echo "      ${test_name[$run]}          ${variables_name[$i]} ${variables_components[$i]}                $value                    $status              "
      done # loop over variables
      echo "--------------------------------------------------------------------------------------------" 
  fi
done # loop over tests

  

