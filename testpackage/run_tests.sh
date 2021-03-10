##--------------------------------------------------##
## code, no need to touch ##
##--------------------------------------------------##


## add absolute paths to folder names, filenames
reference_dir=$( readlink -f $reference_dir )
run_dir=$( readlink -f $run_dir )_$( date +%Y.%m.%d_%H.%M.%S)

bin=$( readlink -f $bin )
test_dir=$( readlink -f $test_dir)

# for run in ${run_tests[*]}
#   do
#   test_dir[$run]=$( readlink -f runs/${test_name[$run]} )
#   test_cfg[$run]=$( readlink -f runs/${test_name[$run]}/${test_name[$run]}.cfg )
# done

if [[ ! $small_run_command ]]; then
	echo "No small_run_command provided in machine config, please update it!"
	exit
fi 
 
flags=$(  $run_command $bin  --version |grep CXXFLAGS)
solveropts=$(echo $flags|sed 's/[-+]//g' | gawk '{for(i = 1;i<=NF;i++) { if( $i=="DDP" || $i=="DFP" || index($i,"PF")|| index($i,"DVEC") || index($i,"SEMILAG") ) printf "__%s", $(i) }}')
revision=$( $run_command $bin --version |gawk '{if(flag==1) {print $1;flag=0}if ($3=="log") flag=1;}' )

if [ $create_verification_files == 1 ]
then
    #if we create the references, then lets simply run in the reference dir and turn off tests below. Revision is 
    #automatically obtained from the --version output
    reference_revision=${revision}${solveropts}
    echo "Computing reference results into ${reference_dir}/${reference_revision}"
fi


if [ -d $run_dir ]
then
    echo $run_dir exists?
    exit
fi
mkdir -p $run_dir 

# loop over different test cases
for run in ${run_tests[*]}
do
    echo "running ${test_name[$run]} "
# directory for test results
    vlsv_dir=${run_dir}/${test_name[$run]}
    cfg_dir=${test_dir}/${test_name[$run]}
    
# Check if folder for new run exists, if not create them, otherwise delete old results
    if [ ! -d ${vlsv_dir} ]; then
        mkdir -p ${vlsv_dir}
    else
        rm -f ${vlsv_dir}/*
    fi
    
# change to run directory of the test case, e.g. test_Fluctuations
    cd ${vlsv_dir}
    cp ${cfg_dir}/* .
    
    export OMP_NUM_THREADS=$t
    export MPICH_MAX_THREAD_SAFETY=funneled

    # Run prerequisite script, if it exists
    test -e test_prelude.sh && ./test_prelude.sh

    # Run the actual simulation
    if [[ ${single_cell[$run]} ]]; then
    	$small_run_command $bin --version  > VERSION.txt
	$small_run_command $bin --run_config=${test_name[$run]}.cfg
    else
	$run_command $bin --version  > VERSION.txt
	$run_command $bin --run_config=${test_name[$run]}.cfg
    fi

    # Run postprocessing script, if it exists
    test -e test_postproc.sh && ./test_postproc.sh


  ###copy new reference data to correct folder
    if [ $create_verification_files == 1 ]
    then
        result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}
        if [ -e  $result_dir ]
        then
            echo "remove old results"
            rm -rf $result_dir
        fi

        mkdir -p $result_dir
        cp * $result_dir      
    fi

    cd $base_dir



### TESTS #####
    if [ ! $create_verification_files == 1 ]
    then
##Compare test case with right solutions
        echo "--------------------------------------------------------------------------------------------" 
        echo "${test_name[$run]}  -  Verifying ${revision}_$solveropts against $reference_revision"    
        echo "--------------------------------------------------------------------------------------------" 
        result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}

     #print header


        echo "------------------------------------------------------------"
        echo " ref-time     |   new-time       |  speedup                |"
        echo "------------------------------------------------------------"
	if [ -e  ${result_dir}/${comparison_phiprof[$run]} ] 
	then
            refPerf=$(grep "Propagate   " ${result_dir}/${comparison_phiprof[$run]} |gawk  '(NR==1){print $11}')
	else
	    refPerf="NA"
	fi
	if [ -e ${vlsv_dir}/${comparison_phiprof[$run]} ] 
	then
            newPerf=$(grep "Propagate   " ${vlsv_dir}/${comparison_phiprof[$run]}  |gawk  '(NR==1){print $11}')
	else
	    newPerf="NA"
	fi
	#print speedup if both refPerf and newPerf are numerical values
        speedup=$( echo $refPerf $newPerf |gawk '{if($2 == $2 + 0 && $1 == $1 + 0 ) print $1/$2; else print "NA"}')
        echo  "$refPerf        $newPerf         $speedup"
        echo "------------------------------------------------------------"
        echo "  variable     |     absolute diff     |     relative diff | "
        echo "------------------------------------------------------------"

	variables=(${variable_names[$run]// / })
	indices=(${variable_components[$run]// / })
        for i in ${!variables[*]}
        do
            if [ "${variables[$i]}" == "fg_e" ] || [ "${variables[$i]}" == "fg_b" ]
            then
                relativeValue=$($run_command_tools vlsvdiff_DP --meshname=fsgrid  ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
                absoluteValue=$($run_command_tools vlsvdiff_DP --meshname=fsgrid  ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
#print the results      
                echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
            
            elif [ ! "${variables[$i]}" == "proton" ]
            then
                relativeValue=$($run_command_tools vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
                absoluteValue=$($run_command_tools vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
#print the results      
                echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
            elif [ "${variables[$i]}" == "proton" ]
            then
                echo "--------------------------------------------------------------------------------------------" 
                echo "   Distribution function diff                                                               "
                echo "--------------------------------------------------------------------------------------------" 
                $run_command_tools vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} proton 0
            fi 
        done # loop over variables

        echo "--------------------------------------------------------------------------------------------" 
    fi
done # loop over tests



