##--------------------------------------------------##
## code, no need to touch ##
##--------------------------------------------------##


## add absolute paths to folder names, filenames
reference_dir=$( readlink -f $reference_dir )
run_dir=$( readlink -f $run_dir )
bin=$( readlink -f $bin )
test_dir=$( readlink -f $test_dir)

# for run in ${run_tests[*]}
#   do
#   test_dir[$run]=$( readlink -f runs/${test_name[$run]} )
#   test_cfg[$run]=$( readlink -f runs/${test_name[$run]}/${test_name[$run]}.cfg )
# done  


if [ $create_verification_files == 1 ]
then
    #if we create the references, then lets simply run in the reference dir and turn off tests below. Revision is 
    #automatically obtained from the --version output
    reference_revision=$( $run_command $bin --version |gawk '{if(flag==1) {print $1;flag=0}if ($3=="log") flag=1;}' )
    flags=$(  $run_command $bin  --version |grep CXXFLAGS)
    solveropts=$(echo $flags|sed 's/[-+]//g' | gawk '{for(i = 1;i<=NF;i++) { if( $i=="DDP" || $i=="DFP" || index($i,"PF")|| index($i,"PV") || index($i,"SEMILAG") ) printf "__%s", $(i) }}')
    reference_revision=${reference_revision}${solveropts}
    echo "Computing reference results into ${reference_dir}/${reference_revision}"

fi


#create main run folder if it doesn not exist
if [ !  -d $run_dir ]
then
    mkdir -p $run_dir 
fi


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

    $run_command $bin --run_config=${test_name[$run]}.cfg

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
        echo "${test_name[$run]}  -  Verifying $bin against revision $reference_revision"    
        echo "--------------------------------------------------------------------------------------------" 
        result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}

     #print header


        echo "------------------------------------------------------------"
        echo " ref-perf     |   new-perf       |  speedup                |"
        echo "------------------------------------------------------------"
        refPerf=$(grep "Propagate   " ${result_dir}/${comparison_phiprof[$run]}  |gawk '{print $12}')
        newPerf=$(grep "Propagate   " ${vlsv_dir}/${comparison_phiprof[$run]}  |gawk '{print $12}')
        speedup=$( echo $refPerf $newPerf |gawk '{print $2/$1}')
        echo  "$refPerf        $newPerf         $speedup"
        echo "------------------------------------------------------------"
        echo "  variable     |     absolute diff     |     relative diff | "
        echo "------------------------------------------------------------"
        for i in ${!variables_name[*]}
        do
            if [ ! "${variables_name[$i]}" == "avgs" ]
            then
                relativeValue=$( vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables_name[$i]} ${variables_components[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
                absoluteValue=$( vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables_name[$i]} ${variables_components[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
#print the results      
                echo "${variables_name[$i]}_${variables_components[$i]}                $absoluteValue                 $relativeValue    "
            fi

        done # loop over variables


        for i in ${!variables_name[*]}
        do
            if [ "${variables_name[$i]}" == "avgs" ]
            then
                echo "--------------------------------------------------------------------------------------------" 
                echo "   Distribution function diff                                                               "
                echo "--------------------------------------------------------------------------------------------" 
                vlsvdiff_DP ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} avgs 0
            fi 
        done # loop over variables

        echo "--------------------------------------------------------------------------------------------" 
    fi
done # loop over tests



