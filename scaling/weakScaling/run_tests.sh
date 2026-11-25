##--------------------------------------------------##
## code, no need to touch ##
##--------------------------------------------------##

# define tab interval sequence so that we have aligned output
# this is now covering at least up proton/vg_ptensor_nonthermal_offdiagonal_0 and numbers printed at setprecision(3) with negative mantissa and exponent
tabseq="1,46,62,78,94,110"
tabs $tabseq &> /dev/null # suppress special character output, list matches expand below


## add absolute paths to folder names, filenames
reference_dir=$( readlink -f $reference_dir )
run_dir=$( readlink -f $run_dir )_$( date +%Y.%m.%d_%H.%M.%S)
reference_revision_parsed=$( readlink -f $reference_dir/$reference_revision )
bin=$( readlink -f $bin )
diffbin=$( readlink -f $diffbin )
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
    #automatically obtained from the --version output (this overwrites the reverence_revision variable from a launch script)
    reference_revision=${revision}${solveropts}
    echo "Computing reference results into ${reference_dir}/${reference_revision}"
    if [[ -v GITHUB_ENV && -z $GITHUB_ENV ]]
    then
       echo "REFERENCE_REVISION=${reference_dir}/${reference_revision}" >> "$GITHUB_ENV"
    fi
else
    echo "----------"
    echo "This will be verifying ${revision}_$solveropts against $reference_revision_parsed"
    echo "----------"
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
    echo -e "\n"
    echo "----------"
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
        reference_result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}
        if [ -e  $reference_result_dir ]
        then
            echo "Removing previous reference results"
            rm -rf $reference_result_dir
        fi

        mkdir -p $reference_result_dir
        cp * $reference_result_dir
    fi

    cd $base_dir



### TESTS #####
    if [ ! $create_verification_files == 1 ]
    then
##Compare test case with right solutions
        echo "--------------------------------------------------------------------------------------------"
        echo "${test_name[$run]}  -  Verifying ${revision}_$solveropts against $reference_revision"
        echo "--------------------------------------------------------------------------------------------"
        reference_result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}

     #print header


        echo "------------------------------------------------------------"
        echo " ref-time     |   new-time       |  speedup                |"
        echo "------------------------------------------------------------"
	if [ -e  ${reference_result_dir}/${comparison_phiprof[$run]} ]
	then
            refPerf=$(grep "Propagate   " ${reference_result_dir}/${comparison_phiprof[$run]} |gawk  '(NR==1){print $11}')
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

        tabs 1,14,33,59 &> /dev/null # match next line
        echo  -e " $refPerf\t|  $newPerf\t|  $speedup\t|" | expand -t 1,14,33,59 # match previous line
        echo "------------------------------------------------------------"
        tabs $tabseq &> /dev/null # reset for other printouts
        echo -e " variable\t| absolute diff\t| relative diff |" | expand -t $tabseq # list matches tabs above
        echo "------------------------------------------------------------"

	variables=(${variable_names[$run]// / })
	indices=(${variable_components[$run]// / })
        for vlsv in ${comparison_vlsv[$run]}
        do
            if [ ! -f "${vlsv_dir}/${vlsv}" ]; then
                echo "Output file ${vlsv_dir}/${vlsv} not found!"
                echo "--------------------------------------------------------------------------------------------"
                continue
            fi
            if [ ! -f "${reference_result_dir}/${vlsv}" ]; then
                echo "Reference file ${reference_result_dir}/${vlsv} not found!"
                echo "--------------------------------------------------------------------------------------------"
                continue
            fi
            echo "Comparing file ${vlsv_dir}/${vlsv} against reference"
            for i in ${!variables[*]}
            do
                if [[ "${variables[$i]}" == "fg_"* ]]
                then
                    A=$( $run_command_tools $diffbin --meshname=fsgrid  ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} ${variables[$i]} ${indices[$i]} )
                    relativeValue=$(grep "The relative 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
                    absoluteValue=$(grep "The absolute 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
                    #print the results
                    echo -e " ${variables[$i]}_${indices[$i]}\t  ${absoluteValue}\t  ${relativeValue}" | expand -t $tabseq #list matches tabs above
                elif [[ "${variables[$i]}" == "ig_"* ]]
                then
                    B=$( $run_command_tools $diffbin --meshname=ionosphere  ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} ${variables[$i]} ${indices[$i]} )
                    relativeValue=$(grep "The relative 0-distance between both datasets" <<< $B |gawk '{print $8}'  )
                    absoluteValue=$(grep "The absolute 0-distance between both datasets" <<< $B |gawk '{print $8}'  )
                    #print the results
                    echo -e " ${variables[$i]}_${indices[$i]}\t  ${absoluteValue}\t  ${relativeValue}" | expand -t $tabseq # list matches tabs above
                elif [ ! "${variables[$i]}" == "proton" ]
                then # Regular vg_ variable
                    C=$( $run_command_tools $diffbin ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} ${variables[$i]} ${indices[$i]} )
                    relativeValue=$(grep "The relative 0-distance between both datasets" <<< $C |gawk '{print $8}'  )
                    absoluteValue=$(grep "The absolute 0-distance between both datasets" <<< $C |gawk '{print $8}'  )
                    #print the results
                    echo -e " ${variables[$i]}_${indices[$i]}\t  ${absoluteValue}\t  ${relativeValue}" | expand -t $tabseq # list matches tabs above
                elif [ "${variables[$i]}" == "proton" ]
                then
                    echo "--------------------------------------------------------------------------------------------"
                    echo "   Distribution function diff                                                               "
                    echo "--------------------------------------------------------------------------------------------"
                    $run_command_tools $diffbin ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} proton 0
                fi
            done # loop over variables

            # Print also time difference, if it is not zero
            timeDiff=$(grep "delta t" <<< $C |gawk '{print $8}'  )
            if (( $(awk 'BEGIN{print ('$timeDiff'!= 0.0)?1:0}') ))
            then
                echo "WARNING! VLSV file timestamps differ by ${timeDiff}s."
            else
                echo "VLSV file timestamps match."
            fi
            echo "--------------------------------------------------------------------------------------------"
        done # loop over vlsv files to compare
        echo "--------------------------------------------------------------------------------------------"
    fi
done # loop over tests
