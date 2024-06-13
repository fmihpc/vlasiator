#!/bin/bash
# (SBATCH headers are removed from this file, as this is dispatched directly from the github actions script)

#If 1, the reference vlsv files are generated
# if 0 then we check the v1
create_verification_files=0

# folder for all reference data
reference_dir="/proj/group/spacephysics/vlasiator_testpackage/"
cd $SLURM_SUBMIT_DIR
#cd $reference_dir # don't run on /proj

bin="$GITHUB_WORKSPACE/vlasiator"
diffbin="$GITHUB_WORKSPACE/vlsvdiff_DP"

#compare agains which revision
reference_revision="CI_reference"

# threads per job (equal to -c )
t=4
module purge
#module load gnu9/9.3.0
#module load openmpi4/4.0.5
#module load pmix/3.1.4

module purge
module load GCC/11.2.0
module load OpenMPI/4.1.1-GCC-11.2.0
module load PMIx/4.1.0-GCCcore-11.2.0
module load PAPI/6.0.0.1-GCCcore-11.2.0

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Carrington has 2 x 16 cores
cores_per_node=32
# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

# With this the code won't print the warning, so we have a shorter report
export OMPI_MCA_io="^ompio"
export UCX_NET_DEVICES="eth0"

#command for running stuff
run_command="srun --job-name 'CI_run' -M carrington --nodes=1 -n 16 -c 4 --mem=0 -p short -t 0:30:00 --mpi=pmix -n $tasks"
small_run_command="srun --job-name 'CI_smallrun' -M carrington --nodes=1 -n 1 -c 4 --mem=0 -p short -t 0:30:00 --mpi=pmix -n 1"
run_command_tools="srun --job-name 'CI_tool' -M carrington --nodes=1 -n 1 -c 4 --mem=0 -p short -t 0:30:00 --mpi=pmix -n 1"

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Print the used node
hostname

# Define test
source test_definitions_small.sh
wait


if [ $create_verification_files == 1 ]; then
   echo "ERROR: creating reference data with the CI script is not supported."
   exit 1
fi

# define tab interval sequence so that we have aligned output
# this is now covering at least up proton/vg_ptensor_nonthermal_offdiagonal_0 and numbers printed at setprecision(3) with negative mantissa and exponent
tabseq="1,46,62,78,94,110"
tabs $tabseq &> /dev/null # suppress special character output, list matches expand below

# Note we are *not* using run_tests.sh here, as we are creating JUnit XML output.

# Get absolute paths
reference_dir=$( readlink -f $reference_dir )
reference_revision_full=$( readlink $reference_dir/$reference_revision )
run_dir=$( readlink -f $run_dir )_$( date +%Y.%m.%d_%H.%M.%S )
bin=$( readlink -f $bin )
diffbin=$( readlink -f $diffbin )
test_dir=$( readlink -f $test_dir)

flags=$(  $run_command $bin  --version |grep CXXFLAGS)
solveropts=$(echo $flags|sed 's/[-+]//g' | gawk '{for(i = 1;i<=NF;i++) { if( $i=="DDP" || $i=="DFP" || index($i,"PF")|| index($i,"DVEC") || index($i,"SEMILAG") ) printf "__%s", $(i) }}')
revision=$( $run_command $bin --version |gawk '{if(flag==1) {print $1;flag=0}if ($3=="log") flag=1;}' )

echo "----------"
echo "This will be verifying ${revision}_$solveropts against $reference_revision"
echo "----------"

#$small_run_command $bin --version > VERSION.txt 2> $GITHUB_WORKSPACE/stderr.txt

echo -e "### Testpackage output:\n" >> $GITHUB_STEP_SUMMARY
echo "CI_reference pointed to $reference_revision_full" >> $GITHUB_STEP_SUMMARY

NONZEROTESTS=0
ZEROTESTS=0
FAILEDTESTS=0

# loop over different test cases
for run in ${run_tests[*]}; do
   # directory for test results
   vlsv_dir=${run_dir}/${test_name[$run]}
   cfg_dir=${test_dir}/${test_name[$run]}

   # Check if folder for new run exists, if not create them, otherwise delete old results
   if [ ! -d ${vlsv_dir} ]; then
      mkdir -p ${vlsv_dir}
   else
      rm -f ${vlsv_dir}/*
   fi

   cd ${vlsv_dir}
   cp ${cfg_dir}/* .

   export OMP_NUM_THREADS=$t

   # Run prerequisite script, if it exists
   test -e test_prelude.sh && ./test_prelude.sh

   # Run the actual simulation
   if [[ ${single_cell[$run]} ]]; then
      { $small_run_command $bin --run_config=${test_name[$run]}.cfg 2>&1 1>&3 3>&- | tee $GITHUB_WORKSPACE/stderr.txt; exit ${PIPESTATUS[0]}; } 3>&1 1>&2 | tee $GITHUB_WORKSPACE/stdout.txt
   else
      { $run_command $bin --run_config=${test_name[$run]}.cfg 2>&1 1>&3 3>&- | tee $GITHUB_WORKSPACE/stderr.txt; exit ${PIPESTATUS[0]}; } 3>&1 1>&2 | tee $GITHUB_WORKSPACE/stdout.txt
   fi

   # Store error return value
   RUN_ERROR=${PIPESTATUS[0]}

   if [[ $RUN_ERROR != 0 ]]; then
      echo -e "<details><summary>:red_circle: ${test_name[$run]}: Failed to run or died with an error.</summary>\n"  >> $GITHUB_STEP_SUMMARY
      echo -e "Stdout:\n \`\`\`\n" >> $GITHUB_STEP_SUMMARY
      cat $GITHUB_WORKSPACE/stdout.txt >> $GITHUB_STEP_SUMMARY
      echo -e "\`\`\`\nStderr:\n \`\`\`\n" >> $GITHUB_STEP_SUMMARY
      cat $GITHUB_WORKSPACE/stderr.txt >> $GITHUB_STEP_SUMMARY
      echo -e "\`\`\`\n</details>" >> $GITHUB_STEP_SUMMARY
      FAILEDTESTS=$((FAILEDTESTS+1))

      touch $GITHUB_WORKSPACE/testpackage_failed
      continue
   fi

    # Run postprocessing script, if it exists
   test -e test_postproc.sh && ./test_postproc.sh

   ##Compare test case with right solutions
   { {
   echo -e "\n"
   echo "----------"
   echo "running ${test_name[$run]} "
   echo "----------"
   } 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt;} 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt
   reference_result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}

   { {
   echo "--------------------------------"
   echo " ref-time | new-time | speedup |"
   echo "--------------------------------"
   } 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt;} 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt
   if [ -e  ${reference_result_dir}/${comparison_phiprof[$run]} ]; then
      refPerf=$(grep "Propagate   " ${reference_result_dir}/${comparison_phiprof[$run]} | gawk '(NR==1){print $11}')
   else
      refPerf="NA"
   fi
   if [ -e ${vlsv_dir}/${comparison_phiprof[$run]} ]; then
      newPerf=$(grep "Propagate   " ${vlsv_dir}/${comparison_phiprof[$run]} | gawk '(NR==1){print $11}')
   else
      newPerf="NA"
   fi

   #print speedup if both refPerf and newPerf are numerical values
   speedup=$( echo $refPerf $newPerf |gawk '{if($2 == $2 + 0 && $1 == $1 + 0 ) print $1/$2; else print "NA"}')
   { {
   tabs 1,12,23 &> /dev/null # match next line
   echo  -e " $refPerf\t$newPerf\t$speedup" | expand -t 1,12,23 # match previous line
   echo "-------------------------------------------"
   tabs $tabseq &> /dev/null # reset for other printouts
   echo -e " variable\t| absolute diff\t| relative diff |" | expand -t $tabseq # list matches tabs above
   echo "-------------------------------------------"
   } 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt;} 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt

   {
   MAXERR=0.  # Absolute error
   MAXREL=0.  # Relative error
   MAXDT=0.   # Output time difference
   MAXERRVAR=""  # Variable with max absolute error
   MAXRELVAR=""  # Variable with max relative error
   COMPAREDFILES=0 # How many files we successfully compare
   TOCOMPAREFILES=0 # How many files we were supposed to compare
   # Save initial values in case we call continue in the loops below
   echo $COMPAREDFILES > $RUNNER_TEMP/COMPAREDFILES.txt
   echo $TOCOMPAREFILES > $RUNNER_TEMP/TOCOMPAREFILES.txt
   variables=(${variable_names[$run]// / })
   indices=(${variable_components[$run]// / })

   for vlsv in ${comparison_vlsv[$run]}
   do
       TOCOMPAREFILES=$((TOCOMPAREFILES+1))
       echo $TOCOMPAREFILES > $RUNNER_TEMP/TOCOMPAREFILES.txt
       if [ ! -f "${vlsv_dir}/${vlsv}" ]; then
           echo "Output file ${vlsv_dir}/${vlsv} not found!"
           echo "----------"
           continue
       fi
       if [ ! -f "${reference_result_dir}/${vlsv}" ]; then
            echo "Reference file ${reference_result_dir}/${vlsv} not found!"
           echo "----------"
           continue
       fi
       echo "Comparing file ${vlsv_dir}/${vlsv} against reference"
       COMPAREDFILES=$((COMPAREDFILES+1))
       echo $COMPAREDFILES > $RUNNER_TEMP/COMPAREDFILES.txt
       
       for i in ${!variables[*]}
       do
           if [[ "${variables[$i]}" == "fg_"* ]]
           then
               A=$( $run_command_tools $diffbin --meshname=fsgrid  ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} ${variables[$i]} ${indices[$i]} )
               relativeValue=$(grep "The relative 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
               absoluteValue=$(grep "The absolute 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
               #print the results
               echo -e " ${variables[$i]}_${indices[$i]}\t  ${absoluteValue}\t  ${relativeValue}" | expand -t $tabseq #list matches tabs above

               # Also log to metrics file
               echo "test_carrington{test=\"${test_name[$run]}\",var=\"${variables[$i]}\",index=\"${indices[$i]}\",diff=\"absolute\"} $absoluteValue" >> $GITHUB_WORKSPACE/metrics.txt
               echo "test_carrington{test=\"${test_name[$run]}\",var=\"${variables[$i]}\",index=\"${indices[$i]}\",diff=\"relative\"} $relativeValue" >> $GITHUB_WORKSPACE/metrics.txt

               # Check if we have a new maximum error
               if (( $( echo "$absoluteValue $MAXERR" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
                   MAXERR=$absoluteValue
                   MAXERRVAR=${variables[$i]}
               fi
               # ... or new max relative error
               if (( $( echo "$relativeValue $MAXREL" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
                   MAXREL=$relativeValue
                   MAXRELVAR=${variables[$i]}
               fi

           elif [[ "${variables[$i]}" == "ig_"* ]]
           then
               A=$( $run_command_tools $diffbin --meshname=ionosphere  ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} ${variables[$i]} ${indices[$i]} )
               relativeValue=$(grep "The relative 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
               absoluteValue=$(grep "The absolute 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
               # print the results
               echo -e " ${variables[$i]}_${indices[$i]}\t  ${absoluteValue}\t  ${relativeValue}" | expand -t $tabseq # list matches tabs above

               # Also log to metrics file
               echo "test_carrington{test=\"${test_name[$run]}\",var=\"${variables[$i]}\",index=\"${indices[$i]}\",diff=\"absolute\"} $absoluteValue" >> $GITHUB_WORKSPACE/metrics.txt
               echo "test_carrington{test=\"${test_name[$run]}\",var=\"${variables[$i]}\",index=\"${indices[$i]}\",diff=\"relative\"} $relativeValue" >> $GITHUB_WORKSPACE/metrics.txt

               # Check if we have a new maximum error
               if (( $( echo "$absoluteValue $MAXERR" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
                   MAXERR=$absoluteValue
                   MAXERRVAR=${variables[$i]}
               fi
               # ... or new max relative error
               if (( $( echo "$relativeValue $MAXREL" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
                   MAXREL=$relativeValue
                   MAXRELVAR=${variables[$i]}
               fi

           elif [ ! "${variables[$i]}" == "proton" ]
           then # Regular vg_ variable
               A=$( $run_command_tools $diffbin ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} ${variables[$i]} ${indices[$i]} )
               relativeValue=$(grep "The relative 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
               absoluteValue=$(grep "The absolute 0-distance between both datasets" <<< $A |gawk '{print $8}'  )
               #print the results
               echo -e " ${variables[$i]}_${indices[$i]}\t  ${absoluteValue}\t  ${relativeValue}" | expand -t $tabseq # list matches tabs above

               # Also log to metrics file
               echo "test_carrington{test=\"${test_name[$run]}\",var=\"${variables[$i]}\",index=\"${indices[$i]}\",diff=\"absolute\"} $absoluteValue" >> $GITHUB_WORKSPACE/metrics.txt
               echo "test_carrington{test=\"${test_name[$run]}\",var=\"${variables[$i]}\",index=\"${indices[$i]}\",diff=\"relative\"} $relativeValue" >> $GITHUB_WORKSPACE/metrics.txt

               # Check if we have a new maximum error
               if (( $( echo "$absoluteValue $MAXERR" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
                   MAXERR=$absoluteValue
                   MAXERRVAR=${variables[$i]}
               fi
               # ... or new max relative error
               if (( $( echo "$relativeValue $MAXREL" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
                   MAXREL=$relativeValue
                   MAXRELVAR=${variables[$i]}
               fi

           elif [ "${variables[$i]}" == "proton" ]
           then
               echo "----------"
               echo "Distribution function diff"
               echo "----------"
               $run_command_tools $diffbin ${reference_result_dir}/${vlsv} ${vlsv_dir}/${vlsv} proton 0
           fi

       done # loop over variables

       # Check if dt is nonzero
       timeDiff=$(grep "delta t" <<< $A |gawk '{print $8}'  )
       if (( $(awk 'BEGIN{print ('$timeDiff'!= 0.0)?1:0}') )); then
           if (( $( echo "${timeDiff#-} $MAXDT" | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
               MAXDT=$timeDiff
           fi
           echo "WARNING! VLSV file timestamps differ by ${timeDiff}s."
       else
           echo "VLSV file timestamps match."
       fi
       echo "----------"

       # This loop runs in a subshell (because of the stdout and stderr capture below),
       # so we save the variables to temp files
       echo $MAXERR > $RUNNER_TEMP/MAXERR.txt
       echo $MAXREL > $RUNNER_TEMP/MAXREL.txt
       echo $MAXDT > $RUNNER_TEMP/MAXDT.txt
       echo $MAXERRVAR > $RUNNER_TEMP/MAXERRVAR.txt
       echo $MAXRELVAR > $RUNNER_TEMP/MAXRELVAR.txt
       echo $speedup > $RUNNER_TEMP/speedup.txt
       echo $COMPAREDFILES > $RUNNER_TEMP/COMPAREDFILES.txt
       echo $TOCOMPAREFILES > $RUNNER_TEMP/TOCOMPAREFILES.txt

   done 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt; } 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt
   # end loop over vlsvfiles

   # Recover error variables
   COMPAREDFILES=`cat $RUNNER_TEMP/COMPAREDFILES.txt`
   TOCOMPAREFILES=`cat $RUNNER_TEMP/TOCOMPAREFILES.txt`
   if [[ $COMPAREDFILES -eq 0 ]]; then
       MAXERR=-42
       MAXERRVAR="n/a"
       MAXREL=42
       MAXRELVAR="n/a"
       MAXDT=0
       speedup=0
   else
       MAXERR=`cat $RUNNER_TEMP/MAXERR.txt`
       MAXERRVAR=`cat $RUNNER_TEMP/MAXERRVAR.txt`
       MAXREL=`cat $RUNNER_TEMP/MAXREL.txt`
       MAXRELVAR=`cat $RUNNER_TEMP/MAXRELVAR.txt`
       MAXDT=`cat $RUNNER_TEMP/MAXDT.txt`
       speedup=`cat $RUNNER_TEMP/speedup.txt`
   fi

   # Output CI step annotation
   if [[ $COMPAREDFILES -ne $TOCOMPAREFILES ]]; then
      echo -e "<details><summary>:red_square: ${test_name[$run]}: Comparison failure, accessed \`$COMPAREDFILES\` out of \`$TOCOMPAREFILES\` files: \`$MAXERRVAR\` has absolute error $MAXERR, \`$MAXRELVAR\` has relative error $MAXREL. Max timestamp difference is $MAXDT.   Speedup: $speedup</summary>\n" >> $GITHUB_STEP_SUMMARY
      NONZEROTESTS=$((NONZEROTESTS+1))
   elif (( $( echo "$MAXERR 0." | awk '{ if($1 > $2) print 1; else print 0 }' ) )) || (( $( echo "$MAXDT 0." | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
      echo -e "<details><summary>:large_orange_diamond: ${test_name[$run]}: Nonzero diffs: \`$MAXERRVAR\` has absolute error $MAXERR, \`$MAXRELVAR\` has relative error $MAXREL. Max timestamp difference is $MAXDT.   Speedup: $speedup</summary>\n" >> $GITHUB_STEP_SUMMARY
      NONZEROTESTS=$((NONZEROTESTS+1))
   else
      echo -e "<details><summary>:heavy_check_mark: ${test_name[$run]}: Ran with zero diffs. Speedup: $speedup</summary>\n" >> $GITHUB_STEP_SUMMARY
      ZEROTESTS=$((ZEROTESTS+1))
   fi

   echo -e "Stdout:\n \`\`\`\n" >> $GITHUB_STEP_SUMMARY
   cat $GITHUB_WORKSPACE/stdout.txt >> $GITHUB_STEP_SUMMARY
   echo -e "\`\`\`\nStderr:\n \`\`\`\n" >> $GITHUB_STEP_SUMMARY
   cat $GITHUB_WORKSPACE/stderr.txt >> $GITHUB_STEP_SUMMARY
   echo -e "\`\`\`\n</details>" >> $GITHUB_STEP_SUMMARY

   # -------- Generate JUnit output ----------
   JUNIT_FILE=$GITHUB_WORKSPACE/testpackage_output_junit_${test_name[$run]}.xml
   cat > $JUNIT_FILE <<HEADER
<?xml version="1.0" encoding="UTF-8"?>
<testsuites>
HEADER
if [[ $RUN_ERROR == 0 ]]; then
   cat >> $JUNIT_FILE <<HEADER
   <testsuite name="Testpackage on carrington" tests="1" errors="0">
      <testcase classname="vlasiator_carrington.${test_name[$run]}" name="${test_name[$run]}" time="$newPerf">
HEADER
else
   cat >> $JUNIT_FILE <<HEADER
   <testsuite name="Vlasiator testpackage" tests="1" errors="1">
      <testcase classname="vlasiator_carrington.${test_name[$run]}" name="${test_name[$run]}" time="$newPerf">
         <error message="Running test failed">The test failed to run. Stdout was:

HEADER
   cat $GITHUB_WORKSPACE/stdout.txt >> $JUNIT_FILE
   cat >> $JUNIT_FILE <<HEADER

Stderr was:

HEADER
   cat $GITHUB_WORKSPACE/stderr.txt >> $JUNIT_FILE
   cat >> $JUNIT_FILE <<HEADER

         </error>
HEADER
fi
echo -e '         <system-out><![CDATA[' >> $JUNIT_FILE
cat $GITHUB_WORKSPACE/stdout.txt >> $JUNIT_FILE
echo -e ']]>\n         </system-out>\n         <system-err><![CDATA[' >> $JUNIT_FILE
cat $GITHUB_WORKSPACE/stderr.txt >> $JUNIT_FILE
cat >> $JUNIT_FILE <<FOOTER
]]>
         </system-err>
      </testcase>
   </testsuite>
</testsuites>
FOOTER
done


# -- Write summary for github PR annotation --
echo "summary=$ZEROTESTS tests with zero diffs, $NONZEROTESTS tests with diffs, $FAILEDTESTS tests failed." >> $GITHUB_WORKSPACE/testpackage_output_variables.txt
if [[ $FAILEDTESTS > 0 ]]; then
   echo "conclusion=failure" >> $GITHUB_WORKSPACE/testpackage_output_variables.txt
elif [[ $NONZEROTESTS > 0 ]]; then
   echo "conclusion=neutral" >> $GITHUB_WORKSPACE/testpackage_output_variables.txt
else
   echo "conclusion=success" >> $GITHUB_WORKSPACE/testpackage_output_variables.txt
fi
