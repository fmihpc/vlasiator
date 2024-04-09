#!/bin/bash
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --job-name=CI_testpackage
##SBATCH -A spacephysics 
#SBATCH -M carrington
# test short medium 20min1d 3d
#SBATCH -p short
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -c 4                 # CPU cores per task
#SBATCH -n 16                  # number of tasks
#SBATCH --mem=0
##SBATCH -x carrington-[801-808]

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

#command for running stuff
run_command="mpirun --mca btl self -mca pml ^vader,tcp,openib,uct,yalla -x UCX_NET_DEVICES=mlx5_0:1 -x UCX_TLS=rc,sm -x UCX_IB_ADDR_TYPE=ib_global -np $tasks"
small_run_command="mpirun --mca btl self -mca pml ^vader,tcp,openib,uct,yalla -x UCX_NET_DEVICES=mlx5_0:1 -x UCX_TLS=rc,sm -x UCX_IB_ADDR_TYPE=ib_global -n 1 -N 1"
run_command_tools="mpirun -np 1 "

umask 007
# Launch the OpenMP job to the allocated compute node
echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"

# Print the used node
hostname

# Define test
source small_test_definitions.sh
wait


if [ $create_verification_files == 1 ]; then
   echo "ERROR: creating reference data with the CI script is not supported."
   exit 1
fi

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

#$small_run_command $bin --version > VERSION.txt 2> $GITHUB_WORKSPACE/stderr.txt

echo -e "### Testpackage output:\n" >> $GITHUB_STEP_SUMMARY
echo "CI_reference pointed to $reference_revision_full" >> $GITHUB_STEP_SUMMARY

NONZEROTESTS=0
ZEROTESTS=0
FAILEDTESTS=0

# loop over different test cases
for run in ${run_tests[*]}; do
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
   # Fore set to error if output file does not exist
   if [ ! -f ${vlsv_dir}/${comparison_vlsv[$run]} ]; then
       RUN_ERROR=1
   fi

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
   echo "--------------------------------------------------------------------------------------------"
   echo "${test_name[$run]}  -  Verifying ${revision}_$solveropts against $reference_revision"
   echo "--------------------------------------------------------------------------------------------"
   } 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt;} 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt
   result_dir=${reference_dir}/${reference_revision}/${test_name[$run]}

   { {
   echo "------------------------------------------------------------"
   echo " ref-time     |   new-time       |  speedup                |"
   echo "------------------------------------------------------------"
   } 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt;} 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt
   if [ -e  ${result_dir}/${comparison_phiprof[$run]} ]; then
      refPerf=$(grep "Propagate   " ${result_dir}/${comparison_phiprof[$run]} | gawk '(NR==1){print $11}')
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
   echo  "$refPerf        $newPerf         $speedup"
   echo "------------------------------------------------------------"
   echo "  variable     |     absolute diff     |     relative diff | "
   echo "------------------------------------------------------------"
   } 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt;} 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt

   {
   MAXERR=0.  # Absolute error
   MAXREL=0.  # Relative error
   MAXERRVAR=""  # Variable with max absolute error
   MAXRELVAR=""  # Variable with max relative error

   variables=(${variable_names[$run]// / })
   indices=(${variable_components[$run]// / })

   for i in ${!variables[*]}
   do
       if [[ "${variables[$i]}" == "fg_"* ]]
       then
           relativeValue=$($run_command_tools $diffbin --meshname=fsgrid  ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
           absoluteValue=$($run_command_tools $diffbin --meshname=fsgrid  ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
           #print the results
           echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
   
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
           relativeValue=$($run_command_tools $diffbin --meshname=ionosphere  ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
           absoluteValue=$($run_command_tools $diffbin --meshname=ionosphere  ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
           # print the results
           echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
   
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
       then
           relativeValue=$($run_command_tools $diffbin ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )
           absoluteValue=$($run_command_tools $diffbin ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} ${variables[$i]} ${indices[$i]} |grep "The absolute 0-distance between both datasets" |gawk '{print $8}'  )
           #print the results
           echo "${variables[$i]}_${indices[$i]}                $absoluteValue                 $relativeValue    "
   
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
           echo "--------------------------------------------------------------------------------------------"
           echo "   Distribution function diff                                                               "
           echo "--------------------------------------------------------------------------------------------"
           $run_command_tools $diffbin ${result_dir}/${comparison_vlsv[$run]} ${vlsv_dir}/${comparison_vlsv[$run]} proton 0
       fi

       # This loop runs in a subshell (because of the stdout and stderr capture below),
       # so we save the variables to temp files
       echo $MAXERR > $RUNNER_TEMP/MAXERR.txt
       echo $MAXREL > $RUNNER_TEMP/MAXREL.txt
       echo $MAXERRVAR > $RUNNER_TEMP/MAXERRVAR.txt
       echo $MAXRELVAR > $RUNNER_TEMP/MAXRELVAR.txt
       echo $speedup > $RUNNER_TEMP/speedup.txt
   done 2>&1 1>&3 3>&- | tee -a $GITHUB_WORKSPACE/stderr.txt; } 3>&1 1>&2 | tee -a $GITHUB_WORKSPACE/stdout.txt

   # Recover error variables
   MAXERR=`cat $RUNNER_TEMP/MAXERR.txt`
   MAXERRVAR=`cat $RUNNER_TEMP/MAXERRVAR.txt`
   MAXREL=`cat $RUNNER_TEMP/MAXREL.txt`
   MAXRELVAR=`cat $RUNNER_TEMP/MAXRELVAR.txt`
   speedup=`cat $RUNNER_TEMP/speedup.txt`

   # Output CI step annotation
   if (( $( echo "$MAXERR 0." | awk '{ if($1 > $2) print 1; else print 0 }' ) )); then
      echo -e "<details><summary>:large_orange_diamond: ${test_name[$run]}: Nonzero diffs: : \`$MAXERRVAR\` has absolute error $MAXERR, \`$MAXRELVAR\` has relative error $MAXREL. Speedup: $speedup</summary>\n" >> $GITHUB_STEP_SUMMARY
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

