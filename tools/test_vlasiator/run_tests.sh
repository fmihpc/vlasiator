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

create_verification_files=0

# choose tests to run
#run_tests=(1 2
run_tests=(1 )


# test 1
test_name[1]="Fluctuations"
test_cfg[1]="data/Fluctuations.cfg"

# test 1
test_name[2]="test_fp"
test_cfg[2]="data/test_fp.cfg"

#folder for reference data
reference_dir="reference"
#where the tests are run
run_dir="run"

#vlasiator binary
bin=vlasiator

# give here t%he variables you want to be tested
variables_name=( "rho" "rho_v" "rho_v" "rho_v" "B" "B" "B" "E" "E" "E" )
# and the corresponding components to variables, 
variables_components=( 0 0 1 2 0 1 2 0 1 2)
#arrays variables_name and variables_components should have same number of elements, e.g., 4th variable is variables_name[4] variables_components[4]=  

# choose the length of the test run: short/long
#short_long="short"





## code, no need to touch ##

#Go to PBS_O_WORKDIR if set (batch job), otherwise go to current folder
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
    # directory where the good results are for comparison
    if [ !  -d $reference_dir ]
        then
        mkdir $reference_dir 
    fi

    
    for run in ${run_tests[*]}
      do
# directory where the good results are for comparison
      result_dir=${reference_dir}/${test_name[$run]}
      if [ -d $result_dir ]
          then
          #we do not want to save old stuff
          #rm -rf $result_dir
          echo "delete"
      fi
      
      mkdir $result_dir
      cd  $result_dir
      echo "Compute reference for ${test_name[$run]} into $( pwd )"
      
      export OMP_NUM_THREADS=$t
      export MPICH_MAX_THREAD_SAFETY=funneled
      $run_command -n $p -d $t $bin --run_config=${test_cfg[$run]}
      
      #down to basic folder, strictly taken not needed
      cd $base_dir
    done
else

    if [ !  -d $run_dir ]
        then
        mkdir $run_dir 
    fi

#print header
    echo "------------------------------------------------------------------------------------------------------------" 
    echo "      Test case         |     variable     |          value           |           ok/bad                    "  


# loop over different test cases
    for run in ${run_tests[*]}
      do
      
      
# directory for test results
      vlsv_dir=${run_dir}/${test_name[$run]}

# Check if files for new vlsv files exists, if not create them
      if [ ! -d ${vlsv_dir} ]; then
          echo creating ${vlsv_dir}
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

##Compare test case with right solutions

# directory where the good results are for comparison
      result_dir=${reference_dir}/${test_name[$run]}
# loop over variables

      echo "------------------------------------------------------------------------------------------------------------"
      for i in ${!variables_name[*]}
        do

        value=$( vlsvdiff_DP ${result_dir}/grid.0000001.vlsv ${vlsv_dir}/grid.0000001.vlsv ${variables_name[$i]} ${variables_components[$i]} |grep "The relative 0-distance between both datasets" |gawk '{print $8}'  )


        a=$(awk -vx="$value" 'BEGIN{print (x > 0 ? x:-x)}')
        b=$(awk -vn="$a" 'BEGIN{print ( n<=10^(-14))?1:0 }')


        if [ "$b" == "1" ]; then
            status="ok"
        else
            status="bad"
        fi
#print the results      
        echo "      ${test_name[$run]}          ${variables_name[$i]} ${variables_components[$i]}            $value                $status              "
      done # loop over variables
      echo "------------------------------------------------------------------------------------------------------------"
    done # loop over tests

fi

