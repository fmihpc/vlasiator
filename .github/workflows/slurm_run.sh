#!/usr/bin/env bash

#This will echo the bash commands in the github view
set -x

if [[ -z $VLASIATOR_ARCH ]]; then
  echo "VLASIATOR_ARCH not set!"
  exit 1
fi

#Flags for core/node counts etc for compiling/heavier srun calls
declare -A core_flags
core_flags["turso-amd_GNU_MPICH"]="-n 1 -c 16"
core_flags["turso-v100_GNU_MPICH"]="-n 1 -c 16"
core_flags["pioneer"]="-n 1 -c 64"
core_flags["hile_gpu"]="-n 1 -c 16"
core_flags["hile_cpu"]="-n 1 -c 16"
core_flags["lumi_2503"]="-n 1 -c 16"
core_flags["lumi_2503_gpu"]="-n 1 -c 16"

if [[ -z ${core_flags[$VLASIATOR_ARCH]} ]]; then
  echo "No matching core_flags (and likely no constraints, etc) for VLASIATOR_ARCH $VLASIATOR_ARCH!"
  exit 1
fi

declare -A lib_flags
lib_flags["turso-amd_GNU_MPICH"]="-n 1 -c 16 --mem=64G"
lib_flags["turso-v100_GNU_MPICH"]="-n 1 -c 16 --mem=64G"
lib_flags["pioneer"]="-n 1 -c 16 --mem=64G"
lib_flags["hile_gpu"]="-n 1 -c 16 --mem=64G"
lib_flags["hile_cpu"]="-n 1 -c 16 --mem=64G"
lib_flags["lumi_2503"]="-n 1 -c 16 --mem=64G"
lib_flags["lumi_2503_gpu"]="-n 1 -c 16 --mem=64G"

if [[ -z ${lib_flags[$VLASIATOR_ARCH]} ]]; then
  echo "No matching lib_flags for VLASIATOR_ARCH $VLASIATOR_ARCH!"
  exit 1
fi

#Constraints for compiling stuff
declare -A constraint
constraint["turso-amd_GNU_MPICH"]="--constraint=amd"
constraint["turso-v100_GNU_MPICH"]="--gres=gpu:V100:1"
constraint["pioneer"]="-p pioneer -t 01:00:00" #not sure if pty needed for pioneer
constraint["hile_gpu"]="-C g"
constraint["hile_cpu"]="-C c"
constraint["lumi_2503"]="--partition=standard --account=project_462001599 -t 01:00:00"
# GPU node for a guaranteed ROCm install; the build itself needs no GPU.
# TMPDIR off the /tmp ramdisk: small-g nodes are shared, and when one of our jobs ends /tmp is cleaned under the others.
constraint["lumi_2503_gpu"]="--partition=small-g --account=project_462001599 --gpus-per-node=1 -t 01:00:00 --export=ALL,TMPDIR=$GITHUB_WORKSPACE"

if [[ -z ${constraint[$VLASIATOR_ARCH]} ]]; then
  echo "No matching constraint for VLASIATOR_ARCH $VLASIATOR_ARCH!"
  exit 1
fi

#Constraints used for smaller jobs like compiling/removing files/catting etc
declare -A constraint_small
constraint_small["turso-amd_GNU_MPICH"]="--constraint=amd -t 0:30:00"
constraint_small["turso-v100_GNU_MPICH"]="--constraint=V100"  # we need to build on a machine with a GPU because linking with the CUDA-aware MPICH needs the NVIDIA driver
constraint_small["pioneer"]="-p pioneer -t 01:00:00"
constraint_small["hile_gpu"]="-C g"
constraint_small["hile_cpu"]="-C c"
constraint_small["lumi_2503"]="--partition=standard --account=project_462001599 -t 00:15:00"
constraint_small["lumi_2503_gpu"]="--partition=small-g --account=project_462001599 --gpus-per-node=1 -t 00:15:00"

if [[ -z ${constraint_small[$VLASIATOR_ARCH]} ]]; then
  echo "No matching constraint_small for VLASIATOR_ARCH $VLASIATOR_ARCH!"
  exit 1
fi

#Memory flags for compiling, note that with --exclusive it is better to use --mem since --mem-per-cpu counts the whole node apparently
declare -A mem_flags

#Could be unified with core_flags but it's nice to have them here separately and maybe this separation has a future use
mem_flags["turso-amd_GNU_MPICH"]="--mem=40G"
mem_flags["turso-v100_GNU_MPICH"]="--mem=64G"  # ~4G per parallel nvcc job
mem_flags["pioneer"]=" "
mem_flags["hile_gpu"]="--mem=32G"
mem_flags["hile_cpu"]="--mem=32G"
mem_flags["lumi_2503"]="--mem=220G"
mem_flags["lumi_2503_gpu"]="--mem=64G"

if [[ -z ${mem_flags[$VLASIATOR_ARCH]} ]]; then
  echo "No matching mem_flags for VLASIATOR_ARCH $VLASIATOR_ARCH!"
  exit 1
fi

#Testpackage job time limits for quick (TP_QUICK=true) and full runs, overriding the #SBATCH -t of the CI test script.
#Arches without an entry keep the script's limit. Keep the full limit below GitHub's 6 h job timeout.
declare -A tp_time_quick
declare -A tp_time_full
tp_time_quick["turso-v100_GNU_MPICH"]="01:00:00"
tp_time_full["turso-v100_GNU_MPICH"]="05:00:00"

#Production compile flags
declare -A compile_flags_prod
declare -A compile_flags_tp

#This is for `debug: true` entry in workflow
if [[ $2 == "true" ]]; then
  #testpackage/debug compile flags
  compile_flags_prod["turso-amd_GNU_MPICH"]='COMPFLAGS="-DDEBUG_VLASIATOR -DDEBUG_SOLVERS -DDEBUG_IONOSPHERE -DHASHINATOR_DEBUG -DDEBUG_SPATIAL_CELL -DDEBUG_VMESH -DDEBUG_VBC -DDEBUG_ACC "'
  #Warp accessors only in tp, so both USE_WARPACCESSORS branches get a debug build.
  compile_flags_prod["turso-v100_GNU_MPICH"]='COMPFLAGS="-DDEBUG_VLASIATOR -DDEBUG_SOLVERS -DDEBUG_IONOSPHERE -DHASHINATOR_DEBUG -DDEBUG_SPATIAL_CELL -DDEBUG_VMESH -DDEBUG_VBC -DDEBUG_ACC "'
  compile_flags_tp["turso-v100_GNU_MPICH"]='COMPFLAGS="-DDEBUG_VLASIATOR -DDEBUG_SOLVERS -DDEBUG_IONOSPHERE -DHASHINATOR_DEBUG -DDEBUG_SPATIAL_CELL -DDEBUG_VMESH -DDEBUG_VBC -DDEBUG_ACC -DUSE_WARPACCESSORS "'
fi
#Module load part
modules="hostname; realpath modules/modules.sh; source modules/$VLASIATOR_ARCH.sh"

#--------------------------------different srun calls------------------------------

#0++++++++++++++++++++++++++++++0
#|         CLEANUP              |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "CLEANUP" ]]; then
  CLEAN_STRING=$(
    cat <<MORO
rm -rf libraries-$VLASIATOR_ARCH library-build testpackage libraries submodules testpackage_logs
rm -f libraries-$VLASIATOR_ARCH.tar.gz testpackage_check_description.txt testpackage-output-$VLASIATOR_ARCH.tar.gz metrics.txt stdout.txt stderr.txt testpackage_output_variables.txt
rm -f *.xml
echo "Cleaned up workspace"
MORO
  )
  srun ${constraint_small[$VLASIATOR_ARCH]} bash -c "$CLEAN_STRING"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         BUILD LIBS           |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "BUILD_LIBS" ]]; then
  srun ${constraint_small[$VLASIATOR_ARCH]} ${lib_flags[$VLASIATOR_ARCH]} -J build_libraries_CI bash -lc "$modules ; ./fetch_and_build_libraries.sh $VLASIATOR_ARCH"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         BUILD TOOLS          |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "BUILD_TOOLS" ]]; then
  #not constraint small since some platforms used gpu partition here but not build libs previously
  srun ${constraint[$VLASIATOR_ARCH]} --job-name CI_TOOLS_COMPILE --interactive --nodes=1 -n 1 -c 1 --mem=4G -t 0:10:0 bash -lc "$modules ; make vlsvextract vlsvdiff fluxfunction ; RC=\$?; sleep 10s; exit \$RC"
  exit $?
fi

COMPILE_STRING="make -j $(echo ${core_flags[$VLASIATOR_ARCH]} | grep -Po '\-c \d+' | grep -Po '\d+')"

#Escaped \$? and \$RC: must be expanded by the srun shell, not this one.
#compile_flags_* must directly prefix make, otherwise COMPFLAGS only reaches the first command of $modules.

#0++++++++++++++++++++++++++++++0
#|         COMPILE PROD         |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "COMPILE_PROD" ]]; then
  srun ${constraint[$VLASIATOR_ARCH]} --job-name CI_PROD_COMPILE --interactive ${mem_flags[$VLASIATOR_ARCH]} ${core_flags[$VLASIATOR_ARCH]} -t 0:30:0 bash -lc "$modules ; ${compile_flags_prod[$VLASIATOR_ARCH]} $COMPILE_STRING ; RC=\$?; sleep 10s; exit \$RC"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         COMPILE TP           |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "COMPILE_TP" ]]; then
  srun ${constraint[$VLASIATOR_ARCH]} --job-name CI_TP_COMPILE --interactive ${mem_flags[$VLASIATOR_ARCH]} ${core_flags[$VLASIATOR_ARCH]} -t 0:30:0 bash -lc "$modules ; ${compile_flags_tp[$VLASIATOR_ARCH]} $COMPILE_STRING testpackage ; RC=\$?; sleep 10s; exit \$RC"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         RUN TP               |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "RUN_TP" ]]; then
  if [[ "$VLASIATOR_ARCH" == "turso-amd_GNU_MPICH" || "$VLASIATOR_ARCH" == "turso-v100_GNU_MPICH" || "$VLASIATOR_ARCH" == "hile_cpu" || "$VLASIATOR_ARCH" == "lumi_2503" ]]; then

    #Platform specific expections can be added here
    if [[ "$VLASIATOR_ARCH" == "turso-amd_GNU_MPICH" || "$VLASIATOR_ARCH" == "turso-v100_GNU_MPICH" ]]; then
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GITHUB_WORKSPACE/libraries-$VLASIATOR_ARCH/lib
    fi

    chmod +x $GITHUB_WORKSPACE/vlasiator
    chmod +x $GITHUB_WORKSPACE/vlsv*_DP
    echo $PWD
    cd testpackage
    # Delete any old run directories
    rm -rf run_20*
    echo "TP_PID=$$" >> $GITHUB_OUTPUT
    if [[ $TP_QUICK == true ]]; then
      tp_time=${tp_time_quick[$VLASIATOR_ARCH]}
    else
      tp_time=${tp_time_full[$VLASIATOR_ARCH]}
    fi
    tp_jobid=$(sbatch --parsable -W ${tp_time:+-t $tp_time} -o testpackage_run_output.txt -J "CI_TP_$$" ./small_test_${VLASIATOR_ARCH}_github_ci.sh) || true
    tp_jobid=${tp_jobid%%;*}
    # A job that Slurm stops (time limit, node failure, memory) never gets to record its own failure.
    tp_state=$(sacct -n -X -P -o State -j "$tp_jobid" 2>/dev/null | head -n 1)
    echo "Testpackage job ${tp_jobid:-?} ended in state ${tp_state:-unknown}"
    tp_problem=""
    if [[ -z $tp_jobid ]]; then
      tp_problem="could not be submitted"
    elif [[ $tp_state =~ ^(TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|DEADLINE|BOOT_FAIL) ]]; then
      tp_problem="was stopped by Slurm ($tp_state)"
    fi
    if [[ -n $tp_problem ]]; then
      echo -e "\n:red_circle: The testpackage job $tp_problem; tests it did not reach have no results." >> $GITHUB_STEP_SUMMARY
      if [[ ! -f $GITHUB_WORKSPACE/testpackage_output_variables.txt ]]; then
        echo "summary=The testpackage job $tp_problem." >> $GITHUB_WORKSPACE/testpackage_output_variables.txt
        echo "conclusion=failure" >> $GITHUB_WORKSPACE/testpackage_output_variables.txt
      fi
      touch $GITHUB_WORKSPACE/testpackage_failed
    fi
    PARSE_OUTPUT_CMD=$(
      cat <<MORO
echo "Job finished, checking output."
cat testpackage_run_output.txt
cat $GITHUB_STEP_SUMMARY > $GITHUB_WORKSPACE/testpackage_check_description.txt
cd $GITHUB_WORKSPACE
ls -halB testpackage_check_description.txt
tar -czf testpackage-output-$VLASIATOR_ARCH.tar.gz testpackage_check_description.txt testpackage_output_variables.txt
MORO
    )
    srun --job-name CI_package_results ${constraint_small[$VLASIATOR_ARCH]} -N 1 -n 1 -c 1 --mem=3G bash -c "$PARSE_OUTPUT_CMD"
    if [ -f $GITHUB_WORKSPACE/testpackage_failed ]; then
      # Fail this step if any test failed.
      exit 1
    fi
  else
    echo "RUN_TP for $VLASIATOR_ARCH not implemented"
    exit 1
  fi
fi

#0++++++++++++++++++++++++++++++0
#|         RUN FLUXTEST         |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "FLUXTEST" ]]; then
  FLUXTEST_STRING="$modules ; $GITHUB_WORKSPACE/fluxfunction testpackage/run_*/Magnetosphere_small/bulk.0000001.vlsv equatorial.bin ; $GITHUB_WORKSPACE/fluxfunction testpackage/run_*/Magnetosphere_polar_small/bulk.0000001.vlsv polar.bin"
  chmod +x $GITHUB_WORKSPACE/fluxfunction
  srun ${constraint_small[$VLASIATOR_ARCH]} --job-name CI_FLUXTEST -N 1 -n 1 -c 1 --mem=2G -t 0:10:0 bash -lc "$FLUXTEST_STRING"

  #Platform specific expections can be added here
  if [[ "$VLASIATOR_ARCH" == "lumi_2503" ]]; then
    export CI_REFERENCE_LOCATION=/projappl/project_462001599/testpackage/CI_reference/;
  elif [[ "$VLASIATOR_ARCH" == "turso-v100_GNU_MPICH" ]]; then
    export CI_REFERENCE_LOCATION=/turso/group/spacephysics/vlasiator/testpackage/CI_gpu_reference/;
  else
    export CI_REFERENCE_LOCATION=/turso/group/spacephysics/vlasiator/testpackage/CI_reference/;
  fi

  diff -q equatorial.bin $CI_REFERENCE_LOCATION/equatorial.bin || if [ $? -eq 1 ]; then true; else false; fi
  diff -q polar.bin      $CI_REFERENCE_LOCATION/polar.bin      || if [ $? -eq 1 ]; then true; else false; fi
fi
