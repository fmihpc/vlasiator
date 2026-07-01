#!/usr/bin/env bash

set -x

if [[ -z $VLASIATOR_ARCH ]]; then
  echo "VLASIATOR_ARCH not set!"
  exit 1
fi
#Flags for core/node counts etc for compiling/heavier srun calls
declare -A core_flags
core_flags["carrington_gcc_openmpi"]="-n 1 -c 16"
core_flags["ukko_dgx"]="-n 1 -c 64"
core_flags["pioneer"]="-n 1 -c 64"
core_flags["hile_gpu"]="-n 1 -c 16"
core_flags["hile_cpu"]="-n 1 -c 16"

#Constraints for compiling stuff
declare -A constraint
constraint["carrington_gcc_openmpi"]="--constraint=ukko|carrington -p short"
# constraint["ukko_dgx"]="--constraint="v100" -gres="gpu:V100" -p gpu"
constraint["ukko_dgx"]="--constraint=ukko -p gpu"
constraint["pioneer"]="--constraint=pioneer -p pioneer" #not sure if pty needed for pioneer
constraint["hile_gpu"]="-C g"
constraint["hile_cpu"]="-C c"

#Constraints used for smaller jobs like compiling/removing files/catting etc
declare -A constraint_small
constraint_small["carrington_gcc_openmpi"]="--constraint=ukko|carrington"
constraint_small["ukko_dgx"]="--constraint=ukko"
constraint_small["pioneer"]="--constraint=pioneer"
constraint_small["hile_gpu"]="-C g"
constraint_small["hile_cpu"]="-C c"

#Memory flags for compiling, note that with --exclusive it is better to use --mem since --mem-per-cpu counts the whole node apparently
declare -A mem_flags
mem_flags["carrington_gcc_openmpi"]="--mem=40G"
mem_flags["ukko_dgx"]="--mem=64G"
mem_flags["pioneer"]=""
mem_flags["hile_gpu"]="--mem=32G"
mem_flags["hile_cpu"]="--mem=32G"

declare -A compile_flags_prod
compile_flags_prod["ukko_dgx"]='COMPFLAGS="-DDEBUG_VLASIATOR -DDEBUG_SOLVERS -DDEBUG_IONOSPHERE -DHASHINATOR_DEBUG -DDEBUG_SPATIAL_CELL -DDEBUG_VMESH -DDEBUG_VBC -DDEBUG_ACC "'

declare -A compile_flags_tp
compile_flags_tp["ukko_dgx"]='COMPFLAGS="-DDEBUG_VLASIATOR -DDEBUG_SOLVERS -DDEBUG_IONOSPHERE -DHASHINATOR_DEBUG -DDEBUG_SPATIAL_CELL -DDEBUG_VMESH -DDEBUG_VBC -DDEBUG_ACC -DUSE_WARPACCESSOR "'

#Flags for bigger sruns, like compiling
# flags="--interactive -n 1 ${platform_flags[$VLASIATOR_ARCH]} ${constraint[$VLASIATOR_ARCH]}"
#flags for smaller sruns, like removing files/small tests
# small_flags="--interactive ${constraint_small[$VLASIATOR_ARCH]} -n 1 -c 1"
#Module load part
modules="hostname; realpath modules/modules.sh; source modules/$VLASIATOR_ARCH.sh"

#--------------------------------different srun calls------------------------------

#todo add echo of the constraints etc

#0++++++++++++++++++++++++++++++0
#|         CLEANUP              |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "CLEANUP" ]]; then
  CLEAN_STRING=$(
    cat <<MORO
rm -rf libraries-$VLASIATOR_ARCH library-build testpackage
rm -f libraries-$VLASIATOR_ARCH.tar.gz testpackage_check_description.txt testpackage-output.tar.gz metrics.txt stdout.txt stderr.txt testpackage_output_variables.txt
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
  srun ${constraint_small[$VLASIATOR_ARCH]} ${platform_flags[$VLASIATOR_ARCH]} --mem=8G -J build_libraries_CI bash -lc "$modules ; ./fetch_and_build_libraries.sh $VLASIATOR_ARCH"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         BUILD TOOLS          |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "BUILD_TOOLS" ]]; then
  srun ${constraint[$VLASIATOR_ARCH]} --job-name CI_TOOLS_COMPILE --interactive --nodes=1 -n 1 -c 1 --mem=4G -t 0:10:0 bash -c "$modules ; make vlsvextract vlsvdiff fluxfunction ; sleep 10s"
  exit $?
fi

COMPILE_STRING="$modules ; make -j $(echo ${core_flags[$VLASIATOR_ARCH]} | grep -Po '\-c \d+' | grep -Po '\d+')"

#0++++++++++++++++++++++++++++++0
#|         COMPILE PROD         |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "COMPILE_PROD" ]]; then
  srun ${constraint[$VLASIATOR_ARCH]} --job-name CI_PROD_COMPILE --interactive ${mem_flags[$VLASIATOR_ARCH]} ${core_flags[$VLASIATOR_ARCH]} -t 0:10:0 bash -c "${compile_flags_prod[$VLASIATOR_ARCH]} $COMPILE_STRING ; sleep 10s"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         COMPILE TP           |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "COMPILE_TP" ]]; then
  srun ${constraint[$VLASIATOR_ARCH]} --job-name CI_TP_COMPILE --interactive ${mem_flags[$VLASIATOR_ARCH]} ${core_flags[$VLASIATOR_ARCH]} -t 0:10:0 bash -c "${compile_flags_tp[$VLASIATOR_ARCH]} $COMPILE_STRING testpackage ; sleep 10s"
  exit $?
fi

#0++++++++++++++++++++++++++++++0
#|         RUN TP               |
#0++++++++++++++++++++++++++++++0
if [[ $1 == "RUN_TP" ]]; then
  if [[ "$VLASIATOR_ARCH" == "carrington_gcc_openmpi" || "$VLASIATOR_ARCH" == "hile_cpu" ]]; then

    #Platform specific expections can be added here
    if [[ "$VLASIATOR_ARCH" == "carrington_gcc_openmpi" ]]; then
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GITHUB_WORKSPACE/libraries-carrington_gcc_openmpi/lib
    fi

    chmod +x $GITHUB_WORKSPACE/vlasiator
    chmod +x $GITHUB_WORKSPACE/vlsv*_DP
    echo $PWD
    cd testpackage
    # Delete any old run directories
    rm -rf run_20*

    sbatch -W -o testpackage_run_output.txt ./small_test_${VLASIATOR_ARCH}_github_ci.sh || true
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
    srun --job-name CI_package_results ${constraint_small[$VLASIATOR_ARCH]} -N 1 -c 1 --mem=3G bash -c "$PARSE_OUTPUT_CMD"
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
  srun ${constraint_small[$VLASIATOR_ARCH]} --job-name CI_FLUXTEST -N 1 -c 1 --mem=2G -t 0:10:0 bash -c "$FLUXTEST_STRING"

  diff -q equatorial.bin /turso/group/spacephysics/vlasiator/testpackage/CI_reference/equatorial.bin || if [ $? -eq 1 ]; then true; else false; fi
  diff -q polar.bin /turso/group/spacephysics/vlasiator/testpackage/CI_reference/polar.bin || if [ $? -eq 1 ]; then true; else false; fi
fi

##0++++++++++++++++++++++++++++++0
##|         RESULTS              |
##0++++++++++++++++++++++++++++++0
#if [[ $1 == $PARSE_OUTPUT_CMD ]]; then
#  srun --job-name CI_package_results ${constraint_small} -c 1 -N 1 --mem=2G bash -c "$PARSE_OUTPUT_CMD"
#  exit $?
#fi
