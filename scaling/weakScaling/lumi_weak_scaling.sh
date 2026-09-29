#!/bin/bash
#SBATCH --nodes=__NODES__
#SBATCH --ntasks-per-node=__TASKS_PER_NODE__
#SBATCH --gpus-per-node=__GPUS_PER_NODE__

#SBATCH --job-name=weak_scaling___NTASKS__
#SBATCH --output=weak_scaling___GPUS__.out
#SBATCH --partition=__PARTITION__

##SBATCH --time=24:00:00
##SBATCH --time=03:00:00
#SBATCH --time=00:10:00
#SBATCH --account=project_462000358
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --hint=nomultithread

# If 1, the reference vlsv files are generated
# If 0 then we check the v1
create_verification_files=0

# Folder for all reference data 
reference_dir="/scratch/project_462000358/testpackage_2025_06/"
cd $SLURM_SUBMIT_DIR

bin="../../vlasiator"
#bin="./vlasiator_WID8"
#diffbin="/scratch/project_462000358/testpackage_2025_06/vlsvdiff_gpu_DP"
diffbin="../../vlsvdiff_DP"

# compare agains which revision?
reference_revision="current"
t=6
export OMP_NUM_THREADS=6

# place before exec
#LD_PRELOAD=/users/marbat/git/vlasiator-mempool/libpreload-me_622.so

# set up GPU/CPU bindings
cat << EOT > select_gpu_${SLURM_JOB_ID}
#!/bin/bash
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
export OMP_NUM_THREADS=6
#LD_PRELOAD="__CWD__/../../libhooks.so:__CWD__/../../libpreload-me_622.so" exec \$*
LD_PRELOAD=__CWD__/../../libpreload-me_622.so exec \$*
exec \$*
EOT
chmod +x ./select_gpu_${SLURM_JOB_ID}
# this should set the ordering correctly: "4 5 2 3 6 7 0 1"
CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"


module load LUMI/24.03; module load partition/G; module load cpeAMD; module load rocm/6.2.2; module load Boost/1.83.0-cpeAMD-24.03; module load papi/7.1.0.1

# module load LUMI/24.03
# module load partition/G
# module load cpeAMD
# #module load rocm/6.0.3
# module load rocm/6.2.2
# module load Boost/1.83.0-cpeGNU-24.03
# module load papi/7.1.0.1
module list

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_SUPPORT_ENABLED=1
# Turn off forced managed memory paging
export HSA_XNACK=1
# use extra threads for MPI in background
#export MPICH_ASYNC_PROGRESS=1
# allow 16 in-parallel queues
export GPU_MAX_HW_QUEUES=16
# No MPI testing for now
run_command="srun -n __NTASKS__ --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
small_run_command="srun -n __NTASKS__ --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
run_command_tools="srun -n __NTASKS__ --cpu-bind=${CPU_BIND} "
# no mempool or GPU binding for tools

umask 007

echo "Running on ${SLURM_NTASKS} mpi tasks with ${OMP_NUM_THREADS} threads per task on ${SLURM_JOB_NUM_NODES} nodes"

# Define test

## Define test and runs

if [ ! -f $bin ]
then
   echo Executable $bin does not exist
   exit
fi

# where the tests are run
run_dir="./__FOLDER__/__CFG__/run"

# where the directories for different tests, including cfg and other needed data files are located 
test_dir="./__FOLDER__/"

# Counter for creating tests
index=1

#######
# ACCELERATION TESTS (1..5)
#######

# 13 Large AMR translation flowthrough test
test_name[${index}]="__CFG__"
comparison_vlsv[${index}]="bulk.0000001.vlsv"
comparison_phiprof[${index}]="phiprof_0.txt"
variable_names[${index}]="proton/vg_rho proton/vg_v proton/vg_v proton/vg_v fg_b fg_b fg_b fg_e fg_e fg_e"
variable_components[${index}]="0 0 1 2 0 1 2 0 1 2"

# Alternatively, set tests manually, e.g.
run_tests=( 1 )

wait
# Run tests
source run_tests.sh
wait 2

# hello jobstep verification
# srun --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu ${SLURM_SUBMIT_DIR}/hello_jobstep

# cleanup
rm -f ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID}

EOF
