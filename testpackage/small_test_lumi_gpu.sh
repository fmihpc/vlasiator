#!/bin/bash -l
#SBATCH --job-name=vlasi_g_tp
#SBATCH --partition=dev-g
##SBATCH --reservation=training_gpu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8

#SBATCH --time=02:00:00
#SBATCH --account=project_465000538
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --hint=nomultithread

# If 1, the reference vlsv files are generated
# If 0 then we check the v1
create_verification_files=0

# Folder for all reference data 
reference_dir="/scratch/project_465000538/battarbe/testpackage/"
cd $SLURM_SUBMIT_DIR

#bin="/scratch/project_465000538/battarbe/testpackage/vlasiator_NOMAD_WID8_tp"
bin="/scratch/project_465000538/battarbe/testpackage/vlasiator_NOMAD2_WID4_tp"
diffbin="/scratch/project_465000538/battarbe/testpackage/vlsvdiff_DP_NOMAD"

# compare agains which revision?
reference_revision="current"

# set up GPU/CPU bindings
cat << EOF > select_gpu
#!/bin/bash
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
export OMP_NUM_THREADS=7
exec \$*
EOF
chmod +x ./select_gpu
# this should set the ordering correctly: "4 5 2 3 6 7 0 1"
CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
CPU_BIND="${CPU_BIND},7e0000,7e000000"
CPU_BIND="${CPU_BIND},7e,7e00"
CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

module load LUMI/22.08 
module load partition/G
module load cpeAMD
module load rocm/5.3.3
module load Boost/1.79.0-cpeAMD-22.08
module list

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_SUPPORT_ENABLED=1
#Turn on managed memory paging
export HSA_XNACK=1
# use extra threads for MPI in background
export MPICH_ASYNC_PROGRESS=1


#Command for running tests and diffs
run_command="srun --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu "
small_run_command="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu "
run_command_tools="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu"

umask 007

echo "Running on ${SLURM_NTASKS} mpi tasks with ${OMP_NUM_THREADS} threads per task on ${SLURM_JOB_NUM_NODES} nodes"

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
wait 2

# hello jobstep verification
srun --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu ${SLURM_SUBMIT_DIR}/hello_jobstep
#cleanup
rm -f ${SLURM_SUBMIT_DIR}/select_gpu
