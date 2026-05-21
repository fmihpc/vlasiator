#!/bin/bash -l
#SBATCH --job-name=gpuTest
#SBATCH --account=project_2017838
#SBATCH --time=00:10:00
#SBATCH --partition=gpupilot
#SBATCH --ntasks-per-node=1 --cpus-per-task=7
#SBATCH --nodes=1
#SBATCH --dependency=singleton
#SBATCH --gres=gpu:gh200:1

# If 1, the reference vlsv files are generated
# If 0 then we check the v1
create_verification_files=0

# Folder for all reference data 
reference_dir="/scratch/project_2017838/mikael/testpackage/dev"
cd $SLURM_SUBMIT_DIR

bin="../vlasiator"
diffbin="../vlsvdiff_DP"

# compare agains which revision?
reference_revision="9e172e4b2937f9faca1afb0fabfc85a5408fb08a__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER__DACC_SEMILAG_PQM__DTRANS_SEMILAG_PPM__DDP__DDPF__DVEC4D_AGNER"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
module load boost papi gcc/13.4.0 cuda/12.6.3
export LD_LIBRARY_PATH=/scratch/project_2017838/mikael/vlasiator/libraries-roihu-gpu/lib:$LD_LIBRARY_PATH
ht=1    #hyper threads per physical core
t=$OMP_NUM_THREADS     #threads per process
nodes=$SLURM_NNODES

# set up GPU/CPU bindings
cat << EOF > select_gpu_${SLURM_JOB_ID}
#!/bin/bash
export CUDA_VISIBLE_DEVICES=\$SLURM_LOCALID
export OMP_NUM_THREADS=7
LD_PRELOAD=/scratch/project_2017838/mikael/MPI_Pancake/libmpipancake.so exec \$*
EOF
chmod +x ./select_gpu_${SLURM_JOB_ID}

run_command="srun -n 1 ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
small_run_command="srun -n 1 ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
run_command_tools="srun -n 1 ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID}"

#run_command="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
#small_run_command="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "
#run_command_tools="srun -n 1 --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID}"

umask 007

echo "Running on ${SLURM_NTASKS} mpi tasks with ${OMP_NUM_THREADS} threads per task on ${SLURM_JOB_NUM_NODES} nodes"

# Define test
source test_definitions_small.sh
wait
# Run tests
source run_tests.sh
wait 2

# cleanup
rm -f ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID}	
