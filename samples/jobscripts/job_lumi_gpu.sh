#!/bin/bash -l
#SBATCH --job-name=vla_weakscaling
##SBATCH --partition=small-g
##SBATCH --time=24:00:00
#SBATCH --partition=dev-g
#SBATCH --time=3:00:00

#SBATCH --nodes=1 # max 32 nodes on dev-g
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:8 #GPUs per node

#SBATCH --account=project_462000358
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --hint=nomultithread

# broken on LUMI
##SBATCH --cpus-per-task=64 # this is broken with GPUs! Do not activate.
##SBATCH --gpu-bind=per_task:1 # also broken

export actualntasks=1

echo "using ${SLURM_CPUS_PER_TASK} CPUs per task"
ulimit -c unlimited
export PHIPROF_PRINTS=detailed,full
umask 007

module load LUMI/24.03; module load partition/G; module load cpeAMD; module load rocm/6.2.2; module load Boost/1.83.0-cpeAMD-24.03; module load papi/7.1.0.1
module list
export OMP_NUM_THREADS=6

cd $SLURM_SUBMIT_DIR

export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_SUPPORT_ENABLED=1
# Turn on XNACK support
export HSA_XNACK=1
# use extra threads for MPI in background
export MPICH_ASYNC_PROGRESS=1
export GPU_MAX_HW_QUEUES=14

export BIN="../vlasiator"
export CFG="../Flowthrough_amr.cfg"


# set up GPU/CPU bindings
cat << EOF > select_gpu_${SLURM_JOB_ID}
#!/bin/bash
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
export OMP_NUM_THREADS=6
LD_PRELOAD=/users/marbat/git/vlasiator-mempool/libpreload-me_622.so exec ${BIN} \$*
EOF


# # Alternatively, one can run profiling on one task only (for a limited number of timesteps) with this replacement:
# export CALL_REG="${BIN} --gridbuilder.timestep_max=20"
# export CALL_PROF="rocprof --trace-period 150s:15s:300s --roctx-trace --hip-trace --parallel-kernels --stats --basenames on ${BIN} --gridbuilder.timestep_max=20"
# # trace period: initial wait : trace period : interval between trace periods
# # other useful traces:
# # --sys-trace --hsa-trace
# export OMNITRACE_PERFETTO_ANNOTATIONS="false"

# cat << EOF > select_gpu_${SLURM_JOB_ID}
# #!/bin/bash
# export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
# if [ \$SLURM_PROCID == 0 ]; then
#     # Here profiler call
#     exec \${CALL_PROF} \$*
# else
#     # Regular launch call
#     exec \${CALL_REG} \$*
# fi
# EOF

chmod +x ./select_gpu_${SLURM_JOB_ID}

# this should set the ordering correctly: "4 5 2 3 6 7 0 1"
export CPU_BIND="mask_cpu:7e000000000000,7e00000000000000,,7e0000,7e000000,7e,7e00,7e00000000,7e0000000000"
ASRUN="srun --ntasks=${actualntasks} --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu_${SLURM_JOB_ID} "

# optional debug information
#srun rocm-smi --showtopo

# Store the version to the slurm output
srun -n 1 --cpu-bind=${CPU_BIND} ${BIN} --version

# Run the simulation
${ASRUN} --run_config ${CFG}

# cleanup
date
rm ./select_gpu_${SLURM_JOB_ID}


