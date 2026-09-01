#!/bin/bash -l
#SBATCH --time=00:15:00
#SBATCH --job-name=vla_gpu_weakscale_prof
#SBATCH --account=project_2004522
#SBATCH --gres=gpu:a100:4  #set last value to match GPUs per node
##SBATCH --partition=gpumedium
##SBATCH --partition=gpusmall
#SBATCH --partition=gputest
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --exclusive

export actualntasks=1
# Partition options:
# gputest: max 15 minutes, max 4 GPUs, max 1 node
# gpusmall: max 36 hours, max 2 GPUs, max 1/2 node
# gpumedium: max 36 hours, max 24 GPUs, max 6 node
# always reserve no more than 32 CPU cores per GPU
#
# Single GPU is ok to run with 1 task and non-exclusive,
# but then turn off CPU core binds.
# To use CPU core binds, one needs whole nodes, so
# that means whole nodes with 2-4 GPUs.

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores
export PHIPROF_PRINTS=compact,full

nodes=$SLURM_NNODES
export UCX_TLS=ud,ud_v
export OMPI_MCA_coll=^hcoll
umask 007
ulimit -c unlimited

module purge
module load gcc/10.4.0 openmpi/4.1.5-cuda cuda/12.1.1 boost/1.82.0-mpi papi/7.1.0

BIN="./vlasiator"
CFG="./Flowthrough_amr.cfg"

# store version information
srun -n 1 ${BIN} --version

# set up CPU mask
CPU_BIND="mask_cpu"
CPU_BIND="${CPU_BIND}:000f000f"
CPU_BIND="${CPU_BIND},00f000f0"
CPU_BIND="${CPU_BIND},0f000f00"
CPU_BIND="${CPU_BIND},f000f000"

# # set up GPU/CPU bindings
cat << EOF > select_gpu_${SLURM_JOB_ID}
#!/bin/bash
export CUDA_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF


# Alternatively, one can run profiling on one task only (for a limited number of timesteps) with this replacement:
# export CALL_REG="${BIN} --gridbuilder.timestep_max=30"
# export CALL_PROF="nsys profile -y 30 -d 20 -w true -t nvtx,cuda -s none --cuda-um-gpu-page-faults=true --stats=true ${BIN} --gridbuilder.timestep_max=30"
# cat << EOF > select_gpu_${SLURM_JOB_ID}
# #!/bin/bash
# export CUDA_VISIBLE_DEVICES=\$SLURM_LOCALID
# if [ \$SLURM_PROCID == 0 ]; then
#     # Here profiler call
#     exec \${CALL_PROF} \$*
# else
#     # Regular launch call
#     exec \${CALL_REG} \$*
# fi
# EOF


chmod +x ./select_gpu_${SLURM_JOB_ID}

srun -n $actualntasks --cpu-bind=${CPU_BIND} ./select_gpu_${SLURM_JOB_ID} --run_config ${CFG}

#clean up
rm select_gpu_${SLURM_JOB_ID}

sleep 1

