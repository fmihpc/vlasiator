/*
 * This file is part of Vlasiator.
 * Copyright 2010-2022 Finnish Meteorological Institute and University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <omp.h>
#include <stdio.h>

#include "common.h"

#include "cuda_context.cuh"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

//CUcontext cuda_acc_context;
cudaStream_t cudaStreamList[MAXCPUTHREADS];

__host__ void cuda_set_device() {

   const uint maxThreads = omp_get_max_threads();

   int deviceCount;
   cudaGetDeviceCount(&deviceCount);
   printf("CUDA device count %d with %d threads/streams\n",deviceCount,maxThreads);

   int gpuid = 0;
   if (gpuid >= deviceCount) {
      printf("Error, attempting to use CUDA device beyond available count!\n");
   }
   cudaSetDevice(gpuid);

   // Pre-gnerate streams
   for (uint i=0; i<maxThreads; ++i) {
      cudaStreamCreate(&(cudaStreamList[i]));
   }

   // Using just a single context for whole MPI task

   // CUdevice cuDevice;
   // CUresult result;
   // result = cuDeviceGet(&cuDevice,gpuid);
   // printf("Active CUDA device %d\n",cuDevice);
   // result = cuCtxCreate(&cuda_acc_context, 0,cuDevice);
   // printf("Created CUDA context %ld \n",cuda_acc_context);

   // Loop to create one CUDA context per thread
   // for (uint i=0; i<omp_get_max_threads(); ++i) {
   //    result = cuCtxCreate(&cuda_thread_context[i], 0,cuDevice);
   //    printf("Created CUDA context %d %d \n",i,cuda_thread_context[i]);
   // }
   // const uint cuda_async_queue_id = omp_get_thread_num();
   // CUcontext cuContext;
   //result = cuDeviceGet(&cuDevice,gpuid);
   //checkError(result);
   //result = cuCtxCreate(&cuContext, 0,cuDevice);
   //checkError(result);
   // result = cuCtxCreate(&cuContext, CU_CTX_MAP_HOST|CU_CTX_BLOCKING_SYNC, cuDevice);
   // checkError(result);
}

__host__ void cuda_clear_device() {
   // Destroy streams
   const uint maxThreads = omp_get_max_threads();
   for (uint i=0; i<maxThreads; ++i) {
      cudaStreamDestroy(cudaStreamList[i]);
   }
}

__host__ void cudaAllocateBlockData(
   Realf** dev_blockData,
   Real** dev_parameters,
   vmesh::LocalID blockCount
   ) {
   HANDLE_ERROR( cudaMalloc((void**)dev_blockData, blockCount*WID3*sizeof(Realf)) );
   HANDLE_ERROR( cudaMalloc((void**)dev_parameters, blockCount*BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Realf)) );
}

__host__ void cudaDeallocateBlockData(
   Realf** dev_blockData,
   Real** dev_parameters
   ) {
   HANDLE_ERROR( cudaFree(*dev_blockData) );
   HANDLE_ERROR( cudaFree(*dev_parameters) );
}

/** In the following functions we do not yet copy over blockparameters as they are only edited / used on the host for now.
 **/

__host__ void cuda_HtoD_BlockData(
   Realf* dev_blockData,
   Realf* blockData,
   Real* dev_parameters,
   Real* parameters,
   vmesh::LocalID blockCount
   ) {
   const uint thread_id = omp_get_thread_num();
   cudaHostRegister(blockData, blockCount*WID3*sizeof(Realf),cudaHostRegisterDefault);
   //cudaHostRegister(parameters, blockCount*BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Realf),cudaHostRegisterDefault);

   cudaMemcpyAsync(dev_blockData, blockData, blockCount*WID3*sizeof(Realf), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
   //cudaMemcpyAsync(dev_parameters, parameters, blockCount*BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Realf), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
}

__host__ void cuda_DtoH_BlockData(
   Realf* dev_blockData,
   Realf* blockData,
   Real* dev_parameters,
   Real* parameters,
   vmesh::LocalID blockCount
   ) {
   const uint thread_id = omp_get_thread_num();
   cudaHostRegister(blockData, blockCount*WID3*sizeof(Realf),cudaHostRegisterDefault);
   //cudaHostRegister(parameters, blockCount*BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Realf),cudaHostRegisterDefault);

   cudaMemcpyAsync(blockData, dev_blockData, blockCount*WID3*sizeof(Realf), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]);
   //cudaMemcpyAsync(parameters, dev_parameters, blockCount*BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Realf), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]);
}

__host__ void cuda_unregister_BlockData(
   Realf* blockData,
   Real* parameters
   ) {
   cudaHostUnregister(blockData);
   cudaHostUnregister(parameters);
}
