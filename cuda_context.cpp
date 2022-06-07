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
//#include <stdlib.h>
//#include <string.h>

#include "common.h"

#include "cuda_context.cuh"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"


CUcontext cuda_acc_context;

__host__ void cuda_set_device() {

   int deviceCount;
   cudaGetDeviceCount(&deviceCount);
   printf("CUDA device count %d\n",deviceCount);
   
   int gpuid = 0;
   if (gpuid >= deviceCount) {
      printf("Error, attempting to use CUDA device beyond available count!\n");
   }
   cudaSetDevice(gpuid);

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

__host__ void cudaAllocateBlockData(
   Realf* dev_blockData,
   vmesh::LocalID blockCount
   ) {
   HANDLE_ERROR( cudaMalloc((void**)&dev_blockData, blockCount*WID3*sizeof(Realf)) );
}

__host__ void cudaDeallocateBlockData(
   Realf* dev_blockData
   ) {
   HANDLE_ERROR( cudaFree(dev_blockData) );
}
