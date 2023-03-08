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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <iostream>

#include "common.h"
#include "mpi.h"

#include "cuda_context.cuh"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

//CUcontext cuda_acc_context;
cudaStream_t cudaStreamList[MAXCPUTHREADS];
Realf *returnRealf[MAXCPUTHREADS];

__host__ void cuda_set_device() {

#ifdef _OPENMP
   const uint maxThreads = omp_get_max_threads();
#else
   const uint maxThreads = 1
#endif

   int deviceCount;
   cudaGetDeviceCount(&deviceCount);
   printf("CUDA device count %d with %d threads/streams\n",deviceCount,maxThreads);

   /* Create communicator with one rank per compute node to identify which GPU to use */
   int amps_size;
   int amps_rank;
   int amps_node_rank;
   int amps_node_size;
   // int amps_write_rank;
   // int amps_write_size;
   MPI_Comm amps_CommWorld = MPI_COMM_NULL;
   MPI_Comm amps_CommNode = MPI_COMM_NULL;

   MPI_Comm_dup(MPI_COMM_WORLD, &amps_CommWorld);
   MPI_Comm_size(amps_CommWorld, &amps_size);
   MPI_Comm_rank(amps_CommWorld, &amps_rank);

   /* Create communicator with one rank per compute node */
#if MPI_VERSION >= 3
   MPI_Comm_split_type(amps_CommWorld, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &amps_CommNode);
#else
   /* Split the node level communicator based on Adler32 hash keys of processor name */
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int namelen;
   MPI_Get_processor_name(processor_name, &namelen);
   uint32_t checkSum = Adler32((unsigned char*)processor_name, namelen);
   /* Comm split only accepts non-negative numbers */
   /* Not super great for hashing purposes but hoping MPI-3 code will be used on most cases */
   checkSum &= INT_MAX;
   MPI_Comm_split(amps_CommWorld, checkSum, amps_rank, &amps_CommNode);
#endif

   MPI_Comm_rank(amps_CommNode, &amps_node_rank);
   MPI_Comm_size(amps_CommNode, &amps_node_size);
   //cerr << "(Grid) rank " << amps_rank << " is noderank "<< amps_node_rank << " of "<< amps_node_size << endl;

   if (amps_node_rank >= deviceCount) {
      std::cerr<<"Error, attempting to use CUDA device beyond available count!"<<std::endl;
      abort();
   }
   if (amps_node_size > deviceCount) {
      std::cerr<<"Error, MPI tasks per node exceeds available CUDA device count!"<<std::endl;
      abort();
   }
   cudaSetDevice(amps_node_rank);

   // Pre-generate streams, allocate return pointers
   for (uint i=0; i<maxThreads; ++i) {
      cudaStreamCreate(&(cudaStreamList[i]));
      HANDLE_ERROR( cudaMalloc((void**)&returnRealf[i], sizeof(Realf)) );
   }

   // Using just a single context for whole MPI task
}

__host__ void cuda_clear_device() {
   // Destroy streams
#ifdef _OPENMP
   const uint maxThreads = omp_get_max_threads();
#else
   const uint maxThreads = 1
#endif
   for (uint i=0; i<maxThreads; ++i) {
      cudaStreamDestroy(cudaStreamList[i]);
      HANDLE_ERROR( cudaFree(*returnRealf) );
   }
}

__host__ cudaStream_t cuda_getStream() {
#ifdef _OPENMP
   const uint thread_id = omp_get_thread_num();
#else
   const uint thread_id = 0;
#endif
   return cudaStreamList[thread_id];
}
