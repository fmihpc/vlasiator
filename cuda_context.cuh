/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#ifndef CUDA_CONTEXT_H
#define CUDA_CONTEXT_H

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#ifdef __CUDACC__
#define CUDA_DEV __device__
#else
#define CUDA_DEV
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

// Extra profiling stream synchronizations?
//#define SSYNC HANDLE_ERROR( cudaStreamSynchronize(stream) )
#define SSYNC

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "definitions.h"
#include "include/splitvector/splitvec.h"

#include <stdio.h>

static const double BLOCK_ALLOCATION_PADDING = 2.5;
static const double BLOCK_ALLOCATION_FACTOR = 1.8;

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

class Managed {
public:
   void *operator new(size_t len) {
      void *ptr;
      cudaMallocManaged(&ptr, len);
      cudaDeviceSynchronize();
      return ptr;
   }

   void operator delete(void *ptr) {
      cudaDeviceSynchronize();
      cudaFree(ptr);
   }

   void* operator new[] (size_t len) {
      void *ptr;
      cudaMallocManaged(&ptr, len);
      cudaDeviceSynchronize();
      return ptr;
   }

   void operator delete[] (void* ptr) {
      cudaDeviceSynchronize();
      cudaFree(ptr);
   }

};

#define DIMS 1
#ifndef CUDABLOCKS
#  define CUDABLOCKS (108)
#endif
#ifndef CUDATHREADS
#  define CUDATHREADS (32) // NVIDIA: 32 AMD: 64
#endif

void cuda_set_device();
void cuda_clear_device();
cudaStream_t cuda_getStream();
cudaStream_t cuda_getPriorityStream();
int cuda_getDevice();

#define MAXCPUTHREADS 64

//extern CUcontext cuda_acc_context;
extern cudaStream_t cudaStreamList[];
extern cudaStream_t cudaPriorityStreamList[];
extern cudaStream_t cudaBaseStream;
extern Real *returnReal[];
extern Realf *returnRealf[];
extern vmesh::LocalID *returnLID[];

// Extern flag for stream attaching
extern bool needAttachedStreams;
extern bool doPrefetches;

#endif
