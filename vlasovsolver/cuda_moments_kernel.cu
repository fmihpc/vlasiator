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

#include "cuda_moments_kernel.cuh"
#include "../cuda_context.cuh"
#include "../common.h"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

using namespace std;

MomentInfo *dev_momentInfos[MAXCPUTHREADS];
MomentInfo *host_momentInfos[MAXCPUTHREADS];

Real *dev_momentArrays[MAXCPUTHREADS];
Real *host_momentArrays[MAXCPUTHREADS];

bool isCudaMomentsAllocated = false;

__host__ void cuda_allocateMomentCalculations(
   const uint nPopulations,
   const uint maxThreads
   ) {
   if (isCudaMomentsAllocated) return;
   for (uint cpuThreadID=0; cpuThreadID<maxThreads; ++cpuThreadID) {

      //cudaMalloc
      HANDLE_ERROR( cudaMalloc((void**)&dev_momentInfos[cpuThreadID], nPopulations*sizeof(MomentInfo)) );
      HANDLE_ERROR( cudaMalloc((void**)&dev_momentArrays[cpuThreadID], nMoments*(nPopulations+1)*sizeof(Real)) );

      // Also allocate and pin memory on host for faster transfers
      HANDLE_ERROR( cudaHostAlloc((void**)&host_momentInfos[cpuThreadID], nPopulations*sizeof(MomentInfo), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&host_momentArrays[cpuThreadID], nMoments*(nPopulations+1)*sizeof(Real), cudaHostAllocPortable) );
   }
   isCudaMomentsAllocated = true;
   return;
}


// Define kernel
__global__ void moments_kernel(
   MomentInfo *dev_momentInfos,
   Real* dev_momentArrays,
   const int nPopulations,
   const bool computeSecond
   ){

   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;

   const Real HALF = 0.5;
   __shared__ Real n_sum[WID3];
   __shared__ Real nvx_sum[WID3];
   __shared__ Real nvy_sum[WID3];
   __shared__ Real nvz_sum[WID3];

   for (uint popID=0; popID<nPopulations; ++popID) {
      n_sum[ti] = 0.0;
      nvx_sum[ti] = 0.0;
      nvy_sum[ti] = 0.0;
      nvz_sum[ti] = 0.0;

      const uint nBlocks = dev_momentInfos[popID].blockCount;
      const Real mass = dev_momentInfos[popID].mass;
      const Real charge = dev_momentInfos[popID].charge;
      Real* blockParams = dev_momentInfos[popID].parameterPointer;
      const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];

      for (uint blockLID=0; blockLID<nBlocks; ++blockLID) {
         Real* blockParams = dev_momentInfos[popID].parameterPointer + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
         Realf* avgs = dev_momentInfos[popID].meshDataPointer +blockLID*WID3;
         const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
         const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
         const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];

         n_sum[ti]   += avgs[ti];
         nvx_sum[ti] += avgs[ti]*VX;
         nvy_sum[ti] += avgs[ti]*VY;
         nvz_sum[ti] += avgs[ti]*VZ;
      }
      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            n_sum[ti] += n_sum[ti + s];
            nvx_sum[ti] += nvx_sum[ti + s];
            nvy_sum[ti] += nvy_sum[ti + s];
            nvz_sum[ti] += nvz_sum[ti + s];
         }
         __syncthreads();
      }
      if (ti==0) {
         dev_momentArrays[popID*nMoments + 0] = n_sum[0]   * DV3;
         dev_momentArrays[popID*nMoments + 1] = nvx_sum[0] * DV3;
         dev_momentArrays[popID*nMoments + 2] = nvy_sum[0] * DV3;
         dev_momentArrays[popID*nMoments + 3] = nvz_sum[0] * DV3;

         // Sum over all populations
         dev_momentArrays[nPopulations*nMoments + 0] += n_sum[0]   * DV3 * mass;
         dev_momentArrays[nPopulations*nMoments + 1] += nvx_sum[0] * DV3 * mass;
         dev_momentArrays[nPopulations*nMoments + 2] += nvy_sum[0] * DV3 * mass;
         dev_momentArrays[nPopulations*nMoments + 3] += nvz_sum[0] * DV3 * mass;
         dev_momentArrays[nPopulations*nMoments + 4] += n_sum[0]   * DV3 * charge;
      }
   }

   if (ti==0) {
      dev_momentArrays[nPopulations*nMoments + 1] = divideIfNonZero(dev_momentArrays[nPopulations*nMoments + 1],dev_momentArrays[nPopulations*nMoments + 0]);
      dev_momentArrays[nPopulations*nMoments + 2] = divideIfNonZero(dev_momentArrays[nPopulations*nMoments + 2],dev_momentArrays[nPopulations*nMoments + 0]);
      dev_momentArrays[nPopulations*nMoments + 3] = divideIfNonZero(dev_momentArrays[nPopulations*nMoments + 3],dev_momentArrays[nPopulations*nMoments + 0]);
   }

   if (computeSecond) {
      return;
   }
   // Continue to second moment calculation

   const Real averageVX = dev_momentArrays[nPopulations*nMoments + 1];
   const Real averageVY = dev_momentArrays[nPopulations*nMoments + 2];
   const Real averageVZ = dev_momentArrays[nPopulations*nMoments + 3];

   __shared__ Real nvx2_sum[WID3];
   __shared__ Real nvy2_sum[WID3];
   __shared__ Real nvz2_sum[WID3];

   for (uint popID=0; popID<nPopulations; ++popID) {

      nvx2_sum[ti] = 0.0;
      nvy2_sum[ti] = 0.0;
      nvz2_sum[ti] = 0.0;

      const uint nBlocks = dev_momentInfos[popID].blockCount;
      const Real mass = dev_momentInfos[popID].mass;
      //const Real charge = dev_momentInfos[popID].charge;
      Real* blockParams = dev_momentInfos[popID].parameterPointer;
      const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];

      for (uint blockLID=0; blockLID<nBlocks; ++blockLID) {
         Real* blockParams = dev_momentInfos[popID].parameterPointer + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
         Realf* avgs = dev_momentInfos[popID].meshDataPointer + blockLID*WID3;

         const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
         const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
         const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];

         nvx2_sum[ti] += avgs[ti] * (VX - averageVX) * (VX - averageVX);
         nvy2_sum[ti] += avgs[ti] * (VY - averageVY) * (VY - averageVY);
         nvz2_sum[ti] += avgs[ti] * (VZ - averageVZ) * (VZ - averageVZ);
      }

      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            nvx2_sum[ti] += nvx2_sum[ti + s];
            nvy2_sum[ti] += nvy2_sum[ti + s];
            nvz2_sum[ti] += nvz2_sum[ti + s];
         }
         __syncthreads();
      }
      if (ti==0) {
         dev_momentArrays[popID*nMoments + 5] = nvx2_sum[0] * DV3;
         dev_momentArrays[popID*nMoments + 6] = nvy2_sum[0] * DV3;
         dev_momentArrays[popID*nMoments + 7] = nvz2_sum[0] * DV3;

         // Sum over all populations
         dev_momentArrays[nPopulations*nMoments + 5] += nvx2_sum[0] * DV3 * mass;
         dev_momentArrays[nPopulations*nMoments + 6] += nvy2_sum[0] * DV3 * mass;
         dev_momentArrays[nPopulations*nMoments + 7] += nvz2_sum[0] * DV3 * mass;
      }
   }
}

// Define kernel caller glue
void calculate_moments_glue(
   MomentInfo *dev_momentInfos,
   Real* dev_momentArrays,
   const int nPopulations,
   const bool computeSecond,
   cudaStream_t stream
   ) {
   dim3 block(WID,WID,WID);
   moments_kernel<<<1, block, 0, stream>>> (
      dev_momentInfos,
      dev_momentArrays,
      nPopulations,
      computeSecond);
   return;
}
