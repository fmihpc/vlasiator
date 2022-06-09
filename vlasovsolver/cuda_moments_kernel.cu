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
#include "../fieldsolver/fs_common.h" // divideIfNonZero()

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

using namespace std;

Realf *dev_meshDataPointers[MAXCPUTHREADS];
Real *dev_parameterPointers[MAXCPUTHREADS];
Real *dev_masses[MAXCPUTHREADS];
Real *dev_charges[MAXCPUTHREADS];
uint *dev_blockCounts[MAXCPUTHREADS];
Real *dev_momentArrays[MAXCPUTHREADS];

std::vector<Realf*> *meshDataPointers[MAXCPUTHREADS];
std::vector<Real*> *parameterPointers[MAXCPUTHREADS];
std::vector<Real> *masses[MAXCPUTHREADS];
std::vector<Real> *charges[MAXCPUTHREADS];
std::vector<uint> *blockCounts[MAXCPUTHREADS];
std::vector<std::array<Real,nMoments> > *momentArrays[MAXCPUTHREADS];

isCudaMomentsAllocated = false;

void cuda_allocateMomentCalculations() {
   if (isCudaMomentsAllocated) return;
   const uint nPopulations = getObjectWrapper().particleSpecies.size();
   const uint maxThreads = omp_get_max_threads();
   for (uint cpuThreadID=0; cpuThreadID<maxThreads; ++cpuThreadID) {
   
      //cudaMalloc
      HANDLE_ERROR( cudaMalloc((void**)dev_meshDataPointers[cpuThreadID], nPopulations*sizeof(uint64_t)) );
      HANDLE_ERROR( cudaMalloc((void**)dev_parameterPointers[cpuThreadID], nPopulations*sizeof(uint64_t)) );
      HANDLE_ERROR( cudaMalloc((void**)dev_blockCounts[cpuThreadID], nPopulations*sizeof(uint)) );
      HANDLE_ERROR( cudaMalloc((void**)dev_masses[cpuThreadID], nPopulations*sizeof(Real)) );
      HANDLE_ERROR( cudaMalloc((void**)dev_charges[cpuThreadID], nPopulations*sizeof(Real)) );
      HANDLE_ERROR( cudaMalloc((void**)dev_momentArrays[cpuThreadID], nMoments*(nPopulations+1)*sizeof(Real)) );
   
      // Also allocate and pin memory on host for faster transfers
      HANDLE_ERROR( cudaHostAlloc((void**)&meshDataPointers[cpuThreadID], nPopulations*sizeof(uint64_t), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&parameterPointers[cpuThreadID], nPopulations*sizeof(uint64_t), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&blockCounts[cpuThreadID], nPopulations*sizeof(uint), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&masses[cpuThreadID], nPopulations*sizeof(Real), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&charges[cpuThreadID], nPopulations*sizeof(Real), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&momentArrays[cpuThreadID], nMoments*(nPopulations+1)*sizeof(Real), cudaHostAllocPortable) );
   }
   isCudaMomentsAllocated = true;
   return;
}


// Define kernel
__global__ void moments_kernel(
   Realf* dev_meshDataPointers,
   Real* dev_parameterPointers,
   uint* dev_blockCounts,
   Real* dev_masses,
   Real* dev_charges,
   Real* dev_momentArrays,
   const int nPopulations,
   const bool computeSecond
   ){

   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = cellIndex(i,j,k);
   
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

      const uint nBlocks = dev_blockCounts[popID];
      const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];

      for (uint blockLID=0; blockLID<nBlocks; ++blockLID) {
         Real* blockParams = dev_parameterPointers[popID] + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;;
         Realf* avgs = dev_meshDataPointers[popID] + blockLID*WID3;
         const REAL VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
         const REAL VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
         const REAL VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
         n_sum[ti]   += avgs[ti];
         nvx_sum[ti] += avgs[ti]*VX;
         nvy_sum[ti] += avgs[ti]*VY;
         nvz_sum[ti] += avgs[ti]*VZ;
      }

      // Implemented just a simple non-optimized thread sum
      Real n = 0;
      Real nvx = 0;
      Real nvy = 0;
      Real nvz = 0;
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            n[ti] += n_sum[ti + s];
            nvx[ti] += nvx_sum[ti + s];
            nvy[ti] += nvy_sum[ti + s];
            nvz[ti] += nvz_sum[ti + s];
         }
         __syncwarp();
      }

      dev_momentArrays[popID*nMoments + 0] = n   * DV3;
      dev_momentArrays[popID*nMoments + 1] = nvx * DV3;
      dev_momentArrays[popID*nMoments + 2] = nvy * DV3;
      dev_momentArrays[popID*nMoments + 3] = nvz * DV3;

      // Sum over all populations
      dev_momentArrays[nPopulations*nMoments + 0] += n   * DV3 * masses[popID];
      dev_momentArrays[nPopulations*nMoments + 1] += nvx * DV3 * masses[popID];
      dev_momentArrays[nPopulations*nMoments + 2] += nvy * DV3 * masses[popID];
      dev_momentArrays[nPopulations*nMoments + 3] += nvz * DV3 * masses[popID];
      dev_momentArrays[nPopulations*nMoments + 4] += n   * DV3 * charges[popID];
      
   }
   dev_momentArrays[nPopulations*nMoments + 1] = divideIfNonZero(dev_momentArrays[nPopulations*nMoments + 1],dev_momentArrays[nPopulations*nMoments + 0]);
   dev_momentArrays[nPopulations*nMoments + 2] = divideIfNonZero(dev_momentArrays[nPopulations*nMoments + 2],dev_momentArrays[nPopulations*nMoments + 0]);
   dev_momentArrays[nPopulations*nMoments + 3] = divideIfNonZero(dev_momentArrays[nPopulations*nMoments + 3],dev_momentArrays[nPopulations*nMoments + 0]);

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
      
      const uint nBlocks = dev_blockCounts[popID];
      const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
      
      for (uint blockLID=0; blockLID<nBlocks; ++blockLID) {
         Real* blockParams = dev_parameterPointers[popID] + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;;
         Realf* avgs = dev_meshDataPointers[popID] + blockLID*WID3;

         const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
         const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
         const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
         
         nvx2_sum[ti] += avgs[ti] * (VX - averageVX) * (VX - averageVX);
         nvy2_sum[ti] += avgs[ti] * (VY - averageVY) * (VY - averageVY);
         nvz2_sum[ti] += avgs[ti] * (VZ - averageVZ) * (VZ - averageVZ);
      }

      // Implemented just a simple non-optimized thread sum
      Real nvx2 = 0;
      Real nvy2 = 0;
      Real nvz2 = 0;
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            nvx2[ti] += nvx2_sum[ti + s];
            nvy2[ti] += nvy2_sum[ti + s];
            nvz2[ti] += nvz2_sum[ti + s];
         }
         __syncwarp();
      }
      
      dev_momentArrays[popID*nMoments + 5] = nvx2 * DV3;
      dev_momentArrays[popID*nMoments + 6] = nvy2 * DV3;
      dev_momentArrays[popID*nMoments + 7] = nvz2 * DV3;

      // Sum over all populations
      dev_momentArrays[nPopulations*nMoments + 5] += nvx2 * DV3 * masses[popID];
      dev_momentArrays[nPopulations*nMoments + 6] += nvy2 * DV3 * masses[popID];
      dev_momentArrays[nPopulations*nMoments + 7] += nvz2 * DV3 * masses[popID];
   }
}

// Define kernel caller glue
void calculate_moments_glue(
   Realf* dev_meshDataPointers,
   Real* dev_parameterPointers,
   uint* dev_blockCounts,
   Real* dev_masses,
   Real* dev_charges,
   Real* dev_momentArrays,
   const int nPopulations,
   const bool computeSecond,
   cudaStream_t stream
   ) {
   dim3 block(WID,WID,WID);
   moments_kernel<<<1, block, 0, stream>>> (
      dev_meshDataPointers,
      dev_parameterPointers,
      dev_blockCounts,
      dev_masses,
      dev_charges,
      dev_momentArrays,
      nPopulations,
      computeSecond);
   return;
}
