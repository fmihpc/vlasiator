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

#include <phiprof.hpp>
#include "cuda_moments.h"
#include "../vlasovmover.h"
#include "../object_wrapper.h"

#include "../cuda_context.cuh"

using namespace std;

MomentInfo *dev_momentInfos[MAXCPUTHREADS];
MomentInfo *host_momentInfos[MAXCPUTHREADS];

Real *dev_momentArrays1[MAXCPUTHREADS];
Real *host_momentArrays1[MAXCPUTHREADS];
Real *dev_momentArrays2[MAXCPUTHREADS];
Real *host_momentArrays2[MAXCPUTHREADS];

bool isCudaMomentsAllocated = false;

__host__ void cuda_allocateMomentCalculations(
   const uint nPopulations,
   const uint maxThreads
   ) {
   if (isCudaMomentsAllocated) return;
   for (uint cpuThreadID=0; cpuThreadID<maxThreads; ++cpuThreadID) {

      //cudaMalloc
      HANDLE_ERROR( cudaMalloc((void**)&dev_momentInfos[cpuThreadID], nPopulations*sizeof(MomentInfo)) );
      HANDLE_ERROR( cudaMalloc((void**)&dev_momentArrays1[cpuThreadID], CUDABLOCKS*nMoments1*(nPopulations+1)*sizeof(Real)) );
      HANDLE_ERROR( cudaMalloc((void**)&dev_momentArrays2[cpuThreadID], CUDABLOCKS*nMoments2*(nPopulations+1)*sizeof(Real)) );

      // Also allocate and pin memory on host for faster transfers
      HANDLE_ERROR( cudaHostAlloc((void**)&host_momentInfos[cpuThreadID], nPopulations*sizeof(MomentInfo), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&host_momentArrays1[cpuThreadID], CUDABLOCKS*nMoments1*(nPopulations+1)*sizeof(Real), cudaHostAllocPortable) );
      HANDLE_ERROR( cudaHostAlloc((void**)&host_momentArrays2[cpuThreadID], CUDABLOCKS*nMoments2*(nPopulations+1)*sizeof(Real), cudaHostAllocPortable) );
   }
   isCudaMomentsAllocated = true;
   return;
}


// Define kernel for calculating zeroth and first velocity moments
__global__ void moments_first_kernel(
   MomentInfo *dev_momentInfos,
   Real* dev_momentArrays1,
   const int nPopulations
   ){

   const int cudaBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;

   const Real HALF = 0.5;
   __shared__ Real n_sum[WID3];
   __shared__ Real nvx_sum[WID3];
   __shared__ Real nvy_sum[WID3];
   __shared__ Real nvz_sum[WID3];

   const uint offset = blocki*nMoments1*(nPopulations+1);
   if (ti==0) {
      dev_momentArrays1[offset + nPopulations*nMoments1 + 0] = 0;
      dev_momentArrays1[offset + nPopulations*nMoments1 + 1] = 0;
      dev_momentArrays1[offset + nPopulations*nMoments1 + 2] = 0;
      dev_momentArrays1[offset + nPopulations*nMoments1 + 3] = 0;
      dev_momentArrays1[offset + nPopulations*nMoments1 + 4] = 0;
   }
   
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

      for (uint blockLID=blocki; blockLID<nBlocks; blockLID += cudaBlocks) {
         if (blockLID >= nBlocks) break;
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
         dev_momentArrays1[offset + popID*nMoments1 + 0] = n_sum[0]   * DV3;
         dev_momentArrays1[offset + popID*nMoments1 + 1] = nvx_sum[0] * DV3;
         dev_momentArrays1[offset + popID*nMoments1 + 2] = nvy_sum[0] * DV3;
         dev_momentArrays1[offset + popID*nMoments1 + 3] = nvz_sum[0] * DV3;

         // Sum over all populations
         dev_momentArrays1[offset + nPopulations*nMoments1 + 0] += n_sum[0]   * DV3 * mass;
         dev_momentArrays1[offset + nPopulations*nMoments1 + 1] += nvx_sum[0] * DV3 * mass;
         dev_momentArrays1[offset + nPopulations*nMoments1 + 2] += nvy_sum[0] * DV3 * mass;
         dev_momentArrays1[offset + nPopulations*nMoments1 + 3] += nvz_sum[0] * DV3 * mass;
         dev_momentArrays1[offset + nPopulations*nMoments1 + 4] += n_sum[0]   * DV3 * charge;
      }
   }
   return;
}

// Define kernel for calculating second velocity moments
__global__ void moments_second_kernel(
   MomentInfo *dev_momentInfos,
   Real* dev_momentArrays2,
   const int nPopulations,
   const Real bulkVX,
   const Real bulkVY,
   const Real bulkVZ
   ){

   const int cudaBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;

   const Real HALF = 0.5;
   __shared__ Real nvx2_sum[WID3];
   __shared__ Real nvy2_sum[WID3];
   __shared__ Real nvz2_sum[WID3];

   const uint offset = blocki*nMoments2*(nPopulations+1);
   if (ti==0) {
      dev_momentArrays2[offset + nPopulations*nMoments2 + 0] = 0;
      dev_momentArrays2[offset + nPopulations*nMoments2 + 1] = 0;
      dev_momentArrays2[offset + nPopulations*nMoments2 + 2] = 0;
   }

   for (uint popID=0; popID<nPopulations; ++popID) {
      nvx2_sum[ti] = 0.0;
      nvy2_sum[ti] = 0.0;
      nvz2_sum[ti] = 0.0;

      const uint nBlocks = dev_momentInfos[popID].blockCount;
      const Real mass = dev_momentInfos[popID].mass;
      //const Real charge = dev_momentInfos[popID].charge;
      Real* blockParams = dev_momentInfos[popID].parameterPointer;
      const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];

      for (uint blockLID=blocki; blockLID<nBlocks; blockLID += cudaBlocks) {
         if (blockLID >= nBlocks) break;
         Real* blockParams = dev_momentInfos[popID].parameterPointer + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
         Realf* avgs = dev_momentInfos[popID].meshDataPointer +blockLID*WID3;
         const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
         const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
         const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];

         nvx2_sum[ti] += avgs[ti] * (VX - bulkVX) * (VX - bulkVX);
         nvy2_sum[ti] += avgs[ti] * (VY - bulkVY) * (VY - bulkVY);
         nvz2_sum[ti] += avgs[ti] * (VZ - bulkVZ) * (VZ - bulkVZ);
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
         dev_momentArrays2[offset + popID*nMoments2 + 0] = nvx2_sum[0] * DV3 * mass;
         dev_momentArrays2[offset + popID*nMoments2 + 1] = nvy2_sum[0] * DV3 * mass;
         dev_momentArrays2[offset + popID*nMoments2 + 2] = nvz2_sum[0] * DV3 * mass;

         // Sum over all populations
         dev_momentArrays2[offset + nPopulations*nMoments2 + 0] += nvx2_sum[0] * DV3 * mass;
         dev_momentArrays2[offset + nPopulations*nMoments2 + 1] += nvy2_sum[0] * DV3 * mass;
         dev_momentArrays2[offset + nPopulations*nMoments2 + 2] += nvz2_sum[0] * DV3 * mass;
      }
   }
   return;
}


/** Calls a CUDA kernel for calculating zeroth, first, and (possibly) second
 * bulk velocity moments for the  given spatial cells.
 * Additionally, for each species, calculate the maximum
 * spatial time step so that CFL(spatial)=1. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _V variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void cuda_calculateMoments_V(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond) {

   phiprof::start("CUDA Compute _V moments");

   // Loop over all particle species
#pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const uint thread_id = omp_get_thread_num();
      SpatialCell* cell = mpiGrid[cells[c]];

      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }

      uint nPopulations = getObjectWrapper().particleSpecies.size();
      // Gather values and pointers for each population
      for (uint popID=0; popID<nPopulations; ++popID) {
         // Clear old moments to zero value
         if (popID == 0) {
            cell->parameters[CellParams::RHOM_V  ] = 0.0;
            cell->parameters[CellParams::VX_V] = 0.0;
            cell->parameters[CellParams::VY_V] = 0.0;
            cell->parameters[CellParams::VZ_V] = 0.0;
            cell->parameters[CellParams::RHOQ_V  ] = 0.0;
            cell->parameters[CellParams::P_11_V] = 0.0;
            cell->parameters[CellParams::P_22_V] = 0.0;
            cell->parameters[CellParams::P_33_V] = 0.0;
         }

         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         const uint nBlocks = vmesh->size();
         if (nBlocks == 0) continue;

         host_momentInfos[thread_id][popID].mass = getObjectWrapper().particleSpecies[popID].mass;
         host_momentInfos[thread_id][popID].charge = getObjectWrapper().particleSpecies[popID].charge;
         host_momentInfos[thread_id][popID].blockCount = nBlocks;
         host_momentInfos[thread_id][popID].parameterPointer = blockContainer->getParameters();
         host_momentInfos[thread_id][popID].meshDataPointer = blockContainer->getData();

         // For now: Launch cuda transfers as data isn't yet fully resident
         phiprof::start("CUDA-HtoD");
         blockContainer->dev_prefetchDevice();
         phiprof::stop("CUDA-HtoD");
      }

      // Transfer metadata to device
      phiprof::start("CUDA-HtoD");
      HANDLE_ERROR( cudaMemcpyAsync(dev_momentInfos[thread_id], host_momentInfos[thread_id], nPopulations*sizeof(MomentInfo), cudaMemcpyHostToDevice, cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-HtoD");

      // Now launch kernel for this spatial cell, all populations, zeroth and first moments
      phiprof::start("CUDA-firstMoments");
      dim3 block(WID,WID,WID);
      moments_first_kernel<<<CUDABLOCKS, block, 4*WID3*sizeof(Real), cudaStreamList[thread_id]>>> (
         dev_momentInfos[thread_id],
         dev_momentArrays1[thread_id],
         nPopulations);
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-firstMoments");

      // Transfer momentArrays1 back
      phiprof::start("CUDA-DtoH");
      HANDLE_ERROR( cudaMemcpyAsync(host_momentArrays1[thread_id], dev_momentArrays1[thread_id], CUDABLOCKS*nMoments1*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]) );
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-DtoH");

      for (uint popID=0; popID<nPopulations; ++popID) {
         Population & pop = cell->get_population(popID);
         pop.RHO_V = 0;
         pop.V_V[0] = 0;
         pop.V_V[1] = 0;
         pop.V_V[2] = 0;
         //pop.RHOQ_V = 0;
         for (uint iCuda=0; iCuda<CUDABLOCKS; ++iCuda) {
            const uint offset = iCuda*nMoments1*(nPopulations+1);
            pop.RHO_V  += host_momentArrays1[thread_id][offset + popID*nMoments1 + 0];
            pop.V_V[0] += host_momentArrays1[thread_id][offset + popID*nMoments1 + 1];
            pop.V_V[1] += host_momentArrays1[thread_id][offset + popID*nMoments1 + 2];
            pop.V_V[2] += host_momentArrays1[thread_id][offset + popID*nMoments1 + 3];
            // pop.RHOQ_V += host_momentArrays1[thread_id][offset + popID*nMoments1 + 4];
            if (popID==0) {
               cell->parameters[CellParams::RHOM] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 0];
               cell->parameters[CellParams::VX  ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 1];
               cell->parameters[CellParams::VY  ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 2];
               cell->parameters[CellParams::VZ  ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 3];
               cell->parameters[CellParams::RHOQ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 4];
            }
         }
         pop.V_V[0] = divideIfNonZero(pop.V_V[0], pop.RHO_V);
         pop.V_V[1] = divideIfNonZero(pop.V_V[1], pop.RHO_V);
         pop.V_V[2] = divideIfNonZero(pop.V_V[2], pop.RHO_V);
      }

      cell->parameters[CellParams::VX] = divideIfNonZero(cell->parameters[CellParams::VX], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VY] = divideIfNonZero(cell->parameters[CellParams::VY], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VZ] = divideIfNonZero(cell->parameters[CellParams::VZ], cell->parameters[CellParams::RHOM]);

      if (!computeSecond) {
         continue;
      }

      phiprof::start("CUDA-secondMoments");
      // Now launch kernel for this spatial cell, all populations, second moments
      moments_second_kernel<<<CUDABLOCKS, block, 3*WID3*sizeof(Real), cudaStreamList[thread_id]>>> (
         dev_momentInfos[thread_id],
         dev_momentArrays2[thread_id],
         nPopulations,
         cell->parameters[CellParams::VX],
         cell->parameters[CellParams::VY],
         cell->parameters[CellParams::VZ]);
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-secondMoments");

// Transfer momentArrays2 back
      phiprof::start("CUDA-DtoH");
      HANDLE_ERROR( cudaMemcpyAsync(host_momentArrays2[thread_id], dev_momentArrays2[thread_id], nMoments2*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]) );
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-DtoH");

      for (uint popID=0; popID<nPopulations; ++popID) {
         Population & pop = cell->get_population(popID);
         pop.P_V[0] = 0;
         pop.P_V[1] = 0;
         pop.P_V[2] = 0;
         for (uint iCuda=0; iCuda<CUDABLOCKS; ++iCuda) {
            const uint offset = iCuda*nMoments2*(nPopulations+1);
            pop.P_V[0] += host_momentArrays2[thread_id][offset + popID*nMoments2 + 0];
            pop.P_V[1] += host_momentArrays2[thread_id][offset + popID*nMoments2 + 1];
            pop.P_V[2] += host_momentArrays2[thread_id][offset + popID*nMoments2 + 2];
            if (popID==0) {
               cell->parameters[CellParams::P_11_V] = host_momentArrays2[thread_id][offset + nPopulations*nMoments2 + 0];
               cell->parameters[CellParams::P_22_V] = host_momentArrays2[thread_id][offset + nPopulations*nMoments2 + 1];
               cell->parameters[CellParams::P_33_V] = host_momentArrays2[thread_id][offset + nPopulations*nMoments2 + 2];
            }
         }
      }
   } // for-loop over spatial cells

   phiprof::stop("CUDA Compute _V moments");
   return;
}

/** Calls a CUDA kernel for calculating zeroth, first, and (possibly) second
 * bulk velocity moments for the  given spatial cells.
 * Additionally, for each species, calculate the maximum
 * spatial time step so that CFL(spatial)=1. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _R variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void cuda_calculateMoments_R(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond) {

   phiprof::start("CUDA Compute _R moments");

   // Loop over all particle species
#pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const uint thread_id = omp_get_thread_num();
      SpatialCell* cell = mpiGrid[cells[c]];

      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }

      uint nPopulations = getObjectWrapper().particleSpecies.size();
      // Gather values and pointers for each population
      for (uint popID=0; popID<nPopulations; ++popID) {
         // Clear old moments to zero value
         if (popID == 0) {
            cell->parameters[CellParams::RHOM_R  ] = 0.0;
            cell->parameters[CellParams::VX_R] = 0.0;
            cell->parameters[CellParams::VY_R] = 0.0;
            cell->parameters[CellParams::VZ_R] = 0.0;
            cell->parameters[CellParams::RHOQ_R  ] = 0.0;
            cell->parameters[CellParams::P_11_R] = 0.0;
            cell->parameters[CellParams::P_22_R] = 0.0;
            cell->parameters[CellParams::P_33_R] = 0.0;
         }

         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         const uint nBlocks = vmesh->size();
         if (nBlocks == 0) continue;

         host_momentInfos[thread_id][popID].mass = getObjectWrapper().particleSpecies[popID].mass;
         host_momentInfos[thread_id][popID].charge = getObjectWrapper().particleSpecies[popID].charge;
         host_momentInfos[thread_id][popID].blockCount = nBlocks;
         host_momentInfos[thread_id][popID].parameterPointer = blockContainer->getParameters();
         host_momentInfos[thread_id][popID].meshDataPointer = blockContainer->getData();

         // For now: Launch cuda transfers as data isn't yet fully resident
         phiprof::start("CUDA-HtoD");
         blockContainer->dev_prefetchDevice();
         phiprof::stop("CUDA-HtoD");
      }

      // Transfer metadata to device
      phiprof::start("CUDA-HtoD");
      HANDLE_ERROR( cudaMemcpyAsync(dev_momentInfos[thread_id], host_momentInfos[thread_id], nPopulations*sizeof(MomentInfo), cudaMemcpyHostToDevice, cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-HtoD");

      // Now launch kernel for this spatial cell, all populations, zeroth and first moments
      phiprof::start("CUDA-firstMoments");
      dim3 block(WID,WID,WID);
      moments_first_kernel<<<CUDABLOCKS, block, 4*WID3*sizeof(Real), cudaStreamList[thread_id]>>> (
         dev_momentInfos[thread_id],
         dev_momentArrays1[thread_id],
         nPopulations);
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-firstMoments");

      // Transfer momentArrays1 back
      phiprof::start("CUDA-DtoH");
      HANDLE_ERROR( cudaMemcpyAsync(host_momentArrays1[thread_id], dev_momentArrays1[thread_id], CUDABLOCKS*nMoments1*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]) );
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-DtoH");

      for (uint popID=0; popID<nPopulations; ++popID) {
         Population & pop = cell->get_population(popID);
         pop.RHO_R = 0;
         pop.V_R[0] = 0;
         pop.V_R[1] = 0;
         pop.V_R[2] = 0;
         //pop.RHOQ_R = 0;
         for (uint iCuda=0; iCuda<CUDABLOCKS; ++iCuda) {
            const uint offset = iCuda*nMoments1*(nPopulations+1);
            pop.RHO_R  += host_momentArrays1[thread_id][offset + popID*nMoments1 + 0];
            pop.V_R[0] += host_momentArrays1[thread_id][offset + popID*nMoments1 + 1];
            pop.V_R[1] += host_momentArrays1[thread_id][offset + popID*nMoments1 + 2];
            pop.V_R[2] += host_momentArrays1[thread_id][offset + popID*nMoments1 + 3];
            // pop.RHOQ_R += host_momentArrays1[thread_id][offset + popID*nMoments1 + 4];
            if (popID==0) {
               cell->parameters[CellParams::RHOM] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 0];
               cell->parameters[CellParams::VX  ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 1];
               cell->parameters[CellParams::VY  ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 2];
               cell->parameters[CellParams::VZ  ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 3];
               cell->parameters[CellParams::RHOQ] += host_momentArrays1[thread_id][offset + nPopulations*nMoments1 + 4];
            }
         }
         pop.V_R[0] = divideIfNonZero(pop.V_R[0], pop.RHO_R);
         pop.V_R[1] = divideIfNonZero(pop.V_R[1], pop.RHO_R);
         pop.V_R[2] = divideIfNonZero(pop.V_R[2], pop.RHO_R);
      }

      cell->parameters[CellParams::VX] = divideIfNonZero(cell->parameters[CellParams::VX], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VY] = divideIfNonZero(cell->parameters[CellParams::VY], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VZ] = divideIfNonZero(cell->parameters[CellParams::VZ], cell->parameters[CellParams::RHOM]);

      if (!computeSecond) {
         continue;
      }

      phiprof::start("CUDA-secondMoments");
      // Now launch kernel for this spatial cell, all populations, second moments
      moments_second_kernel<<<CUDABLOCKS, block, 3*WID3*sizeof(Real), cudaStreamList[thread_id]>>> (
         dev_momentInfos[thread_id],
         dev_momentArrays2[thread_id],
         nPopulations,
         cell->parameters[CellParams::VX],
         cell->parameters[CellParams::VY],
         cell->parameters[CellParams::VZ]);
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-secondMoments");

// Transfer momentArrays2 back
      phiprof::start("CUDA-DtoH");
      HANDLE_ERROR( cudaMemcpyAsync(host_momentArrays2[thread_id], dev_momentArrays2[thread_id], nMoments2*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]) );
      HANDLE_ERROR( cudaStreamSynchronize(cudaStreamList[thread_id]) );
      phiprof::stop("CUDA-DtoH");

      for (uint popID=0; popID<nPopulations; ++popID) {
         Population & pop = cell->get_population(popID);
         pop.P_R[0] = 0;
         pop.P_R[1] = 0;
         pop.P_R[2] = 0;
         for (uint iCuda=0; iCuda<CUDABLOCKS; ++iCuda) {
            const uint offset = iCuda*nMoments2*(nPopulations+1);
            pop.P_R[0] += host_momentArrays2[thread_id][offset + popID*nMoments2 + 0];
            pop.P_R[1] += host_momentArrays2[thread_id][offset + popID*nMoments2 + 1];
            pop.P_R[2] += host_momentArrays2[thread_id][offset + popID*nMoments2 + 2];
            if (popID==0) {
               cell->parameters[CellParams::P_11_R] = host_momentArrays2[thread_id][offset + nPopulations*nMoments2 + 0];
               cell->parameters[CellParams::P_22_R] = host_momentArrays2[thread_id][offset + nPopulations*nMoments2 + 1];
               cell->parameters[CellParams::P_33_R] = host_momentArrays2[thread_id][offset + nPopulations*nMoments2 + 2];
            }
         }
      }
   } // for-loop over spatial cells

   phiprof::stop("CUDA Compute _R moments");
   return;
}
