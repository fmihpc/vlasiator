/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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
#include "gpu_moments.h"
#include "vlasovmover.h"
#include "../object_wrapper.h"
#include "../fieldsolver/fs_common.h" // divideIfNonZero(
#include "../arch/gpu_base.hpp"

using namespace std;

// Buffers for moment calculations
vmesh::VelocityBlockContainer** dev_VBC;
vmesh::VelocityBlockContainer** host_VBC;
Real* dev_moments1;
Real* dev_moments2;
Real* host_moments1;
Real* host_moments2;
uint gpu_allocated_moments = 0;

/**
    GPU kernel for calculating first velocity moments from provided
    velocity block containers
*/
__global__ void first_moments_kernel (
   vmesh::VelocityBlockContainer** dev_VBC,
   Real* dev_moments1,
   const uint nAllCells)
{
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z % WID;
   const int bi = threadIdx.z / WID;
   const int bincr = blockDim.z / WID;
   const int blockSize = blockDim.x*blockDim.y*blockDim.z;
   const uint celli = blockIdx.x; // userd for pointer to cell

   __shared__ Real smom[GPUTHREADS*WARPSPERBLOCK]; //==blockSize
   Real myMom[4] = {0};

   vmesh::VelocityBlockContainer* blockContainer = dev_VBC[celli];
   if (blockContainer==0) {
      if (ti==0) {
         for (uint imom=0; imom<4; imom++) {
            dev_moments1[celli*4 + imom] = 0;
         }
      }
      __syncthreads();
      return;
   }
   uint thisVBCSize = blockContainer->size();
   if (thisVBCSize==0) {
      if (ti==0) {
         for (uint imom=0; imom<4; imom++) {
            dev_moments1[celli*4 + imom] = 0;
         }
      }
      __syncthreads();
      return;
   }
   Realf *data = blockContainer->getData();
   Real *blockParameters = blockContainer->getParameters();
   const Real HALF = 0.5;
   const int indx = cellIndex(i,j,k);
   // this assumes that blockSize is divisible with WID3
   for (uint blockIndex = bi; blockIndex < thisVBCSize; blockIndex += bincr) {
       const Realf f = data[blockIndex*WID3+indx];
       const Real* blockParamsZ = &blockParameters[blockIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS];
       const Real DV3 = blockParamsZ[BlockParams::DVX]*blockParamsZ[BlockParams::DVY]*blockParamsZ[BlockParams::DVZ];
       const Real VX = blockParamsZ[BlockParams::VXCRD] + (i+HALF)*blockParamsZ[BlockParams::DVX];
       const Real VY = blockParamsZ[BlockParams::VYCRD] + (j+HALF)*blockParamsZ[BlockParams::DVY];
       const Real VZ = blockParamsZ[BlockParams::VZCRD] + (k+HALF)*blockParamsZ[BlockParams::DVZ];
       myMom[0] += f * DV3;
       myMom[1] += f * VX * DV3;
       myMom[2] += f * VY * DV3;
       myMom[3] += f * VZ * DV3;
   }
   __syncthreads();
   //Now reduce one-by-one for cell
   for (uint imom=0; imom<4; imom++) {
      smom[ti] = myMom[imom];
      for (unsigned int s=blockSize/2; s>0; s>>=1) {
         __syncthreads();
         if (ti < s) {
            smom[ti] += smom[ti + s];
         }
      }
      __syncthreads();
      if (ti==0) {
         dev_moments1[celli*4 + imom] = smom[0];
      }
      __syncthreads();
   }
 }

/**
    GPU kernel for calculating second velocity moments from provided
    velocity block containers
*/
__global__ void second_moments_kernel (
   vmesh::VelocityBlockContainer** dev_VBC,
   Real* dev_moments1,
   Real* dev_moments2,
   const uint nAllCells)
{
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z % WID;
   const int bi = threadIdx.z / WID;
   const int bincr = blockDim.z / WID;
   const int blockSize = blockDim.x*blockDim.y*blockDim.z;
   const uint celli = blockIdx.x; // userd for pointer to cell

   __shared__ Real smom[GPUTHREADS*WARPSPERBLOCK]; //==blockSize
   Real myMom[3] = {0};

   vmesh::VelocityBlockContainer* blockContainer = dev_VBC[celli];
   if (blockContainer==0) {
      if (ti==0) {
         for (uint imom=0; imom<4; imom++) {
            dev_moments1[celli*4 + imom] = 0;
         }
      }
      __syncthreads();
      return;
   }
   uint thisVBCSize = blockContainer->size();
   if (thisVBCSize==0) {
      if (ti==0) {
         for (uint imom=0; imom<4; imom++) {
            dev_moments1[celli*4 + imom] = 0;
         }
      }
      __syncthreads();
      return;
   }
   Realf *data = blockContainer->getData();
   Real *blockParameters = blockContainer->getParameters();

   const Real HALF = 0.5;
   const int indx = cellIndex(i,j,k);

   // index +0 is number density
   const Real averageVX = dev_moments1[celli*4 + 1];
   const Real averageVY = dev_moments1[celli*4 + 2];
   const Real averageVZ = dev_moments1[celli*4 + 3];

   for (uint blockIndex = bi; blockIndex < thisVBCSize; blockIndex += bincr) {
       const Realf f = data[blockIndex*WID3+indx];
       const Real* blockParamsZ = &blockParameters[blockIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS];
       const Real DV3 = blockParamsZ[BlockParams::DVX]*blockParamsZ[BlockParams::DVY]*blockParamsZ[BlockParams::DVZ];
       const Real VX = blockParamsZ[BlockParams::VXCRD] + (i+HALF)*blockParamsZ[BlockParams::DVX];
       const Real VY = blockParamsZ[BlockParams::VYCRD] + (j+HALF)*blockParamsZ[BlockParams::DVY];
       const Real VZ = blockParamsZ[BlockParams::VZCRD] + (k+HALF)*blockParamsZ[BlockParams::DVZ];
       myMom[0] += f * (VX - averageVX) * (VX - averageVX) * DV3;
       myMom[1] += f * (VY - averageVY) * (VY - averageVY) * DV3;
       myMom[2] += f * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
   }
   __syncthreads();
   //Now reduce one-by-one for cell
   for (uint imom=0; imom<3; imom++) {
      smom[ti] = myMom[imom];
      for (unsigned int s=blockSize/2; s>0; s>>=1) {
         __syncthreads();
         if (ti < s) {
            smom[ti] += smom[ti + s];
         }
      }
      __syncthreads();
      if (ti==0) {
         dev_moments2[celli*3 + imom] = smom[0];
      }
      __syncthreads();
   }
}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include contributions from
 * all existing particle populations. This function is AMR safe.
 * @param cell Spatial cell.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param doNotSkip If false, DO_NOT_COMPUTE cells are skipped.*/
// void calculateCellMoments(spatial_cell::SpatialCell* cell,
//                           const bool& computeSecond,
//                           const bool& computePopulationMomentsOnly,
//                           const bool& doNotSkip) {
// GPUTODO? For a single cell, might as well use ARCH.



/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _R variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void gpu_calculateMoments_R(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells_in,
   const bool& computeSecond) {

   phiprof::Timer computeMomentsTimer {"Compute _R moments"};

   // Ensure unique cells
   std::vector<CellID> cells = cells_in;
   std::sort( cells.begin(), cells.end() );
   cells.erase( std::unique( cells.begin(), cells.end() ), cells.end() );

   const uint nAllCells = cells.size();
   if (nAllCells==0) {
      return;
   }
   gpu_moments_allocate(nAllCells);

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Gather VBCs
      #pragma omp parallel for schedule(static)
      for(uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            host_VBC[celli] = 0;
            continue;
         }
         host_VBC[celli] = cell->dev_get_velocity_blocks(popID); // GPU-side VBC

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
      }

      CHK_ERR( gpuMemcpy(dev_VBC, host_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );
      // Launch kernel calculating this species' contribution to first velocity moments
      //const uint blocksPerBlock = GPUTHREADS * WARPSPERBLOCK / WID3;
      const uint blocksPerBlock = 1;
      dim3 launchSize(WID,WID,WID*blocksPerBlock);
      first_moments_kernel<<<nAllCells, launchSize, 0, 0>>> (
         dev_VBC,
         dev_moments1,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments1, dev_moments1, nAllCells*4*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real charge = getObjectWrapper().particleSpecies[popID].charge;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         pop.RHO_R = host_moments1[4*celli];
         pop.V_R[0] = divideIfNonZero(host_moments1[4*celli + 1], host_moments1[4*celli]);
         pop.V_R[1] = divideIfNonZero(host_moments1[4*celli + 2], host_moments1[4*celli]);
         pop.V_R[2] = divideIfNonZero(host_moments1[4*celli + 3], host_moments1[4*celli]);

         cell->parameters[CellParams::RHOM_R  ] += host_moments1[4*celli]*mass;
         cell->parameters[CellParams::VX_R] += host_moments1[4*celli + 1]*mass;
         cell->parameters[CellParams::VY_R] += host_moments1[4*celli + 2]*mass;
         cell->parameters[CellParams::VZ_R] += host_moments1[4*celli + 3]*mass;
         cell->parameters[CellParams::RHOQ_R  ] += host_moments1[4*celli]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

   phiprof::Timer computeMomentsVTimer {"compute-moments-R-cell-bulkV"};
   #pragma omp parallel for schedule(static)
   for (size_t celli=0; celli<nAllCells; ++celli) {
      SpatialCell* cell = mpiGrid[cells[celli]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      cell->parameters[CellParams::VX_R] = divideIfNonZero(cell->parameters[CellParams::VX_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VY_R] = divideIfNonZero(cell->parameters[CellParams::VY_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VZ_R] = divideIfNonZero(cell->parameters[CellParams::VZ_R], cell->parameters[CellParams::RHOM_R]);

      // copy the bulk flow frame back to device
      //host_moments1[4*celli + 0] = cell->parameters[CellParams::RHOM_R];
      host_moments1[4*celli + 1] = cell->parameters[CellParams::VX_R];
      host_moments1[4*celli + 2] = cell->parameters[CellParams::VY_R];
      host_moments1[4*celli + 3] = cell->parameters[CellParams::VZ_R];
   }
   computeMomentsVTimer.stop();

   // Compute second moments only if requested.
   if (computeSecond == false) {
      return;
   }

   CHK_ERR( gpuMemcpy(dev_moments1, host_moments1, nAllCells*4*sizeof(Real), gpuMemcpyHostToDevice) );
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Launch kernel calculating this species' contribution to second velocity moments
      //const uint blocksPerBlock = GPUTHREADS * WARPSPERBLOCK / WID3;
      const uint blocksPerBlock = 1;
      dim3 launchSize(WID,WID,WID*blocksPerBlock);
      second_moments_kernel<<<nAllCells, launchSize, 0, 0>>> (
         dev_VBC,
         dev_moments1,
         dev_moments2,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments2, dev_moments2, nAllCells*3*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         pop.P_R[0] = mass*host_moments2[3*celli + 0];
         pop.P_R[1] = mass*host_moments2[3*celli + 1];
         pop.P_R[2] = mass*host_moments2[3*celli + 2];

         cell->parameters[CellParams::P_11_R] += pop.P_R[0];
         cell->parameters[CellParams::P_22_R] += pop.P_R[1];
         cell->parameters[CellParams::P_33_R] += pop.P_R[2];
      } // for-loop over spatial cells
   } // for-loop over particle species

}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _V variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void gpu_calculateMoments_V(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells_in,
   const bool& computeSecond) {

   phiprof::Timer computeMomentsTimer {"Compute _V moments"};

   // Ensure unique cells
   std::vector<CellID> cells = cells_in;
   std::sort( cells.begin(), cells.end() );
   cells.erase( std::unique( cells.begin(), cells.end() ), cells.end() );

   const uint nAllCells = cells.size();
   if (nAllCells==0) {
      return;
   }
   gpu_moments_allocate(nAllCells);

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Gather VBCs
      #pragma omp parallel for schedule(static)
      for(uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            host_VBC[celli] = 0;
            continue;
         }
         host_VBC[celli] = cell->dev_get_velocity_blocks(popID); // GPU-side VBC

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
      }

      CHK_ERR( gpuMemcpy(dev_VBC, host_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );
      // Launch kernel calculating this species' contribution to first velocity moments
      //const uint blocksPerBlock = GPUTHREADS * WARPSPERBLOCK / WID3;
      const uint blocksPerBlock = 1;
      dim3 launchSize(WID,WID,WID*blocksPerBlock);
      first_moments_kernel<<<nAllCells, launchSize, 0, 0>>> (
         dev_VBC,
         dev_moments1,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments1, dev_moments1, nAllCells*4*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real charge = getObjectWrapper().particleSpecies[popID].charge;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         pop.RHO_V = host_moments1[4*celli];
         pop.V_V[0] = divideIfNonZero(host_moments1[4*celli + 1], host_moments1[4*celli]);
         pop.V_V[1] = divideIfNonZero(host_moments1[4*celli + 2], host_moments1[4*celli]);
         pop.V_V[2] = divideIfNonZero(host_moments1[4*celli + 3], host_moments1[4*celli]);

         cell->parameters[CellParams::RHOM_V  ] += host_moments1[4*celli]*mass;
         cell->parameters[CellParams::VX_V] += host_moments1[4*celli + 1]*mass;
         cell->parameters[CellParams::VY_V] += host_moments1[4*celli + 2]*mass;
         cell->parameters[CellParams::VZ_V] += host_moments1[4*celli + 3]*mass;
         cell->parameters[CellParams::RHOQ_V  ] += host_moments1[4*celli]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

   #pragma omp parallel for schedule(static)
   for (size_t celli=0; celli<nAllCells; ++celli) {
      phiprof::Timer computeMomentsCellTimer {"compute-moments-R-cell-bulkV"};
      SpatialCell* cell = mpiGrid[cells[celli]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      cell->parameters[CellParams::VX_V] = divideIfNonZero(cell->parameters[CellParams::VX_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VY_V] = divideIfNonZero(cell->parameters[CellParams::VY_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VZ_V] = divideIfNonZero(cell->parameters[CellParams::VZ_V], cell->parameters[CellParams::RHOM_V]);

      // copy the bulk flow frame back to device
      //host_moments1[4*celli + 0] = cell->parameters[CellParams::RHOM_V];
      host_moments1[4*celli + 1] = cell->parameters[CellParams::VX_V];
      host_moments1[4*celli + 2] = cell->parameters[CellParams::VY_V];
      host_moments1[4*celli + 3] = cell->parameters[CellParams::VZ_V];
   }

   // Compute second moments only if requested.
   if (computeSecond == false) {
      return;
   }

   CHK_ERR( gpuMemcpy(dev_moments1, host_moments1, nAllCells*4*sizeof(Real), gpuMemcpyHostToDevice) );
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Launch kernel calculating this species' contribution to second velocity moments
      //const uint blocksPerBlock = GPUTHREADS * WARPSPERBLOCK / WID3;
      const uint blocksPerBlock = 1;
      dim3 launchSize(WID,WID,WID*blocksPerBlock);
      second_moments_kernel<<<nAllCells, launchSize, 0, 0>>> (
         dev_VBC,
         dev_moments1,
         dev_moments2,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments2, dev_moments2, nAllCells*3*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         pop.P_V[0] = mass*host_moments2[3*celli + 0];
         pop.P_V[1] = mass*host_moments2[3*celli + 1];
         pop.P_V[2] = mass*host_moments2[3*celli + 2];

         cell->parameters[CellParams::P_11_V] += pop.P_V[0];
         cell->parameters[CellParams::P_22_V] += pop.P_V[1];
         cell->parameters[CellParams::P_33_V] += pop.P_V[2];
      } // for-loop over spatial cells
   } // for-loop over particle species

}


void gpu_moments_allocate(const uint nAllCells) {
   if (nAllCells <= gpu_allocated_moments) {
      return;
   }
   gpu_moments_deallocate();
   gpu_allocated_moments = nAllCells * BLOCK_ALLOCATION_FACTOR;

   // Host memory will be pinned
   CHK_ERR( gpuHostAlloc((void**)&host_VBC, gpu_allocated_moments*sizeof(vmesh::VelocityBlockContainer*)) );
   CHK_ERR( gpuHostAlloc((void**)&host_moments1, gpu_allocated_moments*4*sizeof(Real)) );
   CHK_ERR( gpuHostAlloc((void**)&host_moments2, gpu_allocated_moments*3*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_VBC, gpu_allocated_moments*sizeof(vmesh::VelocityBlockContainer*)) );
   CHK_ERR( gpuMalloc((void**)&dev_moments1, gpu_allocated_moments*4*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_moments2, gpu_allocated_moments*3*sizeof(Real)) );
   CHK_ERR( gpuDeviceSynchronize());
}

/* Deallocation at end of simulation */
void gpu_moments_deallocate() {
   if (gpu_allocated_moments==0) {
      return;
   }
   gpu_allocated_moments = 0;
   CHK_ERR( gpuFreeHost(host_VBC) );
   CHK_ERR( gpuFreeHost(host_moments1) );
   CHK_ERR( gpuFreeHost(host_moments2) );
   CHK_ERR( gpuFree(dev_VBC) );
   CHK_ERR( gpuFree(dev_moments1) );
   CHK_ERR( gpuFree(dev_moments2) );
}
