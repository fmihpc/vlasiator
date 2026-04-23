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

/**
    GPU kernel for calculating first velocity moments from provided
    velocity block containers
*/
__global__ void __launch_bounds__(WID3) first_moments_kernel (
   const vmesh::VelocityBlockContainer* __restrict__ const *dev_VBC,
   Real* dev_moments1,
   const uint nAllCells)
{
   const uint celli = blockIdx.x; // used for pointer to cell
   const uint stride = gridDim.y; // used for faster looping over contents
   const uint strideOffset = blockIdx.y; // used for faster looping over contents
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const int blockSize = blockDim.x*blockDim.y*blockDim.z;

   extern __shared__ Real smom[];
   Real myMom[nMom1] = {0};

   const vmesh::VelocityBlockContainer* __restrict__ blockContainer = dev_VBC[celli];
   if (blockContainer==0) {
      return;
   }
   const uint thisVBCSize = blockContainer->size();
   if (thisVBCSize==0) {
      return;
   }
   const Realf* __restrict__ data = blockContainer->getData();
   const Real* __restrict__ blockParameters = blockContainer->getParameters();
   const Real HALF = 0.5;
   const int indx = ti%WID3;
   for (uint blockIndex = strideOffset; blockIndex < thisVBCSize; blockIndex += stride) {
      const Realf cellValue = data[blockIndex*WID3+indx];
      const Real* __restrict__ blockParamsZ = &blockParameters[blockIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS];

      const Real paramX = blockParamsZ[BlockParams::DVX];
      const Real paramY = blockParamsZ[BlockParams::DVY];
      const Real paramZ = blockParamsZ[BlockParams::DVZ];

      const Real DV3 = paramX*paramY*paramZ;

      const Real VX = blockParamsZ[BlockParams::VXCRD] + (i+HALF)*paramX;
      const Real VY = blockParamsZ[BlockParams::VYCRD] + (j+HALF)*paramY;
      const Real VZ = blockParamsZ[BlockParams::VZCRD] + (k+HALF)*paramZ;
      
      const Real f = (Real)cellValue;
      const Real fDV3 = f * DV3;

      myMom[0] += fDV3;
      myMom[1] += fDV3 * VX;
      myMom[2] += fDV3 * VY;
      myMom[3] += fDV3 * VZ;
   }

   const int indexInsideWarp = ti % GPUTHREADS;
   const int warpIndex = ti / GPUTHREADS;

   //Now reduce one-by-one for cell
   for (int offset = GPUTHREADS/2; offset > 0; offset /= 2) {
      for (uint imom=0; imom<nMom1; imom++) {
         myMom[imom] += gpuKernelShflDown(myMom[imom], offset);
      }
   }

   if (indexInsideWarp == 0) {
      for (uint imom=0; imom<nMom1; imom++) {
         smom[warpIndex*nMom1+imom] = myMom[imom];
      }
   }

   __syncthreads();

   const int warpsPerBlock = blockSize/GPUTHREADS;

   for (uint imom=warpIndex; imom<nMom1; imom +=warpsPerBlock ) {
      myMom[imom] = (indexInsideWarp < warpsPerBlock) ? smom[indexInsideWarp*nMom1+imom] : 0.0;
      for (int offset = (warpsPerBlock)/2; offset > 0; offset /= 2) {
         myMom[imom] += gpuKernelShflDown(myMom[imom], offset);
      }

      if (indexInsideWarp == 0) {
         if (gridDim.y == 1) {
            dev_moments1[celli*nMom1 + imom] = myMom[imom];
         } else {
            atomicAdd(&dev_moments1[celli*nMom1 + imom],myMom[imom]);
         }
      }
   }
 }

/**
    GPU kernel for calculating second velocity moments from provided
    velocity block containers
*/
__global__ void __launch_bounds__(WID3) second_moments_kernel (
   const vmesh::VelocityBlockContainer* __restrict__ const *dev_VBC,
   const Real* __restrict__ dev_moments1,
   Real* dev_moments2,
   const uint nAllCells)
{
   const uint celli = blockIdx.x; // used for pointer to cell
   const uint stride = gridDim.y; // used for faster looping over contents
   const uint strideOffset = blockIdx.y; // used for faster looping over contents
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const int blockSize = blockDim.x*blockDim.y*blockDim.z;

   Real myMom[nMom2] = {0};
   extern __shared__ Real smom[];

   const vmesh::VelocityBlockContainer* __restrict__ blockContainer = dev_VBC[celli];
   if (blockContainer==0) {
      return;
   }
   const uint thisVBCSize = blockContainer->size();
   if (thisVBCSize==0) {
      return;
   }
   const Realf* __restrict__ data = blockContainer->getData();
   const Real* __restrict__ blockParameters = blockContainer->getParameters();

   const Real HALF = 0.5;
   const int indx = ti%WID3;

   // index +0 is number density
   const Real averageVX = dev_moments1[celli*nMom1 + 1];
   const Real averageVY = dev_moments1[celli*nMom1 + 2];
   const Real averageVZ = dev_moments1[celli*nMom1 + 3];

   for (uint blockIndex = strideOffset; blockIndex < thisVBCSize; blockIndex += stride) {
      const Realf cellValue = data[blockIndex*WID3+indx];
      const Real* __restrict__ blockParamsZ = &blockParameters[blockIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS];

      const Real paramX = blockParamsZ[BlockParams::DVX];
      const Real paramY = blockParamsZ[BlockParams::DVY];
      const Real paramZ = blockParamsZ[BlockParams::DVZ];

      const Real DV3 = paramX*paramY*paramZ;

      const Real VX = blockParamsZ[BlockParams::VXCRD] + (i+HALF)*paramX;
      const Real VY = blockParamsZ[BlockParams::VYCRD] + (j+HALF)*paramY;
      const Real VZ = blockParamsZ[BlockParams::VZCRD] + (k+HALF)*paramZ;

      const Real VXDifference = VX - averageVX;
      const Real VYDifference = VY - averageVY;
      const Real VZDifference = VZ - averageVZ;
      
      const Real f = (Real)cellValue;
      const Real fDV3 = f * DV3;

      const Real fDV3VXDifference = fDV3 * VXDifference;
      const Real fDV3VYDifference = fDV3 * VYDifference;

      myMom[0] += fDV3VXDifference * VXDifference;
      myMom[1] += fDV3VYDifference * VYDifference;
      myMom[2] += fDV3 * VZDifference * VZDifference;
      myMom[3] += fDV3VYDifference * VZDifference;
      myMom[4] += fDV3VXDifference * VZDifference;
      myMom[5] += fDV3VXDifference * VYDifference;
   }

   const int indexInsideWarp = ti % GPUTHREADS;
   const int warpIndex = ti / GPUTHREADS;

   //Now reduce one-by-one for cell
   for (int offset = GPUTHREADS/2; offset > 0; offset /= 2) {
      for (uint imom=0; imom<nMom2; imom++) {
         myMom[imom] += gpuKernelShflDown(myMom[imom], offset);
      }
   }

   if (indexInsideWarp == 0) {
      for (uint imom=0; imom<nMom2; imom++) {
         smom[warpIndex*nMom2+imom] = myMom[imom];
      }
   }

   __syncthreads();

   const int warpsPerBlock = blockSize/GPUTHREADS;

   for (uint imom=warpIndex; imom<nMom2; imom +=warpsPerBlock ) {
      myMom[imom] = (indexInsideWarp < warpsPerBlock) ? smom[indexInsideWarp*nMom2+imom] : 0.0;
      for (int offset = (warpsPerBlock)/2; offset > 0; offset /= 2) {
         myMom[imom] += gpuKernelShflDown(myMom[imom], offset);
      }

      if (indexInsideWarp == 0) {
         if (gridDim.y == 1) {
            dev_moments2[celli*nMom2 + imom] = myMom[imom];
         } else {
            atomicAdd(&dev_moments2[celli*nMom2 + imom],myMom[imom]);
         }
      }
   }
}

/**
    Note: there is no single-cell GPU-only version of moments calculations,
    (calculateCellMoments) as that task is achieved through the
    ARCH-interface with no performance loss.
*/

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _R variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param initialCompute If true, force re-calculation of outflow L1 sysboundary cell moments.
  (otherwise skipped as their VDF contents are not kept up to date)
*/
void gpu_calculateMoments_R(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells_in,
   const bool computeSecond,
   const bool initialCompute) {

   phiprof::Timer computeMomentsTimer {"Compute _R moments"};

   // Ensure unique cells
   std::vector<CellID> cells = cells_in;
   std::sort( cells.begin(), cells.end() );
   cells.erase( std::unique( cells.begin(), cells.end() ), cells.end() );

   const uint nAllCells = cells.size();
   if (nAllCells==0) {
      return;
   }

   gpuMemoryManager.startSession(0,0);

   SESSION_HOST_ALLOCATE(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*));
   SESSION_HOST_ALLOCATE(gpuMemoryManager, Real, host_moments1, nAllCells*nMom1*sizeof(Real));
   SESSION_HOST_ALLOCATE(gpuMemoryManager, Real, host_moments2, nAllCells*nMom2*sizeof(Real));
   SESSION_ALLOCATE(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*));
   SESSION_ALLOCATE(gpuMemoryManager, Real, dev_moments1, nAllCells*nMom1*sizeof(Real));
   SESSION_ALLOCATE(gpuMemoryManager, Real, dev_moments2, nAllCells*nMom2*sizeof(Real));

   vmesh::VelocityBlockContainer** host_VBC = GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBC);
   Real* host_moments1 = GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, host_moments1);
   Real* host_moments2 = GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, host_moments2);
   vmesh::VelocityBlockContainer** dev_VBC = GET_SESSION_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBC);
   Real* dev_moments1 = GET_SESSION_POINTER(gpuMemoryManager, Real, dev_moments1);
   Real* dev_moments2 = GET_SESSION_POINTER(gpuMemoryManager, Real, dev_moments2);

   std::vector<vmesh::LocalID> maxVmeshSizes;

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Gather VBCs
      maxVmeshSizes.push_back(0);
      #pragma omp parallel
      {
         vmesh::LocalID threadMaxVmeshSize = 0;
         #pragma omp for schedule(static)
         for(uint celli = 0; celli < nAllCells; celli++){
            SpatialCell* cell = mpiGrid[cells[celli]];
            if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               host_VBC[celli] = 0;
               continue;
            }
            if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
               // these should have been handled by the boundary code
               host_VBC[celli] = 0;
               continue;
            }
            host_VBC[celli] = cell->dev_get_velocity_blocks(popID); // GPU-side VBC
            // Evaluate cached vmesh size
            const vmesh::LocalID meshSize = cell->get_velocity_mesh(popID)->size();
            threadMaxVmeshSize = meshSize > threadMaxVmeshSize ? meshSize : threadMaxVmeshSize;

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
               cell->parameters[CellParams::P_23_R] = 0.0;
               cell->parameters[CellParams::P_13_R] = 0.0;
               cell->parameters[CellParams::P_12_R] = 0.0;
            }
         }
         #pragma omp critical
         {
            maxVmeshSizes.at(popID) = maxVmeshSizes.at(popID) > threadMaxVmeshSize ? maxVmeshSizes.at(popID) : threadMaxVmeshSize;
         }
      }
      if (maxVmeshSizes.at(popID) == 0) {
         continue;
      }

      vmesh::LocalID maxVmeshLaunch = (gpuMultiProcessorCount*min(threadsPerMP/WID3, blocksPerMP))/nAllCells;
      maxVmeshLaunch = maxVmeshLaunch < 1 ? 1 : maxVmeshLaunch;
      // Send pointers, set initial data to zero
      CHK_ERR( gpuMemcpy(dev_VBC, host_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemset(dev_moments1, 0, nAllCells*nMom1*sizeof(Real)) );
      // Launch kernel calculating this species' contribution to first velocity moments
      dim3 blockSize(WID,WID,WID);
      dim3 gridSize(nAllCells,maxVmeshLaunch,1);
      int sharedMemorySizeFirstMoments = nMom1 * (WID3 / GPUTHREADS) * sizeof(Real);
      first_moments_kernel<<<gridSize, blockSize, sharedMemorySizeFirstMoments, 0>>> (
         dev_VBC,
         dev_moments1,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments1, dev_moments1, nAllCells*nMom1*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real charge = getObjectWrapper().particleSpecies[popID].charge;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
            // these should have been handled by the boundary code
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         pop.RHO_R = host_moments1[nMom1*celli];
         pop.V_R[0] = divideIfNonZero(host_moments1[nMom1*celli + 1], host_moments1[nMom1*celli]);
         pop.V_R[1] = divideIfNonZero(host_moments1[nMom1*celli + 2], host_moments1[nMom1*celli]);
         pop.V_R[2] = divideIfNonZero(host_moments1[nMom1*celli + 3], host_moments1[nMom1*celli]);

         cell->parameters[CellParams::RHOM_R  ] += host_moments1[nMom1*celli]*mass;
         cell->parameters[CellParams::VX_R] += host_moments1[nMom1*celli + 1]*mass;
         cell->parameters[CellParams::VY_R] += host_moments1[nMom1*celli + 2]*mass;
         cell->parameters[CellParams::VZ_R] += host_moments1[nMom1*celli + 3]*mass;
         cell->parameters[CellParams::RHOQ_R  ] += host_moments1[nMom1*celli]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

   #pragma omp parallel for schedule(static)
   for (size_t celli=0; celli<nAllCells; ++celli) {
      SpatialCell* cell = mpiGrid[cells[celli]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
         // these should have been handled by the boundary code
         continue;
      }
      cell->parameters[CellParams::VX_R] = divideIfNonZero(cell->parameters[CellParams::VX_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VY_R] = divideIfNonZero(cell->parameters[CellParams::VY_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VZ_R] = divideIfNonZero(cell->parameters[CellParams::VZ_R], cell->parameters[CellParams::RHOM_R]);

      // copy the bulk flow frame back to device
      //host_moments1[nMom1*celli + 0] = cell->parameters[CellParams::RHOM_R];
      host_moments1[nMom1*celli + 1] = cell->parameters[CellParams::VX_R];
      host_moments1[nMom1*celli + 2] = cell->parameters[CellParams::VY_R];
      host_moments1[nMom1*celli + 3] = cell->parameters[CellParams::VZ_R];
   }

   // Compute second moments only if requested.
   if (computeSecond == false) {
      gpuMemoryManager.endSession();
      return;
   }

   CHK_ERR( gpuMemcpy(dev_moments1, host_moments1, nAllCells*nMom1*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemset(dev_moments2, 0, nAllCells*nMom2*sizeof(Real)) );
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      if (maxVmeshSizes.at(popID) == 0) {
         continue;
      }
      // Launch kernel calculating this species' contribution to second velocity moments

      vmesh::LocalID maxVmeshLaunch = (gpuMultiProcessorCount*min(threadsPerMP/WID3, blocksPerMP))/nAllCells;
      maxVmeshLaunch = maxVmeshLaunch < 1 ? 1 : maxVmeshLaunch;

      dim3 blockSize(WID,WID,WID);
      dim3 gridSize(nAllCells,maxVmeshLaunch,1);
      int sharedMemorySizeSecondMoments = nMom2 * (WID3 / GPUTHREADS) * sizeof(Real);
      second_moments_kernel<<<gridSize, blockSize, sharedMemorySizeSecondMoments, 0>>> (
         dev_VBC,
         dev_moments1,
         dev_moments2,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments2, dev_moments2, nAllCells*nMom2*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
            // these should have been handled by the boundary code
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         for (size_t i = 0; i < nMom2; ++i) {
            pop.P_R[i] = mass*host_moments2[nMom2*celli + i];
         }

         cell->parameters[CellParams::P_11_R] += pop.P_R[0];
         cell->parameters[CellParams::P_22_R] += pop.P_R[1];
         cell->parameters[CellParams::P_33_R] += pop.P_R[2];
         cell->parameters[CellParams::P_23_R] += pop.P_R[3];
         cell->parameters[CellParams::P_13_R] += pop.P_R[4];
         cell->parameters[CellParams::P_12_R] += pop.P_R[5];
      } // for-loop over spatial cells
   } // for-loop over particle species

   gpuMemoryManager.endSession();
}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _V variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param initialCompute If true, force re-calculation of outflow L1 sysboundary cell moments.
  (otherwise skipped as their VDF contents are not kept up to date)
 */
void gpu_calculateMoments_V(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells_in,
   const bool computeSecond,
   const bool initialCompute) {

   phiprof::Timer computeMomentsTimer {"Compute _V moments"};

   // Ensure unique cells
   std::vector<CellID> cells = cells_in;
   std::sort( cells.begin(), cells.end() );
   cells.erase( std::unique( cells.begin(), cells.end() ), cells.end() );

   const uint nAllCells = cells.size();
   if (nAllCells==0) {
      return;
   }

   gpuMemoryManager.startSession(0,0);

   SESSION_HOST_ALLOCATE(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*));
   SESSION_HOST_ALLOCATE(gpuMemoryManager, Real, host_moments1, nAllCells*nMom1*sizeof(Real));
   SESSION_HOST_ALLOCATE(gpuMemoryManager, Real, host_moments2, nAllCells*nMom2*sizeof(Real));
   SESSION_ALLOCATE(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*));
   SESSION_ALLOCATE(gpuMemoryManager, Real, dev_moments1, nAllCells*nMom1*sizeof(Real));
   SESSION_ALLOCATE(gpuMemoryManager, Real, dev_moments2, nAllCells*nMom2*sizeof(Real));

   vmesh::VelocityBlockContainer** host_VBC = GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBC);
   Real* host_moments1 = GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, host_moments1);
   Real* host_moments2 = GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, host_moments2);
   vmesh::VelocityBlockContainer** dev_VBC = GET_SESSION_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBC);
   Real* dev_moments1 = GET_SESSION_POINTER(gpuMemoryManager, Real, dev_moments1);
   Real* dev_moments2 = GET_SESSION_POINTER(gpuMemoryManager, Real, dev_moments2);

   std::vector<vmesh::LocalID> maxVmeshSizes;

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Gather VBCs
      maxVmeshSizes.push_back(0);
      #pragma omp parallel
      {
         vmesh::LocalID threadMaxVmeshSize = 0;
         #pragma omp for schedule(static)
         for(uint celli = 0; celli < nAllCells; celli++){
            SpatialCell* cell = mpiGrid[cells[celli]];
            if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               host_VBC[celli] = 0;
               continue;
            }
            if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
               // these should have been handled by the boundary code
               host_VBC[celli] = 0;
               continue;
            }
            host_VBC[celli] = cell->dev_get_velocity_blocks(popID); // GPU-side VBC
            // Evaluate cached vmesh size
            const vmesh::LocalID meshSize = cell->get_velocity_mesh(popID)->size();
            threadMaxVmeshSize = meshSize > threadMaxVmeshSize ? meshSize : threadMaxVmeshSize;

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
               cell->parameters[CellParams::P_23_V] = 0.0;
               cell->parameters[CellParams::P_13_V] = 0.0;
               cell->parameters[CellParams::P_12_V] = 0.0;
            }
         }
#pragma omp critical
         {
            maxVmeshSizes.at(popID) = maxVmeshSizes.at(popID) > threadMaxVmeshSize ? maxVmeshSizes.at(popID) : threadMaxVmeshSize;
         }
      }
      if (maxVmeshSizes.at(popID) == 0) {
         continue;
      }

      vmesh::LocalID maxVmeshLaunch = (gpuMultiProcessorCount*min(threadsPerMP/WID3, blocksPerMP))/nAllCells;
      maxVmeshLaunch =  maxVmeshLaunch < 1 ? 1 : maxVmeshLaunch;
      // Send pointers, set initial data to zero
      CHK_ERR( gpuMemcpy(dev_VBC, host_VBC, nAllCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemset(dev_moments1, 0, nAllCells*nMom1*sizeof(Real)) );
      // Launch kernel calculating this species' contribution to first velocity moments
      dim3 blockSize(WID,WID,WID);
      dim3 gridSize(nAllCells,maxVmeshLaunch,1);
      int sharedMemorySizeFirstMoments = nMom1 * (WID3 / GPUTHREADS) * sizeof(Real);
      first_moments_kernel<<<gridSize, blockSize, sharedMemorySizeFirstMoments, 0>>> (
         dev_VBC,
         dev_moments1,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments1, dev_moments1, nAllCells*nMom1*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real charge = getObjectWrapper().particleSpecies[popID].charge;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
            // these should have been handled by the boundary code
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         pop.RHO_V = host_moments1[nMom1*celli];
         pop.V_V[0] = divideIfNonZero(host_moments1[nMom1*celli + 1], host_moments1[nMom1*celli]);
         pop.V_V[1] = divideIfNonZero(host_moments1[nMom1*celli + 2], host_moments1[nMom1*celli]);
         pop.V_V[2] = divideIfNonZero(host_moments1[nMom1*celli + 3], host_moments1[nMom1*celli]);

         cell->parameters[CellParams::RHOM_V  ] += host_moments1[nMom1*celli]*mass;
         cell->parameters[CellParams::VX_V] += host_moments1[nMom1*celli + 1]*mass;
         cell->parameters[CellParams::VY_V] += host_moments1[nMom1*celli + 2]*mass;
         cell->parameters[CellParams::VZ_V] += host_moments1[nMom1*celli + 3]*mass;
         cell->parameters[CellParams::RHOQ_V  ] += host_moments1[nMom1*celli]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

   #pragma omp parallel for schedule(static)
   for (size_t celli=0; celli<nAllCells; ++celli) {
      SpatialCell* cell = mpiGrid[cells[celli]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
         // these should have been handled by the boundary code
         continue;
      }
      cell->parameters[CellParams::VX_V] = divideIfNonZero(cell->parameters[CellParams::VX_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VY_V] = divideIfNonZero(cell->parameters[CellParams::VY_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VZ_V] = divideIfNonZero(cell->parameters[CellParams::VZ_V], cell->parameters[CellParams::RHOM_V]);

      // copy the bulk flow frame back to device
      //host_moments1[nMom1*celli + 0] = cell->parameters[CellParams::RHOM_V];
      host_moments1[nMom1*celli + 1] = cell->parameters[CellParams::VX_V];
      host_moments1[nMom1*celli + 2] = cell->parameters[CellParams::VY_V];
      host_moments1[nMom1*celli + 3] = cell->parameters[CellParams::VZ_V];
   }

   // Compute second moments only if requested.
   if (computeSecond == false) {
      gpuMemoryManager.endSession();
      return;
   }

   CHK_ERR( gpuMemcpy(dev_moments1, host_moments1, nAllCells*nMom1*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemset(dev_moments2, 0, nAllCells*nMom2*sizeof(Real)) );
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      if (maxVmeshSizes.at(popID) == 0) {
         continue;
      }

      vmesh::LocalID maxVmeshLaunch = (gpuMultiProcessorCount*min(threadsPerMP/WID3, blocksPerMP))/nAllCells;
      maxVmeshLaunch = maxVmeshLaunch < 1 ? 1 : maxVmeshLaunch;

      dim3 blockSize(WID,WID,WID);
      dim3 gridSize(nAllCells,maxVmeshLaunch,1);
      int sharedMemorySizeSecondMoments = nMom2 * (WID3 / GPUTHREADS) * sizeof(Real);
      second_moments_kernel<<<gridSize, blockSize, sharedMemorySizeSecondMoments, 0>>> (
         dev_VBC,
         dev_moments1,
         dev_moments2,
         nAllCells
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize());
      CHK_ERR( gpuMemcpy(host_moments2, dev_moments2, nAllCells*nMom2*sizeof(Real), gpuMemcpyDeviceToHost) );

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      #pragma omp parallel for schedule(static)
      for (uint celli = 0; celli < nAllCells; celli++){
         SpatialCell* cell = mpiGrid[cells[celli]];
         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) {
            // these should have been handled by the boundary code
            continue;
         }

         // Store species' contribution to bulk velocity moments
         Population &pop = cell->get_population(popID);
         for (size_t i = 0; i < nMom2; ++i) {
            pop.P_V[i] = mass*host_moments2[nMom2*celli + i];
         }

         cell->parameters[CellParams::P_11_V] += pop.P_V[0];
         cell->parameters[CellParams::P_22_V] += pop.P_V[1];
         cell->parameters[CellParams::P_33_V] += pop.P_V[2];
         cell->parameters[CellParams::P_23_V] += pop.P_V[3];
         cell->parameters[CellParams::P_13_V] += pop.P_V[4];
         cell->parameters[CellParams::P_12_V] += pop.P_V[5];
      } // for-loop over spatial cells
   } // for-loop over particle species

   gpuMemoryManager.endSession();
}
