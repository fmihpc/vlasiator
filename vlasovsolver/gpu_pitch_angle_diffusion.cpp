/*
 * This file is part of Vlasiator.
 * Copyright 2010-2025 Finnish Meteorological Institute and University of Helsinki
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


// This is the GPU version of cpu_pitch_angle_diffusion

#include "../parameters.h"
#include "../object_wrapper.h"
#include <math.h>
//#include <cmath> // NaN Inf checks
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <Eigen/Geometry>
#include "common_pitch_angle_diffusion.hpp"
#include "../arch/gpu_base.hpp"
#include "../spatial_cells/block_adjust_gpu.hpp"

#define GPUCELLMUSPACE(var,cellIdx,v_ind,mu_ind) var[(cellIdx)*nbins_v*nbins_mu+(mu_ind)*nbins_v + (v_ind)]

__global__ void __launch_bounds__(WID3) build2dArrayOfFvmu_kernel(
   size_t *dev_cellIdxArray,
   size_t *dev_velocityIdxArray,
   const vmesh::VelocityBlockContainer* __restrict__ const *dev_velocityBlockContainer,
   Real *dev_bulkVX,
   Real* dev_bulkVY,
   Real* dev_bulkVZ,
   Real *dev_bValues,
   Realf *dev_fmu,
   int *dev_fcount,
   Real dVbins,
   const Real dmubins,
   int nbins_v,
   int nbins_mu
   ){
   
   int totalBlockIndex = blockIdx.x; // Corresponds to index spatial and velocity blocks

   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   size_t cellIdx = dev_cellIdxArray[totalBlockIndex];
   size_t velocityIdx = dev_velocityIdxArray[totalBlockIndex];

   const Real* __restrict__ blockParameters = dev_velocityBlockContainer[cellIdx]->getParameters(velocityIdx);

   //Get velocity space coordinates
   const Real VX = blockParameters[BlockParams::VXCRD] + (i + 0.5)*blockParameters[BlockParams::DVX];
   const Real VY = blockParameters[BlockParams::VYCRD] + (j + 0.5)*blockParameters[BlockParams::DVY];
   const Real VZ = blockParameters[BlockParams::VZCRD] + (k + 0.5)*blockParameters[BlockParams::DVZ];

   const Real VplasmaX = VX - dev_bulkVX[cellIdx];
   const Real VplasmaY = VY - dev_bulkVY[cellIdx];
   const Real VplasmaZ = VZ - dev_bulkVZ[cellIdx];
   
   const Real normV = sqrt(VplasmaX*VplasmaX + VplasmaY*VplasmaY + VplasmaZ*VplasmaZ);
   const Real Vpara = VplasmaX*dev_bValues[3*cellIdx] + VplasmaY*dev_bValues[3*cellIdx+1] + VplasmaZ*dev_bValues[3*cellIdx+2];
   const Real mu = Vpara/(normV+std::numeric_limits<Real>::min()); // + min value to avoid division by 0.

   const int Vindex = static_cast<int>(std::nearbyint(floor((normV) / dVbins)));
   const Real Vmu = dVbins * (Vindex+0.5); // Take value at the center of the mu cell
   int muindex = static_cast<int>(std::nearbyint(floor((mu+1.0) / dmubins)));

   const Realf cellValue = dev_velocityBlockContainer[cellIdx]->getData()[velocityIdx*WID3+k*WID2+j*WID+i];
   const Real increment = 2.0 * M_PI * Vmu*Vmu * cellValue;
   // Safety check to handle edge case where mu = exactly 1.0
   const int mui = std::max(0,std::min(muindex,nbins_mu-1));
   const int vi = std::max(0,std::min(Vindex,nbins_v-1));

   // TODO: can this be done without atomicAdd while avoiding race conditions?
   atomicAdd(&GPUCELLMUSPACE(dev_fmu,cellIdx,vi,mui), increment);
   atomicAdd(&GPUCELLMUSPACE(dev_fcount,cellIdx,vi,mui), 1);
}

__global__ void __launch_bounds__(WID3) computeNewCellValues_kernel(
   size_t *dev_cellIdxArray,
   size_t *dev_remappedCellIdxArray,
   size_t *dev_velocityIdxArray,
   vmesh::VelocityBlockContainer* __restrict__ *dev_velocityBlockContainer,
   Real *dev_bulkVX,
   Real* dev_bulkVY,
   Real* dev_bulkVZ,
   Real *dev_bValues,
   Realf *dev_dfdt_mu,
   Real *dev_Ddt,
   Real dVbins,
   const Real dmubins,
   int nbins_v,
   int nbins_mu
   ){
   
   int totalBlockIndex = blockIdx.x; // Corresponds to index spatial and velocity blocks

   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   size_t cellIdx = dev_cellIdxArray[totalBlockIndex];
   size_t velocityIdx = dev_velocityIdxArray[totalBlockIndex];

   const Real* __restrict__ blockParameters = dev_velocityBlockContainer[cellIdx]->getParameters(velocityIdx);

   //Get velocity space coordinates
   const Real VX = blockParameters[BlockParams::VXCRD] + (i + 0.5)*blockParameters[BlockParams::DVX];
   const Real VY = blockParameters[BlockParams::VYCRD] + (j + 0.5)*blockParameters[BlockParams::DVY];
   const Real VZ = blockParameters[BlockParams::VZCRD] + (k + 0.5)*blockParameters[BlockParams::DVZ];

   const Real VplasmaX = VX - dev_bulkVX[cellIdx];
   const Real VplasmaY = VY - dev_bulkVY[cellIdx];
   const Real VplasmaZ = VZ - dev_bulkVZ[cellIdx];
   
   const Real normV = sqrt(VplasmaX*VplasmaX + VplasmaY*VplasmaY + VplasmaZ*VplasmaZ);
   const Real Vpara = VplasmaX*dev_bValues[3*cellIdx] + VplasmaY*dev_bValues[3*cellIdx+1] + VplasmaZ*dev_bValues[3*cellIdx+2];
   const Real mu = Vpara/(normV+std::numeric_limits<Real>::min()); // + min value to avoid division by 0.

   const int Vindex = static_cast<int>(std::nearbyint(floor((normV) / dVbins)));
   int muindex = static_cast<int>(std::nearbyint(floor((mu+1.0) / dmubins)));

   Realf dfdt = 0.0;
   // Safety check to handle edge case where mu = exactly 1.0
   const int mui = std::max(0,std::min(muindex,nbins_mu-1));
   const int vi = std::max(0,std::min(Vindex,nbins_v-1));
   dfdt = GPUCELLMUSPACE(dev_dfdt_mu,cellIdx,vi,mui); // dfdt_mu was scaled back down by 2pi*v^2 on creation

   // Update cell value, ensuring result is non-negative
   Realf NewCellValue    = (dev_velocityBlockContainer[cellIdx]->getData()[velocityIdx*WID3+k*WID2+j*WID+i]) + dfdt * dev_Ddt[cellIdx];
   const bool lessZero = (NewCellValue < 0.0);
   NewCellValue = lessZero ? 0.0 : NewCellValue;
   dev_velocityBlockContainer[cellIdx]->getData()[velocityIdx*WID3+k*WID2+j*WID+i] = NewCellValue;
}

__global__ void computeDerivativesCFLDdt_kernel(
   size_t *dev_smallCellIdxArray,
   int *dev_fcount,
   Realf *dev_fmu,
   Real *dev_nu0Values,
   Realf *dev_dfdt_mu,
   int *dev_cellIdxKeys,
   Real *dev_potentialDdtValues,
   Realf *dev_sparsity,
   Real dVbins,
   const Real dmubins,
   const Real epsilon,
   Realf PADCFL,
   int nbins_v,
   int nbins_mu,
   int blocksPerSpatialCell,
   int lastBlockSize
   ){

   int spatialBlockIndex = blockIdx.x/blocksPerSpatialCell; // Corresponds to index spatial and velocity blocks
   int indexInsideBlock = blockIdx.x%blocksPerSpatialCell;
   int nextIndexInsideBlock = (blockIdx.x+1)%blocksPerSpatialCell;
   int idx = indexInsideBlock*blockDim.x+threadIdx.x;
   int threadIndex = threadIdx.x;
   int indv = idx%nbins_v;
   int indmu = (idx/nbins_v)%nbins_mu;
   size_t cellIdx = dev_smallCellIdxArray[spatialBlockIndex];
   
   // Initiate shared memory values
   extern __shared__ Real localDdtValues[];
   localDdtValues[threadIndex] = std::numeric_limits<Real>::max();

   __syncthreads();

   if(nextIndexInsideBlock != 0 || threadIndex < lastBlockSize){

      // Search limits for how many cells in mu-direction should be max evaluated when searching for a near neighbour?
      // Assuming some oversampling; changing these values may result in method breaking at very small plasma frame velocities.
      const int rlimit = nbins_mu-1;
      const int llimit = 0;

      int cLeft = 0;
      int cRight = 0;

      // !!! DANGER, WARP DIVERGENCE AHEAD !!!
      // TODO: Can this be avoided?

      if (indmu != nbins_mu-1) {
         for(cRight = 1; indmu + cRight < rlimit; cRight++){
            if(GPUCELLMUSPACE(dev_fcount,cellIdx,indv,indmu + cRight) != 0){
               break;
            }
         }
         if( (GPUCELLMUSPACE(dev_fcount,cellIdx,indv,indmu + cRight) == 0) && (indmu + cRight == rlimit) ) {
            cRight  = 0;
         }
      }
      if (indmu != 0) {
         for(cLeft = 1; indmu - cLeft > llimit; cLeft++){
            if(GPUCELLMUSPACE(dev_fcount,cellIdx,indv,indmu - cLeft) != 0){
               break;
            }
         }
         if( (GPUCELLMUSPACE(dev_fcount,cellIdx,indv,indmu - cLeft) == 0) && (indmu - cLeft == llimit) ) {
            cLeft  = 0;
         }
      }

      __syncthreads();

      const Real Vmu = dVbins * (float(indv)+0.5);
      Realf dfdmu = 0.0;
      Realf dfdmu2 = 0.0;
      // Compute spatial derivatives
      if( (cRight != 0) || (cLeft != 0)) {
         dfdmu = (GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu + cRight) - GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu - cLeft))
                  /((cRight + cLeft)*dmubins) ;
      }
      if( (cRight != 0) && (cLeft != 0)) {
         dfdmu2 = ( (GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu + cRight) - GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu))/(cRight*dmubins)
                  - (GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu) - GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu - cLeft))
                  /(cLeft*dmubins) ) / (0.5 * dmubins * (cRight + cLeft));
      }

      // Compute time derivative
      const Realf mu    = (indmu+0.5)*dmubins - 1.0;
      const Realf Dmumu = dev_nu0Values[cellIdx]/2.0 * ( abs(mu)/(1.0 + abs(mu)) + epsilon ) * (1.0 - mu*mu);
      const Realf dDmu  = dev_nu0Values[cellIdx]/2.0 * ( (mu/abs(mu)) * ((1.0 - mu*mu)/((1.0 + abs(mu))*(1.0 + abs(mu)))) - 2.0*mu*( abs(mu)/(1.0 + abs(mu)) + epsilon));
      // We divide dfdt_mu by the normalization factor 2pi*v^2 already here.
      const Realf dfdt_mu_val = ( dDmu * dfdmu + Dmumu * dfdmu2 ) / (2.0 * M_PI * Vmu*Vmu);
      GPUCELLMUSPACE(dev_dfdt_mu,cellIdx,indv,indmu) = dfdt_mu_val;

      // Only consider CFL for non-negative phase-space cells above the sparsity threshold
      const Realf CellValue = GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu) / (2.0 * M_PI * Vmu*Vmu);
      const Realf absdfdt = abs(dfdt_mu_val); // Already scaled
      
      // Save calculated Ddt value
      if (absdfdt > 0.0 && CellValue > dev_sparsity[cellIdx]) {
         localDdtValues[threadIndex] = CellValue * PADCFL * (1.0/absdfdt);
      }

      // Reduction in shared memory
      __syncthreads();
      for (int s = blockDim.x / 2; s > GPUTHREADS/2; s >>= 1) {
         if (threadIndex < s) {
            localDdtValues[threadIndex] = min(localDdtValues[threadIndex], localDdtValues[threadIndex + s]);
         }
         __syncthreads();
      }
   }

   // Warp reduction
   if (threadIndex < GPUTHREADS) {
      gpuWarpSync();

      Real val = localDdtValues[threadIndex];
      for (int offset = GPUTHREADS/2; offset > 0; offset /= 2) {
         val = min(val, gpuKernelShflDown(val, offset));
      }

      // Write the result from the first thread of each block
      if (threadIndex == 0) {
         dev_potentialDdtValues[blockIdx.x] = val;
         dev_cellIdxKeys[blockIdx.x] = cellIdx;
      }
   }
}

__global__ void reduceDdtValues_kernel(
   int *dev_cellIdxKeys,
   Real *dev_potentialDdtValues,
   Real *dev_Ddt,
   int blocksPerSpatialCell
   ){

   int startIdx = blockIdx.x*blocksPerSpatialCell;
   int threadIndex = threadIdx.x;
   size_t cellIdx = dev_cellIdxKeys[startIdx];
   
   // Initiate shared memory values
   extern __shared__ Real localDdtValues[];
   localDdtValues[threadIndex] = std::numeric_limits<Real>::max();

   for(int i = threadIndex; i < blocksPerSpatialCell; i+=blockDim.x){
      localDdtValues[threadIndex] = min(localDdtValues[threadIndex], dev_potentialDdtValues[startIdx+i]);
   }

   // Reduction in shared memory
   __syncthreads();
   for (int s = blockDim.x / 2; s > GPUTHREADS/2; s >>= 1) {
      if (threadIndex < s) {
         localDdtValues[threadIndex] = min(localDdtValues[threadIndex], localDdtValues[threadIndex + s]);
      }
      __syncthreads();
   }
   // Warp reduction
   if (threadIndex < GPUTHREADS) {
      gpuWarpSync();

      Real val = localDdtValues[threadIndex];
      for (int offset = GPUTHREADS/2; offset > 0; offset /= 2) {
         val = min(val, gpuKernelShflDown(val, offset));
      }

      // Write the result from the first thread of each block
      if (threadIndex == 0) {
         dev_Ddt[cellIdx] = val;
      }
   }
}

__global__ void __launch_bounds__(Hashinator::defaults::MAX_BLOCKSIZE/2) dividefByCount_kernel(
   size_t *dev_smallCellIdxArray,
   Realf *dev_fmu,
   int *dev_fcount,
   int nbins_v,
   int nbins_mu,
   int maxThreadIndex
   ){

   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if(idx >= maxThreadIndex){return;}

   int indv = idx%nbins_v;
   int indmu = (idx/nbins_v)%nbins_mu;
   int spatialBlockIndex = idx/(nbins_v*nbins_mu); // Corresponds to index spatial and velocity blocks
   size_t cellIdx = dev_smallCellIdxArray[spatialBlockIndex];

   if (GPUCELLMUSPACE(dev_fcount,cellIdx,indv,indmu) == 0 || GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu) <= 0.0) {
      GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu) = 0;
   } else {
      GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu) = GPUCELLMUSPACE(dev_fmu,cellIdx,indv,indmu) / GPUCELLMUSPACE(dev_fcount,cellIdx,indv,indmu);
   }
}

__global__ void __launch_bounds__(Hashinator::defaults::MAX_BLOCKSIZE/2) getCellIndexArray_kernel(
   size_t *dev_cellIdxArray,
   size_t *dev_velocityIdxArray,
   size_t *dev_cellIdxStartCutoff,
   int numberOfComputedVelocityBlocks,
   int maxCellIndex
   ){
   
   size_t totalBlockIndex = blockIdx.x*blockDim.x + threadIdx.x;
   
   if(totalBlockIndex >= (size_t)numberOfComputedVelocityBlocks){return;}
   
   // Binary search
   int left = 0;
   int right = maxCellIndex - 1;
   int cellIndex = 0;

#ifdef DEBUG_SOLVERS
   assert(right>=left);
#endif

   while (left <= right) {
      int mid = (left + right) >> 1;
      if (dev_cellIdxStartCutoff[mid] <= (size_t)totalBlockIndex) {
         cellIndex = mid;
         left = mid + 1;
      } else {
         right = mid - 1;
      }
   }

   dev_cellIdxArray[totalBlockIndex] = cellIndex;
   dev_velocityIdxArray[totalBlockIndex] = totalBlockIndex - dev_cellIdxStartCutoff[cellIndex];
}

__global__ void __launch_bounds__(WID3) calculateDensity_kernel(
   Realf *dev_density,
   const vmesh::VelocityBlockContainer* __restrict__ const *dev_velocityBlockContainer
   ){

   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const int threadIndex = k*WID2+j*WID+i;
   const int cellIdx = blockIdx.x;
   const int blockSize = blockDim.x*blockDim.y*blockDim.z;

   __shared__ Realf localDensity[WID3];
   localDensity[threadIndex] = 0.0;

   const uint numberOfVelocityCells = dev_velocityBlockContainer[cellIdx]->size();
   
   for(uint velocityIdx = 0; velocityIdx < numberOfVelocityCells; velocityIdx++){
      localDensity[threadIndex] += dev_velocityBlockContainer[cellIdx]->getData()[velocityIdx*WID3+k*WID2+j*WID+i];
   }

   // Reduction in shared memory
   __syncthreads();
   for (int s = blockSize / 2; s > 0; s >>= 1) {
      if (threadIndex < s) {
         localDensity[threadIndex] += localDensity[threadIndex + s];
      }
      __syncthreads();
   }

   // Write the result from the first thread of each block
   if (threadIndex == 0) {
      dev_density[cellIdx] = localDensity[threadIndex];
   }
}

__global__ void __launch_bounds__(WID3) conserveMass_kernel(
   Realf *dev_densityPreAdjust,
   Realf *dev_densityPostAdjust,
   vmesh::VelocityBlockContainer* __restrict__ *dev_velocityBlockContainer
   ){
   
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const int cellIdx = blockIdx.x;

   if (dev_densityPostAdjust[cellIdx] == 0.0 || dev_densityPreAdjust[cellIdx] == dev_densityPostAdjust[cellIdx]){ return; }

   Realf adjustRatio = dev_densityPreAdjust[cellIdx]/dev_densityPostAdjust[cellIdx];

   const uint numberOfVelocityCells = dev_velocityBlockContainer[cellIdx]->size();
   
   for(uint velocityIdx = 0; velocityIdx < numberOfVelocityCells; velocityIdx++){
      dev_velocityBlockContainer[cellIdx]->getData()[velocityIdx*WID3+k*WID2+j*WID+i] *= adjustRatio;
   }
}

void pitchAngleDiffusion(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const uint popID){

   // Ensure nu0 dat file is read, if requested
   if (P::PADcoefficient < 0) {
      readNuArrayFromFile();
   }

   int nbins_v  = Parameters::PADvbins;
   int nbins_mu = Parameters::PADmubins;
   const Real dmubins = 2.0/nbins_mu;

   // resonance gap filling coefficient, not needed assuming even number of bins in mu-space
   const Real epsilon = 0.0;

   phiprof::Timer diffusionTimer {"pitch-angle-diffusion"};

   const std::vector<CellID>& LocalCells=getLocalCells();

   size_t numberOfLocalCells = LocalCells.size();

   std::vector<Real> dtTotalDiff(numberOfLocalCells, 0.0); // Diffusion time elapsed for each spatial cells

   int maxThreadsPerBlock = Hashinator::defaults::MAX_BLOCKSIZE;
   int blocksPerSpatialCell = (nbins_v*nbins_mu+maxThreadsPerBlock-1)/maxThreadsPerBlock;

   // Compute total number of velocity blocks
   int totalNumberOfVelocityBlocks = 0;
   for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells

      const auto CellID                  = LocalCells[CellIdx];
      spatial_cell::SpatialCell& cell                  = *mpiGrid[CellID];

      vmesh::LocalID numberOfVelocityBlocks = cell.get_number_of_velocity_blocks(popID);

      totalNumberOfVelocityBlocks += numberOfVelocityBlocks;
   } // End spatial cell loop

   // Allocate or reallocate memory if necessary
   gpu_batch_allocate(numberOfLocalCells, 0);
   gpu_pitch_angle_diffusion_allocate(numberOfLocalCells, nbins_v, nbins_mu, blocksPerSpatialCell, totalNumberOfVelocityBlocks);

   std::vector<bool> spatialLoopComplete(numberOfLocalCells, false);

   // Compute parameters
   #pragma omp parallel for
   for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells
      bool currentSpatialLoopComplete;
      Realf sparsity;
      std::array<Real,3> b;
      Real nu0;

      const auto CellID                  = LocalCells[CellIdx];
      SpatialCell& cell                  = *mpiGrid[CellID];

      computePitchAngleDiffusionParameters(
         cell,
         popID, CellIdx, currentSpatialLoopComplete,
         sparsity, b, nu0
      );

      // Save computed values to host
      spatialLoopComplete[CellIdx] = currentSpatialLoopComplete;
      host_sparsity[CellIdx] = sparsity;
      host_bValues[3*CellIdx] = b[0];
      host_bValues[3*CellIdx+1] = b[1];
      host_bValues[3*CellIdx+2] = b[2];
      host_nu0Values[CellIdx] = nu0;
   } // End spatial cell loop

   bool allSpatialCellTimeLoopsComplete = true;
   // Check if at least one cell needs to be calculated
   for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells
      if(!spatialLoopComplete[CellIdx]){
         allSpatialCellTimeLoopsComplete = false;
         break;
      }
   }

   if(allSpatialCellTimeLoopsComplete){
      // No need for any calculations
      return;
   }

   // Load CPU data
   const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
   const vmesh::MeshParameters& vMeshParams = vmesh::getMeshWrapper()->velocityMeshes->at(meshID);

   const Real Vmax   = 2*sqrt(3)*vMeshParams.meshLimits[1];
   Real dVbins = Vmax/nbins_v;

   #pragma omp parallel for
   for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells
      const auto CellID                  = LocalCells[CellIdx];
      spatial_cell::SpatialCell *cell                  = mpiGrid[CellID];

      host_bulkVX[CellIdx] = cell->parameters[CellParams::VX];
      host_bulkVY[CellIdx] = cell->parameters[CellParams::VY];
      host_bulkVZ[CellIdx] = cell->parameters[CellParams::VZ];

      host_VBCs[CellIdx] = cell->dev_get_velocity_blocks(popID);
   } // End spatial cell loop

   // Copy data to device
   CHK_ERR( gpuMemcpy(dev_bValues, host_bValues, 3*numberOfLocalCells*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_nu0Values, host_nu0Values, numberOfLocalCells*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_sparsity, host_sparsity, numberOfLocalCells*sizeof(Realf), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_bulkVX, host_bulkVX, numberOfLocalCells*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_bulkVY, host_bulkVY, numberOfLocalCells*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_bulkVZ, host_bulkVZ, numberOfLocalCells*sizeof(Real), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_VBCs, host_VBCs, numberOfLocalCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );

   if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
      dim3 threadsPerBlock_massConservation(WID, WID, WID);
      int blocksPerGrid_massConservation = numberOfLocalCells;

      // Ensure mass conservation
      calculateDensity_kernel<<<blocksPerGrid_massConservation, threadsPerBlock_massConservation>>>(
         dev_densityPreAdjust,
         dev_VBCs
      );
      
      CHK_ERR( gpuPeekAtLastError() );
   }

   while (!allSpatialCellTimeLoopsComplete) { // Substep loop

      // Compute maximum indices and construct cellIdx arrays
      int remappedCellIdx = 0;
      int numberOfComputedVelocityBlocks = 0;
      for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells
         // Add elements to remapped cellIdx array
         host_remappedCellIdxArray[CellIdx] = remappedCellIdx;

         if(spatialLoopComplete[CellIdx]){
            continue;
         }

         const auto CellID                  = LocalCells[CellIdx];
         spatial_cell::SpatialCell& cell                  = *mpiGrid[CellID];

         vmesh::LocalID numberOfVelocityBlocks = cell.get_number_of_velocity_blocks(popID);

         // Add elements to cellIdx arrays
         host_cellIdxStartCutoff[remappedCellIdx] = numberOfComputedVelocityBlocks;
         host_smallCellIdxArray[remappedCellIdx] = CellIdx;

         numberOfComputedVelocityBlocks += numberOfVelocityBlocks;
         remappedCellIdx++;
      } // End spatial cell loop
      
      int maxCellIndex = remappedCellIdx;
      int maxBlockIndex = numberOfComputedVelocityBlocks;

      // Copy data to device
      CHK_ERR( gpuMemcpy(dev_cellIdxStartCutoff, host_cellIdxStartCutoff, maxCellIndex*sizeof(size_t), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(dev_smallCellIdxArray, host_smallCellIdxArray, maxCellIndex*sizeof(size_t), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(dev_remappedCellIdxArray, host_remappedCellIdxArray, numberOfLocalCells*sizeof(size_t), gpuMemcpyHostToDevice) );

      // Initialize with zero values
      CHK_ERR( gpuMemset(dev_fmu, 0.0, numberOfLocalCells*nbins_v*nbins_mu*sizeof(Realf)) );
      CHK_ERR( gpuMemset(dev_fcount, 0, numberOfLocalCells*nbins_v*nbins_mu*sizeof(int)) );

      int totalThreadsPerBlock_getCellIndexArray = Hashinator::defaults::MAX_BLOCKSIZE/2; //Using Hashinator::defaults::MAX_BLOCKSIZE/2 = 512 blocks can lead to better streaming multiprocessor occupancy
      int maxThreadIndex_getCellIndexArray = numberOfComputedVelocityBlocks;
      int blocksPerGrid_getCellIndexArray = (maxThreadIndex_getCellIndexArray+totalThreadsPerBlock_getCellIndexArray-1)/totalThreadsPerBlock_getCellIndexArray;

      // Find spatial and velocity cell indices corresponding to each GPU block based on cutoffs,
      // so that each block will know the correct indeces in later kernels
      phiprof::Timer cellIdxArrayTimer {"getCellIndexArray_kernel"};
      getCellIndexArray_kernel<<<blocksPerGrid_getCellIndexArray, totalThreadsPerBlock_getCellIndexArray>>>(
         dev_cellIdxArray,
         dev_velocityIdxArray,
         dev_cellIdxStartCutoff,
         numberOfComputedVelocityBlocks,
         maxCellIndex
      );
      
      CHK_ERR( gpuPeekAtLastError() );
      cellIdxArrayTimer.stop();

      dim3 threadsPerBlock_build2dArrayOfFvmu(WID, WID, WID);
      int blocksPerGrid_build2dArrayOfFvmu = maxBlockIndex;

      // Build Fvmu array by dividing data to bins
      phiprof::Timer builFvmuTimer {"build2dArrayOfFvmu_kernel"};
      build2dArrayOfFvmu_kernel<<<blocksPerGrid_build2dArrayOfFvmu, threadsPerBlock_build2dArrayOfFvmu>>>(
         dev_cellIdxArray,
         dev_velocityIdxArray,
         dev_VBCs,
         dev_bulkVX,
         dev_bulkVY,
         dev_bulkVZ,
         dev_bValues,
         dev_fmu,
         dev_fcount,
         dVbins,
         dmubins,
         nbins_v,
         nbins_mu
      );
      
      CHK_ERR( gpuPeekAtLastError() );
      builFvmuTimer.stop();

      int totalThreadsPerBlock_dividefByCount = Hashinator::defaults::MAX_BLOCKSIZE/2; //Using Hashinator::defaults::MAX_BLOCKSIZE/2 = 512 blocks can lead to better streaming multiprocessor occupancy
      int maxThreadIndex_dividefByCount = numberOfLocalCells*nbins_v*nbins_mu;
      int blocksPerGrid_dividefByCount = (maxThreadIndex_dividefByCount+totalThreadsPerBlock_dividefByCount-1)/totalThreadsPerBlock_dividefByCount;

      // Divide by count
      phiprof::Timer divideFByCountTimer {"dividefByCount_kernel"};
      dividefByCount_kernel<<<blocksPerGrid_dividefByCount, totalThreadsPerBlock_dividefByCount>>>(
         dev_smallCellIdxArray,
         dev_fmu,
         dev_fcount,
         nbins_v,
         nbins_mu,
         maxThreadIndex_dividefByCount
      );
      
      CHK_ERR( gpuPeekAtLastError() );
      divideFByCountTimer.stop();
      
      int lastBlockSize = nbins_v*nbins_mu-(blocksPerSpatialCell-1)*maxThreadsPerBlock;
      int totalThreadsPerBlock_computeDerivativesCFLDdt;
      if(blocksPerSpatialCell == 1){
         totalThreadsPerBlock_computeDerivativesCFLDdt = nextPowerOfTwo(nbins_v*nbins_mu);
      }else{
         totalThreadsPerBlock_computeDerivativesCFLDdt = maxThreadsPerBlock;
      }
      int blocksPerGrid_computeDerivativesCFLDdt = numberOfLocalCells*blocksPerSpatialCell;
      int sharedMemorySize = totalThreadsPerBlock_computeDerivativesCFLDdt * sizeof(Real);

      // Compute derivatives and Ddt
      phiprof::Timer computeDerivativesTimer {"computeDerivativesCFLDdt_kernel"};
      computeDerivativesCFLDdt_kernel<<<blocksPerGrid_computeDerivativesCFLDdt, totalThreadsPerBlock_computeDerivativesCFLDdt, sharedMemorySize>>>(
         dev_smallCellIdxArray,
         dev_fcount,
         dev_fmu,
         dev_nu0Values,
         dev_dfdt_mu,
         dev_cellIdxKeys,
         dev_potentialDdtValues,
         dev_sparsity,
         dVbins,
         dmubins,
         epsilon,
         Parameters::PADCFL,
         nbins_v,
         nbins_mu,
         blocksPerSpatialCell,
         lastBlockSize
      );
      
      CHK_ERR( gpuPeekAtLastError() );
      computeDerivativesTimer.stop();

      int totalThreadsPerBlock_reduceDdtValues;
      if(blocksPerSpatialCell < maxThreadsPerBlock){
         totalThreadsPerBlock_reduceDdtValues = max(nextPowerOfTwo(blocksPerSpatialCell), GPUTHREADS);
      }else{
         totalThreadsPerBlock_reduceDdtValues = maxThreadsPerBlock;
      }
      int blocksPerGrid_reduceDdtValues = numberOfLocalCells;
      int sharedMemorySize_reduceDdtValues = totalThreadsPerBlock_reduceDdtValues * sizeof(Real);

      // Find minimum values of the calculated potential ddt values for each spatial cell
      phiprof::Timer reduceDdtValuesTimer {"reduceDdtValues_kernel"};
      reduceDdtValues_kernel<<<blocksPerGrid_reduceDdtValues, totalThreadsPerBlock_reduceDdtValues, sharedMemorySize_reduceDdtValues>>>(
         dev_cellIdxKeys,
         dev_potentialDdtValues,
         dev_Ddt,
         blocksPerSpatialCell
      );
      
      CHK_ERR( gpuPeekAtLastError() );
      reduceDdtValuesTimer.stop();
      CHK_ERR( gpuDeviceSynchronize() );

      CHK_ERR( gpuMemcpy(host_Ddt, dev_Ddt, numberOfLocalCells * sizeof(Real), gpuMemcpyDeviceToHost) );

      // Compute Ddt
      
      remappedCellIdx = 0;
      for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells
         if(spatialLoopComplete[CellIdx]){
            continue;
         }
         const Real RemainT  = Parameters::dt - dtTotalDiff[CellIdx]; //Remaining time before reaching simulation time step
         if (host_Ddt[CellIdx] > RemainT) {
            host_Ddt[CellIdx] = RemainT;
         }
         dtTotalDiff[CellIdx] += host_Ddt[CellIdx];
         remappedCellIdx++;
      } // End spatial cell loop

      CHK_ERR( gpuMemcpy(dev_Ddt, host_Ddt, numberOfLocalCells*sizeof(Real), gpuMemcpyHostToDevice) );

      dim3 threadsPerBlock_computeNewCellValues(WID, WID, WID);
      int blocksPerGrid_computeNewCellValues = maxBlockIndex;

      // Get new cell values
      phiprof::Timer newCellValuesTimer {"computeNewCellValues_kernel"};
      computeNewCellValues_kernel<<<blocksPerGrid_computeNewCellValues, threadsPerBlock_computeNewCellValues>>>(
         dev_cellIdxArray,
         dev_remappedCellIdxArray,
         dev_velocityIdxArray,
         dev_VBCs,
         dev_bulkVX,
         dev_bulkVY,
         dev_bulkVZ,
         dev_bValues,
         dev_dfdt_mu,
         dev_Ddt,
         dVbins,
         dmubins,
         nbins_v,
         nbins_mu
      );
      
      CHK_ERR( gpuPeekAtLastError() );
      newCellValuesTimer.stop();
      CHK_ERR( gpuDeviceSynchronize() );

      // Check if all cells are done
      for (size_t CellIdx = 0; CellIdx < numberOfLocalCells; CellIdx++) { // Iterate over all spatial cells
         allSpatialCellTimeLoopsComplete = true;
         if(dtTotalDiff[CellIdx] < Parameters::dt){
            allSpatialCellTimeLoopsComplete = false;
         }else{
            spatialLoopComplete[CellIdx] = true;
         }
      }

   } // End Time loop

   if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
      dim3 threadsPerBlock_massConservation(WID, WID, WID);
      int blocksPerGrid_massConservation = numberOfLocalCells;

      // Ensure mass conservation
      calculateDensity_kernel<<<blocksPerGrid_massConservation, threadsPerBlock_massConservation>>>(
         dev_densityPostAdjust,
         dev_VBCs
      );
         
      CHK_ERR( gpuPeekAtLastError() );

      conserveMass_kernel<<<blocksPerGrid_massConservation, threadsPerBlock_massConservation>>>(
         dev_densityPreAdjust,
         dev_densityPostAdjust,
         dev_VBCs
      );
         
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize() );
   }

   diffusionTimer.stop();
} // End function
