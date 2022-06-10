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
#include "cuda_moments_kernel.cuh"

using namespace std;

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
         vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         const uint nBlocks = vmesh.size();
         if (nBlocks == 0) continue;

         host_momentInfos[thread_id][popID].mass = getObjectWrapper().particleSpecies[popID].mass;
         host_momentInfos[thread_id][popID].charge = getObjectWrapper().particleSpecies[popID].charge;
         host_momentInfos[thread_id][popID].blockCount = nBlocks;
         host_momentInfos[thread_id][popID].parameterPointer = blockContainer.dev_getParameters();
         host_momentInfos[thread_id][popID].meshDataPointer = blockContainer.dev_getData();

         // For now: Launch cuda transfers as data isn't yet fully resident
         phiprof::start("CUDA-HtoD");
         blockContainer.dev_Allocate(vmesh.size());
         blockContainer.dev_syncBlocksToDevice();
         // Also sync blockparams
         blockContainer.dev_syncParametersToDevice();
         phiprof::stop("CUDA-HtoD");
      }

      // Transfer metadata over
      phiprof::start("CUDA-HtoD");
      cudaMemcpyAsync(dev_momentInfos[thread_id], host_momentInfos[thread_id], nPopulations*sizeof(MomentInfo), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemsetAsync(dev_momentArrays[thread_id], 0, nMoments*(nPopulations+1)*sizeof(Real), cudaStreamList[thread_id]);
      phiprof::stop("CUDA-HtoD");

      // Now launch kernel for this spatial cell, all populations
      calculate_moments_glue(
         dev_momentInfos[thread_id],
         dev_momentArrays[thread_id],
         nPopulations,
         computeSecond,
         cudaStreamList[thread_id]
         );

      // Transfer momentArrays back
      phiprof::start("CUDA-DtoH");
      cudaMemcpyAsync(host_momentArrays[thread_id], dev_momentArrays[thread_id], nMoments*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]);
      phiprof::stop("CUDA-DtoH");

      // Store species' contributions to velocity moments
      for (uint popID=0; popID<nPopulations; ++popID) {
         Population & pop = cell->get_population(popID);
         pop.RHO_V = host_momentArrays[thread_id][popID*nMoments + 0];
         pop.V_V[0] = host_momentArrays[thread_id][popID*nMoments + 1];
         pop.V_V[1] = host_momentArrays[thread_id][popID*nMoments + 2];
         pop.V_V[2] = host_momentArrays[thread_id][popID*nMoments + 3];
         // pop.RHOQ_V = host_momentArrays[thread_id][popID*nMoments + 4];
         pop.P_V[0] = host_momentArrays[thread_id][popID*nMoments + 5];
         pop.P_V[1] = host_momentArrays[thread_id][popID*nMoments + 6];
         pop.P_V[2] = host_momentArrays[thread_id][popID*nMoments + 7];
      }
      cell->parameters[CellParams::RHOM_V] = host_momentArrays[thread_id][nPopulations*nMoments + 0];
      cell->parameters[CellParams::VX_V  ] = host_momentArrays[thread_id][nPopulations*nMoments + 1];
      cell->parameters[CellParams::VY_V  ] = host_momentArrays[thread_id][nPopulations*nMoments + 2];
      cell->parameters[CellParams::VZ_V  ] = host_momentArrays[thread_id][nPopulations*nMoments + 3];
      cell->parameters[CellParams::RHOQ_V] = host_momentArrays[thread_id][nPopulations*nMoments + 4];
      cell->parameters[CellParams::P_11_V] = host_momentArrays[thread_id][nPopulations*nMoments + 5];
      cell->parameters[CellParams::P_22_V] = host_momentArrays[thread_id][nPopulations*nMoments + 6];
      cell->parameters[CellParams::P_33_V] = host_momentArrays[thread_id][nPopulations*nMoments + 7];

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
         vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         const uint nBlocks = vmesh.size();
         if (nBlocks == 0) continue;

         host_momentInfos[thread_id][popID].mass = getObjectWrapper().particleSpecies[popID].mass;
         host_momentInfos[thread_id][popID].charge = getObjectWrapper().particleSpecies[popID].charge;
         host_momentInfos[thread_id][popID].blockCount = nBlocks;
         host_momentInfos[thread_id][popID].parameterPointer = blockContainer.dev_getParameters();
         host_momentInfos[thread_id][popID].meshDataPointer = blockContainer.dev_getData();

         // For now: Launch cuda transfers as data isn't yet fully resident
         phiprof::start("CUDA-HtoD");
         blockContainer.dev_Allocate(vmesh.size());
         blockContainer.dev_syncBlocksToDevice();
         // Also sync blockparams
         blockContainer.dev_syncParametersToDevice();
         phiprof::stop("CUDA-HtoD");
      }

      // Transfer metadata over
      phiprof::start("CUDA-HtoD");
      cudaMemcpyAsync(dev_momentInfos[thread_id], host_momentInfos[thread_id], nPopulations*sizeof(MomentInfo), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemsetAsync(dev_momentArrays[thread_id], 0, nMoments*(nPopulations+1)*sizeof(Real), cudaStreamList[thread_id]);
      phiprof::stop("CUDA-HtoD");

      // Now launch kernel for this spatial cell, all populations
      calculate_moments_glue(
         dev_momentInfos[thread_id],
         dev_momentArrays[thread_id],
         nPopulations,
         computeSecond,
         cudaStreamList[thread_id]
         );

      // Transfer momentArrays back
      phiprof::start("CUDA-DtoH");
      cudaMemcpyAsync(host_momentArrays[thread_id], dev_momentArrays[thread_id], nMoments*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]);
      phiprof::stop("CUDA-DtoH");

      // Store species' contributions to velocity moments
      for (uint popID=0; popID<nPopulations; ++popID) {
         Population & pop = cell->get_population(popID);
         pop.RHO_R = host_momentArrays[thread_id][popID*nMoments + 0];
         pop.V_R[0] = host_momentArrays[thread_id][popID*nMoments + 1];
         pop.V_R[1] = host_momentArrays[thread_id][popID*nMoments + 2];
         pop.V_R[2] = host_momentArrays[thread_id][popID*nMoments + 3];
         // pop.RHOQ_R = host_momentArrays[thread_id][popID*nMoments + 4];
         pop.P_R[0] = host_momentArrays[thread_id][popID*nMoments + 5];
         pop.P_R[1] = host_momentArrays[thread_id][popID*nMoments + 6];
         pop.P_R[2] = host_momentArrays[thread_id][popID*nMoments + 7];
      }
      cell->parameters[CellParams::RHOM_R] = host_momentArrays[thread_id][nPopulations*nMoments + 0];
      cell->parameters[CellParams::VX_R  ] = host_momentArrays[thread_id][nPopulations*nMoments + 1];
      cell->parameters[CellParams::VY_R  ] = host_momentArrays[thread_id][nPopulations*nMoments + 2];
      cell->parameters[CellParams::VZ_R  ] = host_momentArrays[thread_id][nPopulations*nMoments + 3];
      cell->parameters[CellParams::RHOQ_R] = host_momentArrays[thread_id][nPopulations*nMoments + 4];
      cell->parameters[CellParams::P_11_R] = host_momentArrays[thread_id][nPopulations*nMoments + 5];
      cell->parameters[CellParams::P_22_R] = host_momentArrays[thread_id][nPopulations*nMoments + 6];
      cell->parameters[CellParams::P_33_R] = host_momentArrays[thread_id][nPopulations*nMoments + 7];

   } // for-loop over spatial cells

   phiprof::stop("CUDA Compute _R moments");
   return;
}
