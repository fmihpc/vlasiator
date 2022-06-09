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
#include "../fieldsolver/fs_common.h" // divideIfNonZero()
#include "../cuda_context.cuh"

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
         vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         const uint nBlocks = vmesh.size();
         if (nBlocks == 0) continue;
         
         masses[thread_id].at(popID) = getObjectWrapper().particleSpecies[popID].mass;
         charges[thread_id].at(popID) = getObjectWrapper().particleSpecies[popID].charge;
         meshDataPointers[thread_id].at(popID) = blockContainer.dev_getData();
         parameterPointers[thread_id].at(popID) = blockContainer.dev_getParameters();
         blockCounts[thread_id].at(popID) = nBlocks;
         // For now: Launch cuda transfers as data isn't yet fully resident
         phiprof::start("CUDA-HtoD");
         blockContainer.dev_Allocate(vmesh.size());
         blockContainer.dev_syncToDevice();
         // Also sync blockparams
         blockContainer.dev_syncParametersToDevice();
         phiprof::stop("CUDA-HtoD");
      }

      // Transfer metadata over
      phiprof::start("CUDA-HtoD");
      cudaMemcpyAsync(dev_meshDataPointers[thread_id], meshDataPointers[thread_id], nPopulations*sizeof(uint64_t), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemcpyAsync(dev_parameterPointers[thread_id], parameterPointers[thread_id], nPopulations*sizeof(uint64_t), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemcpyAsync(dev_blockCounts[thread_id], blockCounts[thread_id], nPopulations*sizeof(uint), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemcpyAsync(dev_masses[thread_id], masses[thread_id], nPopulations*sizeof(Real), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemcpyAsync(dev_charges[thread_id], charges[thread_id], nPopulations*sizeof(Real), cudaMemcpyHostToDevice, cudaStreamList[thread_id]);
      cudaMemsetAsync(dev_momentArrays[thread_id], 0, nMoments*(nPopulations+1)*sizeof(Real), cudaStreamList[thread_id]);

      phiprof::stop("CUDA-HtoD");

      // Now launch kernel for this spatial cell, all populations
      calculate_moments_glue(
         dev_meshDataPointers[thread_id],
         dev_parameterPointers[thread_id],
         dev_blockCounts[thread_id],
         dev_masses[thread_id],
         dev_charges[thread_id],
         dev_momentArrays[thread_id],
         nPopulations,
         computeSecond,
         cudaStreamList[thread_id]
         );
      
      // Transfer momentArrays back
      phiprof::start("CUDA-DtoH");
      cudaMemcpyAsync(momentArrays[thread_id], dev_momentArrays[thread_id], nMoments*(nPopulations+1)*sizeof(Real), cudaMemcpyDeviceToHost, cudaStreamList[thread_id]);
      phiprof::stop("CUDA-DtoH");

      // Store species' contributions to velocity moments
      for (uint popID=0; popID<nPopulations; ++popID) {         
         Population & pop = cell->get_population(popID);
         pop.RHO_V = momentArrays[thread_id].at(popID)[0];
         pop.V_V[0] = momentArrays[thread_id].at(popID)[1];
         pop.V_V[1] = momentArrays[thread_id].at(popID)[2];
         pop.V_V[2] = momentArrays[thread_id].at(popID)[3];
         // pop.RHOQ_V = momentArrays[thread_id].at(popID)[4];
         pop.P_V[0] = momentArrays[thread_id].at(popID)[5];
         pop.P_V[1] = momentArrays[thread_id].at(popID)[6];
         pop.P_V[2] = momentArrays[thread_id].at(popID)[7];         
      }         
      cell->parameters[CellParams::RHOM_V] = momentArrays[thread_id].at(nPopulations)[0];
      cell->parameters[CellParams::VX_V  ] = momentArrays[thread_id].at(nPopulations)[1];
      cell->parameters[CellParams::VY_V  ] = momentArrays[thread_id].at(nPopulations)[2];
      cell->parameters[CellParams::VZ_V  ] = momentArrays[thread_id].at(nPopulations)[3];
      cell->parameters[CellParams::RHOQ_V] = momentArrays[thread_id].at(nPopulations)[4];
      cell->parameters[CellParams::P_11_V] = momentArrays[thread_id].at(nPopulations)[5];
      cell->parameters[CellParams::P_22_V] = momentArrays[thread_id].at(nPopulations)[6];
      cell->parameters[CellParams::P_33_V] = momentArrays[thread_id].at(nPopulations)[7];

   } // for-loop over spatial cells

   phiprof::stop("CUDA Compute _V moments");
   return;
}
