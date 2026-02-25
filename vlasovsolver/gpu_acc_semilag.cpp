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

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <phiprof.hpp>
#include "../definitions.h"

#include "gpu_acc_semilag.hpp"
#include "cpu_acc_intersections.hpp"
#include "gpu_acc_map.hpp"
#include "../spatial_cells/block_adjust_gpu.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
  \brief Propagates the distribution function in velocity space of given list of
  real space cells using a semi-Lagrangian acceleration approach. GPU version.

  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

  @param mpiGrid DCCRG container of spatial cells
  @param acceleratedCells vector of cells for which to perform acceleration
  @param popID ID of the accelerated particle species.
  @param map_order Order [0,2] in which vx,vy,vz mappings are performed.
*/

void gpu_accelerate_cells(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& acceleratedCells,
                          const uint popID,
                          const uint map_order
   ) {

   CHK_ERR( gpuDeviceSynchronize() );
   phiprof::Timer verificationTimer {"gpu ACC allocation verifications"};
   const uint nCells = (uint)acceleratedCells.size();
   gpu_batch_allocate(nCells,0);
   verificationTimer.stop();

   // Calculate intersections (should be constant cost per cell). Also reduces
   // the largest found block count in order to ensure allocations.
   int intersections_id {phiprof::initializeTimer("cell-compute-intersections")};
   uint gpuMaxBlockCount = 0;
   #pragma omp parallel
   {
      uint threadGpuMaxBlockCount = 0;
      #pragma omp for schedule(static)
      for (size_t cellIndex=0; cellIndex<acceleratedCells.size(); ++cellIndex) {
         const CellID cid = acceleratedCells[cellIndex];
         SpatialCell* SC = mpiGrid[cid];
         Population& pop = SC->get_population(popID);
         compute_cell_intersections(SC, popID, map_order, pop.subcycleDt, intersections_id);

         const vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
         const uint blockCount = vmesh->size();
         threadGpuMaxBlockCount = std::max(threadGpuMaxBlockCount,blockCount);
      }
      #pragma omp critical
      {
         gpuMaxBlockCount = std::max(gpuMaxBlockCount,threadGpuMaxBlockCount);
      }
   }

   // Do some overall preparation regarding dimensions and acceleration order
   const uint D0 = (*vmesh::getMeshWrapper()->velocityMeshes)[popID].gridLength[0];
   const uint D1 = (*vmesh::getMeshWrapper()->velocityMeshes)[popID].gridLength[1];
   const uint D2 = (*vmesh::getMeshWrapper()->velocityMeshes)[popID].gridLength[2];

   std::vector<int> dimOrder(3);
   switch(map_order) {
      case 0: { //Map order XYZ
         dimOrder={0,1,2};
         break;
      }
      case 1: { //Map order YZX
         dimOrder={1,2,0};
         break;
      }
      case 2: { //Map order ZXY
         dimOrder={2,0,1};
         break;
      }
      default:
         std::cerr<<"ERROR! Incorrect map_order "<<map_order<<"!"<<std::endl;
         abort();
   }

   /**
      Loop over three velocity dimensions, based on map_order,
      and accelerate all cells for that dimension.
   */
   for (int dimIndex = 0; dimIndex<3; ++dimIndex) {
      int dimension = dimOrder[dimIndex];

      // Gather up-to-date pointers for cell contents
      uint gpuMaxBlockCount = 0;
      #pragma omp parallel
      {
         uint threadGpuMaxBlockCount = 0;
         #pragma omp for schedule(static)
         for (size_t cellIndex=0; cellIndex<acceleratedCells.size(); ++cellIndex) {
            const CellID cid = acceleratedCells[cellIndex];
            SpatialCell* SC = mpiGrid[cid];
            const uint blockCount = SC->get_velocity_mesh(popID)->size();
            // Ensure per-cell allocations
            SC->setReservation(popID,blockCount);
            SC->applyReservation(popID);

            threadGpuMaxBlockCount = std::max(threadGpuMaxBlockCount,blockCount);
            // Store pointers in batch buffers
            (GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes))[cellIndex] = SC->dev_get_velocity_mesh(popID);
            (GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBCs))[cellIndex] = SC->dev_get_velocity_blocks(popID);
            (GET_POINTER(gpuMemoryManager, Real, host_minValues))[cellIndex] = SC->getVelocityBlockMinValue(popID);
            (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec))[cellIndex] = SC->dev_velocity_block_with_content_list;
            (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new))[cellIndex] = SC->dev_list_with_replace_new;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_delete))[cellIndex] = SC->dev_list_delete;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_to_replace))[cellIndex] = SC->dev_list_to_replace;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_with_replace_old))[cellIndex] = SC->dev_list_with_replace_old;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*cellIndex] = SC->dev_velocity_block_with_content_map;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*cellIndex+1] = SC->dev_velocity_block_with_no_content_map;
         }
         #pragma omp critical
         {
            gpuMaxBlockCount = std::max(gpuMaxBlockCount,threadGpuMaxBlockCount);
         }
      }
      // Ensure accelerator has enough temporary memory allocated
      verificationTimer.start();
      gpu_vlasov_allocate(gpuMaxBlockCount);
      gpu_acc_allocate(gpuMaxBlockCount);
      verificationTimer.stop();

      // Copy pointers and counters over to device
      phiprof::Timer copyTimer {"copy pointer addresses to device"};
      CHK_ERR( gpuMemset(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBefore), 0, nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMemset(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nAfter), 0, nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMemset(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBlocksToChange), 0, nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMemset(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess), 0, nCells*sizeof(vmesh::LocalID)) );

      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps), 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes), nCells*sizeof(vmesh::VelocityMesh*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, Real, dev_minValues), GET_POINTER(gpuMemoryManager, Real, host_minValues), nCells*sizeof(Real), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec), nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new), nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete), GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_delete), nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace), GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_to_replace), nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old), GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_with_replace_old), nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs), GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBCs), nCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );
      copyTimer.stop();

      string profName = "accelerate "+getObjectWrapper().particleSpecies[popID].name;
      phiprof::Timer accTimer {profName};

      // used when computing id of target block.
      uint block_indices_to_id[3] = {0, 0, 0};
      uint block_indices_to_probe[3] = {0, 0, 0};
      uint cell_indices_to_id[3] = {0, 0, 0};

      // Find probe cube extents as well
      int Dacc=0, Dother=0;

      switch (dimension) {
         case 0: /* i and k coordinates have been swapped */
            /* set values in array that is used to convert block indices to id using a dot product */
            block_indices_to_id[0] = D0*D1;
            block_indices_to_id[1] = D0;
            block_indices_to_id[2] = 1;

            /* set values in array that is used to convert block indices to position in probe cube
               propagate along X, flatten Y+Z */
            block_indices_to_probe[0] = D1*D2;
            block_indices_to_probe[1] = D2;
            block_indices_to_probe[2] = 1;
            Dacc = D0;
            Dother = D1*D2;

            /* set values in array that is used to convert block indices to id using a dot product */
            cell_indices_to_id[0] = WID2;
            cell_indices_to_id[1] = WID;
            cell_indices_to_id[2] = 1;
            break;
         case 1: /* j and k coordinates have been swapped */
            /* set values in array that is used to convert block indices to id using a dot product */
            block_indices_to_id[0] = 1;
            block_indices_to_id[1] = D0*D1;
            block_indices_to_id[2] = D0;

            /* set values in array that is used to convert block indices to position in probe cube
               propagate along Y, flatten X+Z */
            block_indices_to_probe[0] = D2;
            block_indices_to_probe[1] = D0*D2;
            block_indices_to_probe[2] = 1;
            Dacc = D1;
            Dother = D0*D2;

            /* set values in array that is used to convert block indices to id using a dot product */
            cell_indices_to_id[0] = 1;
            cell_indices_to_id[1] = WID2;
            cell_indices_to_id[2] = WID;
            break;
         case 2:
            /* set values in array that is used to convert block indices to id using a dot product */
            block_indices_to_id[0] = 1;
            block_indices_to_id[1] = D0;
            block_indices_to_id[2] = D0*D1;

            /* set values in array that is used to convert block indices to position in probe cube
               propagate along Z, flatten X+Y */
            block_indices_to_probe[0] = D1;
            block_indices_to_probe[1] = 1;
            block_indices_to_probe[2] = D0*D1;
            Dacc = D2;
            Dother = D0*D1;

            /* set values in array that is used to convert block indices to id using a dot product. */
            cell_indices_to_id[0] = 1;
            cell_indices_to_id[1] = WID;
            cell_indices_to_id[2] = WID2;
            break;
         default:
            std::cerr<<"Invalid dimension "<<dimension<<"!"<<std::endl;
            abort();
      }

      // Copy indexing information to device. To be tested: might be faster to pass a single
      // device-side struct or just 9 plain arguments?
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, uint, gpu_cell_indices_to_id), cell_indices_to_id, 3*sizeof(uint), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, uint, gpu_block_indices_to_id), block_indices_to_id, 3*sizeof(uint), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, uint, gpu_block_indices_to_probe), block_indices_to_probe, 3*sizeof(uint), gpuMemcpyHostToDevice) );

      // Select correct intersections for each mapping order
      #pragma omp parallel for
      for (size_t cellIndex=0; cellIndex<acceleratedCells.size(); ++cellIndex) {
         const CellID cellID = acceleratedCells[cellIndex];
         Population& pop = mpiGrid[cellID]->get_population(popID);
         // Place intersections into array so that propagation direction is k-coordinate ("z")
         switch (dimension) {
            case 0:
               // X: swap intersection i and k coordinates
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+0]=(Realf)pop.intersection_x;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+1]=(Realf)pop.intersection_x_dk;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+2]=(Realf)pop.intersection_x_dj;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+3]=(Realf)pop.intersection_x_di;
               break;
            case 1:
               // Y: swap intersection j and k coordinates
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+0]=(Realf)pop.intersection_y;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+1]=(Realf)pop.intersection_y_di;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+2]=(Realf)pop.intersection_y_dk;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+3]=(Realf)pop.intersection_y_dj;
               break;
            case 2:
               // Z: k remains propagation coordinate, no swaps
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+0]=(Realf)pop.intersection_z;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+1]=(Realf)pop.intersection_z_di;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+2]=(Realf)pop.intersection_z_dj;
               (GET_POINTER(gpuMemoryManager, Realf, host_intersections))[cellIndex*4+3]=(Realf)pop.intersection_z_dk;
               break;
            default:
               std::cerr<<"Invalid dimension "<<dimension<<"!"<<std::endl;
               abort();
         }
      }
      // Send intersection data to device
      CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, Realf, dev_intersections), GET_POINTER(gpuMemoryManager, Realf, host_intersections), 4*nCells*sizeof(Realf), gpuMemcpyHostToDevice) );

      // Call acceleration solver in chunks, the size of which is determined by the GPU
      // Vlasov allocation number.
      const uint maxChunkSize = gpu_getAllocationCount();

      uint queuedCells = 0;
      uint checkedCells = 0;
      size_t cumulativeOffset = 0;
      uint chunk = 0;
      std::vector<CellID> launchCells;

      for (size_t cellIndex=0; cellIndex<nCells; ++cellIndex) {
         CellID cid = acceleratedCells[cellIndex];
         SpatialCell* SC = mpiGrid[cid];
         const uint blockCount = SC->get_velocity_mesh(popID)->size();

         if (blockCount > 0) {
            // Only accelerate non-empty cells
            launchCells.push_back(cid);
            queuedCells++;
         }
         // Keep track of all checked cells for cumulative offset into pointer buffers
         checkedCells++;

         // Once enough cells have been gathered into the chunk, or we have evaluated
         // the last of potential cells, Launch the acceleration solver for this chunk.
         // Will not launch if no cells to be accelerated are left.
         if (queuedCells == maxChunkSize || ( (cellIndex==nCells-1) && (queuedCells > 0) )) {
            // Phiprof timer
            string timerName = "semilag-acc-dim"+std::to_string(dimension)+"-chunk";
            //timerName += "-"+std::to_string(chunk); // Optional: phiprof label for chunk id
            phiprof::Timer accChunkTimer {timerName};

            gpu_acc_map_1d(mpiGrid,
                           launchCells,
                           popID,
                           dimension,
                           Dacc,
                           Dother,
                           cumulativeOffset
               );
            cumulativeOffset += checkedCells;
            queuedCells = 0;
            checkedCells = 0;
            chunk++;
            launchCells.clear();
         }
      }
   }
}
