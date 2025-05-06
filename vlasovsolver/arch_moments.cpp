/*
 * This file is part of Vlasiator.
 * Copyright 2024-2025 University of Helsinki, CSC
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
#include "arch_moments.h"
#include "vlasovmover.h"
#include "../object_wrapper.h"
#include "../fieldsolver/fs_common.h" // divideIfNonZero()

#ifdef USE_GPU
#include "gpu_moments.h"
#endif

using namespace std;

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include contributions from
 * all existing particle populations. This function is AMR safe.
 * @param cell Spatial cell.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param doNotSkip If false, DO_NOT_COMPUTE cells are skipped.*/
void calculateCellMoments(spatial_cell::SpatialCell* cell,
                          const bool& computeSecond,
                          const bool& computePopulationMomentsOnly,
                          const bool& doNotSkip) {

   // Called once per cell. If doNotSkip == true, then DO_NOT_COMPUTE cells aren't skipped.
   if (!doNotSkip && cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
       return;
   }

   // Clear old moments to zero value
   if (computePopulationMomentsOnly == false) {
      cell->parameters[CellParams::RHOM  ] = 0.0;
      cell->parameters[CellParams::VX] = 0.0;
      cell->parameters[CellParams::VY] = 0.0;
      cell->parameters[CellParams::VZ] = 0.0;
      cell->parameters[CellParams::RHOQ  ] = 0.0;
      cell->parameters[CellParams::P_11] = 0.0;
      cell->parameters[CellParams::P_22] = 0.0;
      cell->parameters[CellParams::P_33] = 0.0;
      cell->parameters[CellParams::P_23] = 0.0;
      cell->parameters[CellParams::P_13] = 0.0;
      cell->parameters[CellParams::P_12] = 0.0;
   }

   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
      #endif
      const uint nBlocks = cell->get_velocity_mesh(popID)->size();
      Population &pop = cell->get_population(popID);
      if (nBlocks == 0) {
         pop.RHO = 0;
         for (int i=0; i<3; ++i) {
            pop.V[i]=0;
         }
         for (int i=0; i<nMom2; ++i) {
            pop.P[i]=0;
         }
         continue;
      }

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real charge = getObjectWrapper().particleSpecies[popID].charge;

      // Temporary array for storing moments
      Real array[nMom1] = {0};

      // Calculate species' contribution to first velocity moments
      phiprof::Timer firstMomentsTimer {"calcFirstMoments"};
      blockVelocityFirstMoments(blockContainer,
                                array,
                                nBlocks);
      firstMomentsTimer.stop();
      pop.RHO = array[0];
      pop.V[0] = divideIfNonZero(array[1], array[0]);
      pop.V[1] = divideIfNonZero(array[2], array[0]);
      pop.V[2] = divideIfNonZero(array[3], array[0]);

      if (!computePopulationMomentsOnly) {
         // Store species' contribution to bulk velocity moments
         cell->parameters[CellParams::RHOM  ] += array[0]*mass;
         cell->parameters[CellParams::VX] += array[1]*mass;
         cell->parameters[CellParams::VY] += array[2]*mass;
         cell->parameters[CellParams::VZ] += array[3]*mass;
         cell->parameters[CellParams::RHOQ  ] += array[0]*charge;
      }
   } // for-loop over particle species

   if(!computePopulationMomentsOnly) {
      cell->parameters[CellParams::VX] = divideIfNonZero(cell->parameters[CellParams::VX], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VY] = divideIfNonZero(cell->parameters[CellParams::VY], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VZ] = divideIfNonZero(cell->parameters[CellParams::VZ], cell->parameters[CellParams::RHOM]);
   }

   // Compute second moments only if requested
   if (computeSecond == false) {
      return;
   }

   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
      #endif
      const uint nBlocks = cell->get_velocity_mesh(popID)->size();
      if (nBlocks == 0) {
         continue;
      }

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;

      // Temporary array for storing moments
      Real array[nMom2] = {0};

      // Calculate species' contribution to second velocity moments
      phiprof::Timer secondMomentsTimer {"calcSecondMoments"};
      blockVelocitySecondMoments(blockContainer,
                                 cell->parameters[CellParams::VX],
                                 cell->parameters[CellParams::VY],
                                 cell->parameters[CellParams::VZ],
                                 array,
                                 nBlocks);
      secondMomentsTimer.stop();
      // Store species' contribution to bulk velocity moments
      Population &pop = cell->get_population(popID);
      for (size_t i=0; i<nMom2; ++i) {
         pop.P[i] = mass * array[i];
      }

      if (!computePopulationMomentsOnly) {
         cell->parameters[CellParams::P_11] += pop.P[0];
         cell->parameters[CellParams::P_22] += pop.P[1];
         cell->parameters[CellParams::P_33] += pop.P[2];
         cell->parameters[CellParams::P_23] += pop.P[3];
         cell->parameters[CellParams::P_13] += pop.P[4];
         cell->parameters[CellParams::P_12] += pop.P[5];
      }
   } // for-loop over particle species

}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _R variables. This function is AMR safe.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void calculateMoments_R(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond,
   const bool initialCompute) {

   // override with optimized GPU version to launch
   // single kernel accessing all cells at once (10x faster)
   #ifdef USE_GPU
   gpu_calculateMoments_R(mpiGrid,cells,computeSecond);
   return;
   #endif

   phiprof::Timer computeMomentsTimer {"Compute _R moments"};
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
#pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         phiprof::Timer computeMomentsCellTimer {"compute-moments-R-cell-first"};
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }
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

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         Population &pop = cell->get_population(popID);
         if (nBlocks == 0) {
            pop.RHO_R = 0;
            for (int i=0; i<3; ++i) {
               pop.V_R[i]=0;
            }
            for (int i=0; i<nMom2; ++i) {
               pop.P_R[i]=0;
            }
            continue;
         }
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array where the moments for this species are accumulated
         Real array[nMom1] = {0};

         // Calculate species' contribution to first velocity moments
         phiprof::Timer firstMomentsTimer {"calcFirstMoments_R"};
         blockVelocityFirstMoments(blockContainer,
                                   array,
                                   nBlocks);
         firstMomentsTimer.stop();
         // Store species' contribution to bulk velocity moments
         pop.RHO_R = array[0];
         pop.V_R[0] = divideIfNonZero(array[1], array[0]);
         pop.V_R[1] = divideIfNonZero(array[2], array[0]);
         pop.V_R[2] = divideIfNonZero(array[3], array[0]);

         cell->parameters[CellParams::RHOM_R  ] += array[0]*mass;
         cell->parameters[CellParams::VX_R] += array[1]*mass;
         cell->parameters[CellParams::VY_R] += array[2]*mass;
         cell->parameters[CellParams::VZ_R] += array[3]*mass;
         cell->parameters[CellParams::RHOQ_R  ] += array[0]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

#pragma omp parallel for schedule(static)
   for (size_t c=0; c<cells.size(); ++c) {
      phiprof::Timer computeMomentsCellTimer {"compute-moments-R-cell-bulkV"};
      SpatialCell* cell = mpiGrid[cells[c]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
         continue;
      }
      cell->parameters[CellParams::VX_R] = divideIfNonZero(cell->parameters[CellParams::VX_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VY_R] = divideIfNonZero(cell->parameters[CellParams::VY_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VZ_R] = divideIfNonZero(cell->parameters[CellParams::VZ_R], cell->parameters[CellParams::RHOM_R]);
   }

   // Compute second moments only if requested.
   if (computeSecond == false) {
      return;
   }

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
#pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         phiprof::Timer computeMomentsCellTimer {"compute-moments-R-cell-second"};
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         if (nBlocks == 0) {
            continue;
         }
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;

         // Temporary array where species' contribution to 2nd moments is accumulated
         Real array[nMom2] = {0};

         // Calculate species' contribution to second velocity moments
         phiprof::Timer secondMomentsTimer {"calcSecondMoments_R"};
         blockVelocitySecondMoments(blockContainer,
                                    cell->parameters[CellParams::VX_R],
                                    cell->parameters[CellParams::VY_R],
                                    cell->parameters[CellParams::VZ_R],
                                    array,
                                    nBlocks);
         secondMomentsTimer.stop();
         // Store species' contribution to 2nd bulk velocity moments
         Population &pop = cell->get_population(popID);
         for (size_t i = 0; i < nMom2; ++i) {
            pop.P_R[i] = mass * array[i];
         }

         cell->parameters[CellParams::P_11_R] += pop.P_R[0];
         cell->parameters[CellParams::P_22_R] += pop.P_R[1];
         cell->parameters[CellParams::P_33_R] += pop.P_R[2];
         cell->parameters[CellParams::P_23_R] += pop.P_R[3];
         cell->parameters[CellParams::P_13_R] += pop.P_R[4];
         cell->parameters[CellParams::P_12_R] += pop.P_R[5];
      } // for-loop over spatial cells
   } // for-loop over particle species
}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. Additionally, for each species, calculate the maximum
 * spatial time step so that CFL(spatial)=1. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _V variables. This function is AMR safe.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void calculateMoments_V(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond,
   const bool initialCompute) {

   // override with optimized GPU version to launch
   // single kernel accessing all cells at once (10x faster)
   #ifdef USE_GPU
   gpu_calculateMoments_V(mpiGrid,cells,computeSecond);
   return;
   #endif

   phiprof::Timer computeMomentsTimer {"Compute _V moments"};
   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
#pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         phiprof::Timer computeMomentsCellTimer {"compute-moments-V-cell"};
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }

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

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         Population &pop = cell->get_population(popID);
         if (nBlocks == 0) {
            pop.RHO_V = 0;
            for (int i=0; i<3; ++i) {
               pop.V_V[i]=0;
            }
            for (int i=0; i<nMom2; ++i) {
               pop.P_V[i]=0;
            }
            continue;
         }

         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array for storing moments
         Real array[nMom1] = {0};

         // Calculate species' contribution to first velocity moments
         phiprof::Timer firstMomentsTimer {"calcFirstMoments_V"};
         blockVelocityFirstMoments(blockContainer,
                                   array,
                                   nBlocks);
         firstMomentsTimer.stop();
         // Store species' contribution to bulk velocity moments
         pop.RHO_V = array[0];
         pop.V_V[0] = divideIfNonZero(array[1], array[0]);
         pop.V_V[1] = divideIfNonZero(array[2], array[0]);
         pop.V_V[2] = divideIfNonZero(array[3], array[0]);

         cell->parameters[CellParams::RHOM_V  ] += array[0]*mass;
         cell->parameters[CellParams::VX_V] += array[1]*mass;
         cell->parameters[CellParams::VY_V] += array[2]*mass;
         cell->parameters[CellParams::VZ_V] += array[3]*mass;
         cell->parameters[CellParams::RHOQ_V  ] += array[0]*charge;

      } // for-loop over spatial cells
   } // for-loop over particle species

#pragma omp parallel for schedule(static)
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* cell = mpiGrid[cells[c]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
         continue;
      }
      cell->parameters[CellParams::VX_V] = divideIfNonZero(cell->parameters[CellParams::VX_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VY_V] = divideIfNonZero(cell->parameters[CellParams::VY_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VZ_V] = divideIfNonZero(cell->parameters[CellParams::VZ_V], cell->parameters[CellParams::RHOM_V]);
   }

   // Compute second moments only if requested
   if (computeSecond == false) {
      return;
   }

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
#pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         phiprof::Timer computeMomentsCellTimer {"compute-moments-V-cell"};
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         if (nBlocks == 0) {
            continue;
         }

         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array where moments are stored
         Real array[nMom2] = {0};

         // Calculate species' contribution to second velocity moments
         phiprof::Timer secondMomentsTimer {"calcSecondMoments_V"};
         blockVelocitySecondMoments(blockContainer,
                                    cell->parameters[CellParams::VX_V],
                                    cell->parameters[CellParams::VY_V],
                                    cell->parameters[CellParams::VZ_V],
                                    array,
                                    nBlocks);
         secondMomentsTimer.stop();
         // Store species' contribution to 2nd bulk velocity moments
         Population &pop = cell->get_population(popID);
         for (size_t i = 0; i < nMom2; ++i) {
            pop.P_V[i] = mass * array[i];
         }

         cell->parameters[CellParams::P_11_V] += pop.P_V[0];
         cell->parameters[CellParams::P_22_V] += pop.P_V[1];
         cell->parameters[CellParams::P_33_V] += pop.P_V[2];
         cell->parameters[CellParams::P_23_V] += pop.P_V[3];
         cell->parameters[CellParams::P_13_V] += pop.P_V[4];
         cell->parameters[CellParams::P_12_V] += pop.P_V[5];
      } // for-loop over spatial cells
   } // for-loop over particle species
}
