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
#include "cpu_moments.h"
#include "../vlasovmover.h"
#include "../object_wrapper.h"
#include "../fieldsolver/fs_common.h" // divideIfNonZero()

using namespace std;

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include contributions from
 * all existing particle populations. This function is VAMR safe.
 * @param cell Spatial cell.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param computePopulationMomentsOnly Do not update the combined moments in CellParams if true
 * @param doNotSkip If false, DO_NOT_COMPUTE cells are skipped.*/
void calculateCellMoments(spatial_cell::SpatialCell* cell,
                          const bool& computeSecond,
                          const bool& computePopulationMomentsOnly,
                          const bool& doNotSkip) {
   // if doNotSkip == true then the first clause is false and we will never return,
   // i.e. always compute, otherwise we skip DO_NOT_COMPUTE cells
   bool skipMoments = false;
   if (!doNotSkip && cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
       skipMoments = true;
   }

   // Clear old moments to zero value
   if (skipMoments == false && computePopulationMomentsOnly == false) {
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
   if (skipMoments == false) {
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         Population & pop = cell->get_population(popID);
         if (blockContainer.size() == 0) {
            pop.RHO = 0;
            for (int i=0; i<3; ++i) {
               pop.V[i]=0;
               pop.P[i]=0;
            }
            continue;
         }
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array for storing moments
         Real array[4];
         for (int i=0; i<4; ++i) array[i] = 0.0;

         // Calculate species' contribution to first velocity moments
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocityFirstMoments(data+blockLID*WID3,
                                      blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                      array);
         }

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

      if (!computePopulationMomentsOnly) {
         cell->parameters[CellParams::VX] = divideIfNonZero(cell->parameters[CellParams::VX], cell->parameters[CellParams::RHOM]);
         cell->parameters[CellParams::VY] = divideIfNonZero(cell->parameters[CellParams::VY], cell->parameters[CellParams::RHOM]);
         cell->parameters[CellParams::VZ] = divideIfNonZero(cell->parameters[CellParams::VZ], cell->parameters[CellParams::RHOM]);
      }
   }

    // Compute second moments only if requested
   if (computeSecond == false) {
      return;
   }

   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
      if (blockContainer.size() == 0) {
         continue;
      }
      const Realf* data       = blockContainer.getData();
      const Real* blockParams = blockContainer.getParameters();
      const Real mass = getObjectWrapper().particleSpecies[popID].mass;

      // Temporary array for storing moments
      std::vector<Real> array(6, 0.0);

      // Calculate species' contribution to second velocity moments
      Population & pop = cell->get_population(popID);
      for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
         blockVelocitySecondMoments(data+blockLID*WID3,
                                    blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                    cell->parameters[CellParams::VX],
                                    cell->parameters[CellParams::VY],
                                    cell->parameters[CellParams::VZ],
                                    array);
      }

      // Store species' contribution to bulk velocity moments
      for (size_t i = 0; i < array.size(); ++i) {
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
 * are stored to SpatialCell::parameters in _R variables. This function is VAMR safe.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void calculateMoments_R(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const std::vector<CellID>& cells,
        const bool& computeSecond) {

    phiprof::Timer momentsTimer {"compute-moments-n"};

    for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
       #pragma omp parallel for
       for (size_t c=0; c<cells.size(); ++c) {
          SpatialCell* cell = mpiGrid[cells[c]];

          if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
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

          vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
          Population & pop = cell->get_population(popID);
          if (blockContainer.size() == 0) {
             pop.RHO_R = 0;
            for (int i = 0; i < 3; ++i) {
               pop.V_R[i] = 0;
            }
            for (int i = 0; i < 6; ++i) {
               pop.P_R[i] = 0;
            }
             continue;
          }
          const Realf* data       = blockContainer.getData();
          const Real* blockParams = blockContainer.getParameters();
          const Real mass = getObjectWrapper().particleSpecies[popID].mass;
          const Real charge = getObjectWrapper().particleSpecies[popID].charge;

          #ifdef DEBUG_MOMENTS
          bool ok = true;
          if (data == NULL && blockContainer.size() > 0) ok = false;
          if (blockParams == NULL && blockContainer.size() > 0) ok = false;
          if (ok == false) {
             stringstream ss;
             ss << "ERROR in moment calculation in " << __FILE__ << ":" << __LINE__ << endl;
             ss << "\t &data = " << data << "\t &blockParams = " << blockParams << endl;
             ss << "\t size = " << blockContainer.size() << endl;
             cerr << ss.str();
             exit(1);
          }
          #endif

          // Temporary array where the moments for this species are accumulated
          Real array[4];
          for (int i=0; i<4; ++i) array[i] = 0.0;

          // Calculate species' contribution to first velocity moments
          for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocityFirstMoments(data+blockLID*WID3,
                                      blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                      array);
          } // for-loop over velocity blocks

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

    #pragma omp parallel for
    for (size_t c=0; c<cells.size(); ++c) {
       SpatialCell* cell = mpiGrid[cells[c]];
       if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
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
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         if (blockContainer.size() == 0) {
            continue;
         }
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;

         // Temporary array where species' contribution to 2nd moments is accumulated
         std::vector<Real> array(6, 0.0);

         // Calculate species' contribution to second velocity moments
         Population & pop = cell->get_population(popID);
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocitySecondMoments(data+blockLID*WID3,
                                       blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                       cell->parameters[CellParams::VX_R],
                                       cell->parameters[CellParams::VY_R],
                                       cell->parameters[CellParams::VZ_R],
                                       array);
         } // for-loop over velocity blocks

         // Store species' contribution to 2nd bulk velocity moments
         for (size_t i = 0; i < array.size(); ++i) {
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
 * are stored to SpatialCell::parameters in _V variables. This function is VAMR safe.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void calculateMoments_V(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const std::vector<CellID>& cells,
        const bool& computeSecond) {

   phiprof::Timer momentsTimer {"Compute _V moments"};

   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
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

         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         Population & pop = cell->get_population(popID);
         if (blockContainer.size() == 0) {
            pop.RHO_V = 0;
            for (int i = 0; i < 3; ++i) {
               pop.V_V [i] = 0;
            }
            for (int i = 0; i < 6; ++i) {
               pop.P_V [i] = 0;
            }
            continue;
         }
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array for storing moments
         Real array[4];
         for (int i=0; i<4; ++i) array[i] = 0.0;

         // Calculate species' contribution to first velocity moments
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocityFirstMoments(data+blockLID*WID3,
                                      blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                      array);
         }

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

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* cell = mpiGrid[cells[c]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
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
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }

         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         if (blockContainer.size() == 0) {
            continue;
         }
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;

         // Temporary array where moments are stored
         std::vector<Real> array(6, 0.0);

         // Calculate species' contribution to second velocity moments
         Population & pop = cell->get_population(popID);
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocitySecondMoments(
                                       data+blockLID*WID3,
                                       blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                       cell->parameters[CellParams::VX_V],
                                       cell->parameters[CellParams::VY_V],
                                       cell->parameters[CellParams::VZ_V],
                                       array);
         } // for-loop over velocity blocks

         // Store species' contribution to 2nd bulk velocity moments
         for (size_t i = 0; i < array.size(); ++i) {
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
