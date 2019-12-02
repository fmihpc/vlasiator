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
 * all existing particle populations. This function is AMR safe.
 * @param cell Spatial cell.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param doNotSkip If false, DO_NOT_COMPUTE cells are skipped.*/
void calculateCellMoments(spatial_cell::SpatialCell* cell,
                          const bool& computeSecond,
                          const bool& doNotSkip) {

    // if doNotSkip == true then the first clause is false and we will never return,
    // i.e. always compute, otherwise we skip DO_NOT_COMPUTE cells
    bool skipMoments = false;
    if (!doNotSkip && cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
        skipMoments = true;
    }

    // Clear old moments to zero value
    if (skipMoments == false) {
        cell->parameters[CellParams::RHOM  ] = 0.0;
        cell->parameters[CellParams::VX] = 0.0;
        cell->parameters[CellParams::VY] = 0.0;
        cell->parameters[CellParams::VZ] = 0.0;
        cell->parameters[CellParams::RHOQ  ] = 0.0;
        cell->parameters[CellParams::P_11] = 0.0;
        cell->parameters[CellParams::P_22] = 0.0;
        cell->parameters[CellParams::P_33] = 0.0;
    }

    // Loop over all particle species
    if (skipMoments == false) {
       for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
          vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
          if (blockContainer.size() == 0) continue;
          
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
          
          Population & pop = cell->get_population(popID);
          pop.RHO = array[0];
          pop.V[0] = divideIfNonZero(array[1], array[0]);
          pop.V[1] = divideIfNonZero(array[2], array[0]);
          pop.V[2] = divideIfNonZero(array[3], array[0]);
          
          // Store species' contribution to bulk velocity moments
          cell->parameters[CellParams::RHOM  ] += array[0]*mass;
          cell->parameters[CellParams::VX] += array[1]*mass;
          cell->parameters[CellParams::VY] += array[2]*mass;
          cell->parameters[CellParams::VZ] += array[3]*mass;
          cell->parameters[CellParams::RHOQ  ] += array[0]*charge;
       } // for-loop over particle species
       
       cell->parameters[CellParams::VX] = divideIfNonZero(cell->parameters[CellParams::VX], cell->parameters[CellParams::RHOM]);
       cell->parameters[CellParams::VY] = divideIfNonZero(cell->parameters[CellParams::VY], cell->parameters[CellParams::RHOM]);
       cell->parameters[CellParams::VZ] = divideIfNonZero(cell->parameters[CellParams::VZ], cell->parameters[CellParams::RHOM]);
    }

    // Compute second moments only if requested
    if (computeSecond == false) return;
            
    // Loop over all particle species
    for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
       vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
       if (blockContainer.size() == 0) continue;
       
       const Realf* data       = blockContainer.getData();
       const Real* blockParams = blockContainer.getParameters();
       const Real mass = getObjectWrapper().particleSpecies[popID].mass;
       
       // Temporary array for storing moments
       Real array[3];
       for (int i=0; i<3; ++i) array[i] = 0.0;

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
       pop.P[0] = mass*array[0];
       pop.P[1] = mass*array[1];
       pop.P[2] = mass*array[2];
       
       cell->parameters[CellParams::P_11] += pop.P[0];
       cell->parameters[CellParams::P_22] += pop.P[1];
       cell->parameters[CellParams::P_33] += pop.P[2];
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
        const bool& computeSecond) {
 
    phiprof::start("compute-moments-n");
    creal HALF = 0.5;

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
          }

          const Real dx = cell->parameters[CellParams::DX];
          const Real dy = cell->parameters[CellParams::DY];
          const Real dz = cell->parameters[CellParams::DZ];

          vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
          if (blockContainer.size() == 0) continue;
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
          Population & pop = cell->get_population(popID);
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
      phiprof::stop("compute-moments-n");
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
         if (blockContainer.size() == 0) continue;
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;

         // Temporary array where species' contribution to 2nd moments is accumulated
         Real array[3];
         for (int i=0; i<3; ++i) array[i] = 0.0;

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
         pop.P_R[0] = mass*array[0];
         pop.P_R[1] = mass*array[1];
         pop.P_R[2] = mass*array[2];
         
         cell->parameters[CellParams::P_11_R] += pop.P_R[0];
         cell->parameters[CellParams::P_22_R] += pop.P_R[1];
         cell->parameters[CellParams::P_33_R] += pop.P_R[2];
      } // for-loop over spatial cells
   } // for-loop over particle species

   phiprof::stop("compute-moments-n");
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
        const bool& computeSecond) {
 
   phiprof::start("Compute _V moments");
   
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
         }

         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         if (blockContainer.size() == 0) continue;
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
         Population & pop = cell->get_population(popID);
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
      phiprof::stop("Compute _V moments");
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
         if (blockContainer.size() == 0) continue;
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;

         // Temporary array where moments are stored
         Real array[3];
         for (int i=0; i<3; ++i) array[i] = 0.0;

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
         pop.P_V[0] = mass*array[0];
         pop.P_V[1] = mass*array[1];
         pop.P_V[2] = mass*array[2];
         
         cell->parameters[CellParams::P_11_V] += pop.P_V[0];
         cell->parameters[CellParams::P_22_V] += pop.P_V[1];
         cell->parameters[CellParams::P_33_V] += pop.P_V[2];
         
      } // for-loop over spatial cells
   } // for-loop over particle species

   phiprof::stop("Compute _V moments");
}
