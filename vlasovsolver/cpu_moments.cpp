/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

using namespace std;

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the 
 * given spatial cell. The calculated moments include contributions from 
 * all existing particle populations. This function is AMR safe.
 * @param cell Spatial cell.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param doNotSkip If false, DO_NOT_COMPUTE cells, or boundary cells of layer larger than 1, are skipped.*/
void calculateCellMoments(spatial_cell::SpatialCell* cell,
                          const bool& computeSecond,
                          const bool& doNotSkip) {

    // if doNotSkip == true then the first clause is false and we will never return,
    // i.e. always compute, otherwise we skip DO_NOT_COMPUTE cells
    // or boundary cells of layer larger than 1.
    bool skipMoments = false;
    if (!doNotSkip &&
        (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
        (cell->sysBoundaryLayer != 1  &&
         cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
        ) {
        skipMoments = true;
    }

    // Clear old moments to zero value
    if (skipMoments == false) {
        cell->parameters[CellParams::RHOM  ] = 0.0;
        cell->parameters[CellParams::RHOMVX] = 0.0;
        cell->parameters[CellParams::RHOMVY] = 0.0;
        cell->parameters[CellParams::RHOMVZ] = 0.0;
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
          const Real massRatio = getObjectWrapper().particleSpecies[popID].mass / physicalconstants::MASS_PROTON;
          for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
             blockVelocityFirstMoments(data+blockLID*WID3,
                                       blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                       massRatio,array);
          }
          
          Population & pop = cell->get_population(popID);
          pop.RHO = array[0];
          pop.RHOV[0] = array[1];
          pop.RHOV[1] = array[2];
          pop.RHOV[2] = array[3];
          
          // Store species' contribution to bulk velocity moments
          cell->parameters[CellParams::RHOM  ] += array[0]*mass;
          cell->parameters[CellParams::RHOMVX] += array[1]*mass;
          cell->parameters[CellParams::RHOMVY] += array[2]*mass;
          cell->parameters[CellParams::RHOMVZ] += array[3]*mass;
          cell->parameters[CellParams::RHOQ  ] += array[0]*charge;
       } // for-loop over particle species
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
                                     pop.RHO,
                                     pop.RHOV,
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
 * given spatial cell. Additionally, for each species, calculate the maximum 
 * spatial time step so that CFL(spatial)=1. The calculated moments include 
 * contributions from all existing particle populations. The calculated moments 
 * are stored to SpatialCell::parameters in _R variables. This function is AMR safe.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.*/
void calculateMoments_R_maxdt(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const std::vector<CellID>& cells,
        const bool& computeSecond) {
 
    phiprof::start("compute-moments-n-maxdt");
    creal HALF = 0.5;

    for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
       #pragma omp parallel for
       for (size_t c=0; c<cells.size(); ++c) {
          const CellID cellID = cells[c];
          SpatialCell* cell = mpiGrid[cells[c]];
          
          // Clear old moments to zero value
          if (popID == 0) {
             cell->parameters[CellParams::RHOM_R  ] = 0.0;
             cell->parameters[CellParams::RHOMVX_R] = 0.0;
             cell->parameters[CellParams::RHOMVY_R] = 0.0;
             cell->parameters[CellParams::RHOMVZ_R] = 0.0;
             cell->parameters[CellParams::RHOQ_R  ] = 0.0;
             cell->parameters[CellParams::P_11_R] = 0.0;
             cell->parameters[CellParams::P_22_R] = 0.0;
             cell->parameters[CellParams::P_33_R] = 0.0;
          }

          const Real dx = cell->parameters[CellParams::DX];
          const Real dy = cell->parameters[CellParams::DY];
          const Real dz = cell->parameters[CellParams::DZ];

          // Reset spatial max DT
          if (popID == 0) cell->parameters[CellParams::MAXRDT] = numeric_limits<Real>::max();
          cell->set_max_r_dt(popID,numeric_limits<Real>::max());

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
             // compute maximum dt. Algorithm has a CFL condition, since it
             // is written only for the case where we have a stencil
             // supporting max translation of one cell
             const Real EPS = numeric_limits<Real>::min()*1000;
             for (unsigned int i=0; i<WID;i+=WID-1) {
                const Real Vx 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::VXCRD] 
                  + (i+HALF)*blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVX]
                  + EPS;
                const Real Vy 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::VYCRD] 
                  + (i+HALF)*blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVY]
                  + EPS;
                    const Real Vz 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::VZCRD]
                  + (i+HALF)*blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVZ]
                  + EPS;

                const Real dt_max_cell = min(dx/fabs(Vx),min(dy/fabs(Vy),dz/fabs(Vz)));
                cell->parameters[CellParams::MAXRDT] = min(dt_max_cell,cell->parameters[CellParams::MAXRDT]);
                cell->set_max_r_dt(popID,min(dt_max_cell,cell->get_max_r_dt(popID)));
             }

             const Real massRatio = getObjectWrapper().particleSpecies[popID].mass / physicalconstants::MASS_PROTON;
             blockVelocityFirstMoments(data+blockLID*WID3,
                                       blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                       massRatio,array);
          } // for-loop over velocity blocks

          // Store species' contribution to bulk velocity moments
          Population & pop = cell->get_population(popID);
          pop.RHO_R = array[0];
          pop.RHOV_R[0] = array[1];
          pop.RHOV_R[1] = array[2];
          pop.RHOV_R[2] = array[3];
          
          cell->parameters[CellParams::RHOM_R  ] += array[0]*mass;
          cell->parameters[CellParams::RHOMVX_R] += array[1]*mass;
          cell->parameters[CellParams::RHOMVY_R] += array[2]*mass;
          cell->parameters[CellParams::RHOMVZ_R] += array[3]*mass;
          cell->parameters[CellParams::RHOQ_R  ] += array[0]*charge;
       } // for-loop over spatial cells
    } // for-loop over particle species

   // Compute second moments only if requested.
   if (computeSecond == false) {
      phiprof::stop("compute-moments-n-maxdt");
      return;
   }

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         SpatialCell* cell = mpiGrid[cells[c]];
       
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
                                       pop.RHO_R,
                                       pop.RHOV_R,
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

   phiprof::stop("compute-moments-n-maxdt");
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
         const CellID cellID = cells[c];
         SpatialCell* cell = mpiGrid[cells[c]];
         
         // Clear old moments to zero value
         if (popID == 0) {
             cell->parameters[CellParams::RHOM_V  ] = 0.0;
             cell->parameters[CellParams::RHOMVX_V] = 0.0;
             cell->parameters[CellParams::RHOMVY_V] = 0.0;
             cell->parameters[CellParams::RHOMVZ_V] = 0.0;
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

         const Real massRatio = getObjectWrapper().particleSpecies[popID].mass / physicalconstants::MASS_PROTON;

         // Calculate species' contribution to first velocity moments
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocityFirstMoments(data+blockLID*WID3,
                                      blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                      massRatio,array);
         }
         
         // Store species' contribution to bulk velocity moments
         Population & pop = cell->get_population(popID);
         pop.RHO_V = array[0];
         pop.RHOV_V[0] = array[1];
         pop.RHOV_V[1] = array[2];
         pop.RHOV_V[2] = array[3];
         
         cell->parameters[CellParams::RHOM_V  ] += array[0]*mass;
         cell->parameters[CellParams::RHOMVX_V] += array[1]*mass;
         cell->parameters[CellParams::RHOMVY_V] += array[2]*mass;
         cell->parameters[CellParams::RHOMVZ_V] += array[3]*mass;
         cell->parameters[CellParams::RHOQ_V  ] += array[0]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

   // Compute second moments only if requested
   if (computeSecond == false) {
      phiprof::stop("Compute _V moments");
      return;
   }

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         SpatialCell* cell = mpiGrid[cells[c]];

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
                                       pop.RHO_V,
                                       pop.RHOV_V,
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
