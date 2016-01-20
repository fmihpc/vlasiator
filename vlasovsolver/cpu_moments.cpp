/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute

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
        cell->parameters[CellParams::RHO  ] = 0.0;
        cell->parameters[CellParams::RHOVX] = 0.0;
        cell->parameters[CellParams::RHOVY] = 0.0;
        cell->parameters[CellParams::RHOVZ] = 0.0;
        cell->parameters[CellParams::P_11] = 0.0;
        cell->parameters[CellParams::P_22] = 0.0;
        cell->parameters[CellParams::P_33] = 0.0;
    }

    // Loop over all particle species
    if (skipMoments == false) {
       for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
          vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
          if (blockContainer.size() == 0) continue;
          
          const Realf* data       = blockContainer.getData();
          const Real* blockParams = blockContainer.getParameters();
          
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
          
          // Store species' contribution to bulk velocity moments
          cell->parameters[CellParams::RHO  ] += array[0];
          cell->parameters[CellParams::RHOVX] += array[1];
          cell->parameters[CellParams::RHOVY] += array[2];
          cell->parameters[CellParams::RHOVZ] += array[3];
       } // for-loop over particle species
    }

    // Compute second moments only if requested
    if (computeSecond == false) return;
            
    // Loop over all particle species
    for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
       vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
       if (blockContainer.size() == 0) continue;
       
       const Realf* data       = blockContainer.getData();
       const Real* blockParams = blockContainer.getParameters();
       
       // Temporary array for storing moments
       Real array[3];
       for (int i=0; i<3; ++i) array[i] = 0.0;

       // Calculate species' contribution to second velocity moments
       for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
          blockVelocitySecondMoments(data+blockLID*WID3,
                                     blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                     cell->parameters,
                                     CellParams::RHO,
                                     CellParams::RHOVX,
                                     CellParams::RHOVY,
                                     CellParams::RHOVZ,
                                     array);
       }

       // Store species' contribution to bulk velocity moments
       const Real mass = getObjectWrapper().particleSpecies[popID].mass;
       cell->parameters[CellParams::P_11] += mass*array[0];
       cell->parameters[CellParams::P_22] += mass*array[1];
       cell->parameters[CellParams::P_33] += mass*array[2];
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

    for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
       #pragma omp parallel for
       for (size_t c=0; c<cells.size(); ++c) {
          const CellID cellID = cells[c];
          SpatialCell* cell = mpiGrid[cells[c]];

          // Clear old moments to zero value
          if (popID == 0) {
             cell->parameters[CellParams::RHO_R  ] = 0.0;
             cell->parameters[CellParams::RHOVX_R] = 0.0;
             cell->parameters[CellParams::RHOVY_R] = 0.0;
             cell->parameters[CellParams::RHOVZ_R] = 0.0;
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
          cell->parameters[CellParams::RHO_R  ] += array[0];
          cell->parameters[CellParams::RHOVX_R] += array[1];
          cell->parameters[CellParams::RHOVY_R] += array[2];
          cell->parameters[CellParams::RHOVZ_R] += array[3];
       } // for-loop over spatial cells
    } // for-loop over particle species

   // Compute second moments only if requested.
   if (computeSecond == false) {
      phiprof::stop("compute-moments-n-maxdt");
      return;
   }

   for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         SpatialCell* cell = mpiGrid[cells[c]];
       
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         if (blockContainer.size() == 0) continue;
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();

         // Temporary array where species' contribution to 2nd moments is accumulated
         Real array[3];
         for (int i=0; i<3; ++i) array[i] = 0.0;

         // Calculate species' contribution to second velocity moments
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocitySecondMoments(data+blockLID*WID3,
                                       blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                       cell->parameters,
                                       CellParams::RHO_R,
                                       CellParams::RHOVX_R,
                                       CellParams::RHOVY_R,
                                       CellParams::RHOVZ_R,
                                       array);
         } // for-loop over velocity blocks

         // Store species' contribution to 2nd bulk velocity moments
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         cell->parameters[CellParams::P_11_R] += mass*array[0];
         cell->parameters[CellParams::P_22_R] += mass*array[1];
         cell->parameters[CellParams::P_33_R] += mass*array[2];
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
   for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         SpatialCell* cell = mpiGrid[cells[c]];

         // Clear old moments to zero value
         if (popID == 0) {
            cell->parameters[CellParams::RHO_V  ] = 0.0;
            cell->parameters[CellParams::RHOVX_V] = 0.0;
            cell->parameters[CellParams::RHOVY_V] = 0.0;
            cell->parameters[CellParams::RHOVZ_V] = 0.0;
            cell->parameters[CellParams::P_11_V] = 0.0;
            cell->parameters[CellParams::P_22_V] = 0.0;
            cell->parameters[CellParams::P_33_V] = 0.0;
         }

         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         if (blockContainer.size() == 0) continue;
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();

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
         cell->parameters[CellParams::RHO_V  ] += array[0];
         cell->parameters[CellParams::RHOVX_V] += array[1];
         cell->parameters[CellParams::RHOVY_V] += array[2];
         cell->parameters[CellParams::RHOVZ_V] += array[3];         
      } // for-loop over spatial cells
   } // for-loop over particle species

   // Compute second moments only if requested
   if (computeSecond == false) {
      phiprof::stop("Compute _V moments");
      return;
   }

   for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         SpatialCell* cell = mpiGrid[cells[c]];

         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
         if (blockContainer.size() == 0) continue;
         const Realf* data       = blockContainer.getData();
         const Real* blockParams = blockContainer.getParameters();

         // Temporary array where moments are stored
         Real array[3];
         for (int i=0; i<3; ++i) array[i] = 0.0;

         // Calculate species' contribution to second velocity moments
         for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
            blockVelocitySecondMoments(
                                       data+blockLID*WID3,
                                       blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                       cell->parameters,
                                       CellParams::RHO_V,
                                       CellParams::RHOVX_V,
                                       CellParams::RHOVY_V,
                                       CellParams::RHOVZ_V,
                                       array);
         } // for-loop over velocity blocks
         
         // Store species' contribution to 2nd bulk velocity moments
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         cell->parameters[CellParams::P_11_V] += mass*array[0];
         cell->parameters[CellParams::P_22_V] += mass*array[1];
         cell->parameters[CellParams::P_33_V] += mass*array[2];
      } // for-loop over spatial cells
   } // for-loop over particle species

   phiprof::stop("Compute _V moments");
}
