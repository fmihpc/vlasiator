/*
  This file is part of Vlasiator.

  Copyright 2010-2015 Finnish Meteorological Institute

*/

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdint.h>

#ifdef _OPENMP
   #include <omp.h>
#endif
#include <zoltan.h>
#include <dccrg.hpp>
#include <phiprof.hpp>

#include "../spatial_cell.hpp"
#include "../vlasovmover.h"
#include "../grid.h"
#include "../definitions.h"
#include "../object_wrapper.h"
#include "../mpiconversion.h"

#include "cpu_moments.h"
#include "cpu_acc_semilag.hpp"
#include "cpu_trans_map.hpp"

using namespace std;
using namespace spatial_cell;

creal ZERO    = 0.0;
creal HALF    = 0.5;
creal FOURTH  = 1.0/4.0;
creal SIXTH   = 1.0/6.0;
creal ONE     = 1.0;
creal TWO     = 2.0;
creal EPSILON = 1.0e-25;

/** Propagates the distribution function in spatial space. 
    
    Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
    three‐dimensional monotone and conservative semi‐Lagrangian scheme
    (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
    Meteorological Society 138.667 (2012): 1640-1651.
  
 */
void calculateSpatialTranslation(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const vector<CellID>& localCells,
        const vector<CellID>& local_propagated_cells,
        const vector<CellID>& local_target_cells,
        const vector<CellID>& remoteTargetCellsx,
        const vector<CellID>& remoteTargetCellsy,
        const vector<CellID>& remoteTargetCellsz,
        creal dt,
        const int& popID) {

    int trans_timer;
    bool localTargetGridGenerated = false;

    // ------------- SLICE - map dist function in Z --------------- //
   if(P::zcells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-z","MPI");
      phiprof::start(trans_timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      mpiGrid.start_remote_neighbor_copy_updates(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      // generate target grid in the temporary arrays, same size as
      // original one. We only need to create these in target cells
      createTargetGrid(mpiGrid,remoteTargetCellsz,popID);

      if(!localTargetGridGenerated){ 
         createTargetGrid(mpiGrid,local_target_cells,popID);
         localTargetGridGenerated=true;
      }

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      phiprof::start("compute-mapping-z");
      #pragma omp parallel
      {
         const int tid = omp_get_thread_num();
         no_subnormals();
         for (size_t c=0; c<local_propagated_cells.size(); ++c) {
            Real t_start = 0;
            if (tid == 0) if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();

            trans_map_1d(mpiGrid,local_propagated_cells[c], 2, dt,popID); // map along z//

            if (tid == 0) if (Parameters::prepareForRebalance == true) {
               mpiGrid[local_propagated_cells[c]]->get_cell_parameters()[CellParams::LBWEIGHTCOUNTER] 
                       += (MPI_Wtime()-t_start);
            }
         }
      }
      phiprof::stop("compute-mapping-z");

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(trans_timer);

      trans_timer=phiprof::initializeTimer("update_remote-z","MPI");
      phiprof::start("update_remote-z");
      update_remote_mapping_contribution(mpiGrid, 2,+1,popID);
      update_remote_mapping_contribution(mpiGrid, 2,-1,popID);
      phiprof::stop("update_remote-z");

      clearTargetGrid(mpiGrid,remoteTargetCellsz);
      swapTargetSourceGrid(mpiGrid, local_target_cells,popID);
      zeroTargetGrid(mpiGrid, local_target_cells);
   }

   // ------------- SLICE - map dist function in X --------------- //
   if(P::xcells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-x","MPI");
      phiprof::start(trans_timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      mpiGrid.start_remote_neighbor_copy_updates(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      createTargetGrid(mpiGrid,remoteTargetCellsx,popID);
       if(!localTargetGridGenerated){ 
         createTargetGrid(mpiGrid,local_target_cells,popID);
         localTargetGridGenerated=true;
      }

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);

      phiprof::start("compute-mapping-x");
      #pragma omp parallel
      {
         const int tid = omp_get_thread_num();
         no_subnormals();
         for (size_t c=0; c<local_propagated_cells.size(); ++c) {
            Real t_start = 0;
            if (tid == 0) if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();

            trans_map_1d(mpiGrid,local_propagated_cells[c],0,dt,popID); // map along x//

            if (tid == 0) if (Parameters::prepareForRebalance == true) {
               mpiGrid[local_propagated_cells[c]]->get_cell_parameters()[CellParams::LBWEIGHTCOUNTER] 
                       += (MPI_Wtime()-t_start);
            }
         }
      }
      phiprof::stop("compute-mapping-x");

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(trans_timer);

      trans_timer=phiprof::initializeTimer("update_remote-x","MPI");
      phiprof::start("update_remote-x");
      update_remote_mapping_contribution(mpiGrid, 0,+1,popID);
      update_remote_mapping_contribution(mpiGrid, 0,-1,popID);
      phiprof::stop("update_remote-x");
      clearTargetGrid(mpiGrid,remoteTargetCellsx);
      swapTargetSourceGrid(mpiGrid, local_target_cells,popID);
      zeroTargetGrid(mpiGrid, local_target_cells);
   }
   
   // ------------- SLICE - map dist function in Y --------------- //
   if(P::ycells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-y","MPI");
      phiprof::start(trans_timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      mpiGrid.start_remote_neighbor_copy_updates(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      createTargetGrid(mpiGrid,remoteTargetCellsy,popID);
      if(!localTargetGridGenerated){ 
         createTargetGrid(mpiGrid,local_target_cells,popID);
         localTargetGridGenerated=true;
      }
      
      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);

      phiprof::start("compute-mapping-y");
      #pragma omp parallel
      {
         const int tid = omp_get_thread_num();
         no_subnormals();
         for (size_t c=0; c<local_propagated_cells.size(); ++c) {
            Real t_start = 0;
            if (tid == 0) if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();
            
            trans_map_1d(mpiGrid,local_propagated_cells[c],1,dt,popID); // map along y//
            
            if (tid == 0) if (Parameters::prepareForRebalance == true) {
               mpiGrid[local_propagated_cells[c]]->get_cell_parameters()[CellParams::LBWEIGHTCOUNTER] 
                       += (MPI_Wtime()-t_start);
            }
         }
      }
      phiprof::stop("compute-mapping-y");

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(trans_timer);
      
      trans_timer=phiprof::initializeTimer("update_remote-y","MPI");
      phiprof::start("update_remote-y");
      update_remote_mapping_contribution(mpiGrid, 1,+1,popID);
      update_remote_mapping_contribution(mpiGrid, 1,-1,popID);
      phiprof::stop("update_remote-y");
      clearTargetGrid(mpiGrid,remoteTargetCellsy);
      swapTargetSourceGrid(mpiGrid, local_target_cells,popID);
   }

   clearTargetGrid(mpiGrid,local_target_cells);
}

/*!
  
  Propagates the distribution function in spatial space. 
  
  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

*/
void calculateSpatialTranslation(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        creal dt) {
   typedef Parameters P;
   
   phiprof::start("semilag-trans");

   const vector<CellID>& localCells = getLocalCells();
   vector<CellID> remoteTargetCellsx;
   vector<CellID> remoteTargetCellsy;
   vector<CellID> remoteTargetCellsz;
   vector<CellID> local_propagated_cells;
   vector<CellID> local_target_cells;
   
   // If dt=0 we are either initializing or distribution functions are not translated. 
   // In both cases go to the end of this function and calculate the moments.
   if (dt == 0.0) goto momentCalculation;
   
    phiprof::start("compute_cell_lists");
    remoteTargetCellsx = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID);
    remoteTargetCellsy = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID);
    remoteTargetCellsz = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID);

    // Figure out which spatial cells are translated, 
    // result independent of particle species.
    for (size_t c=0; c<localCells.size(); ++c) {
       if (do_translate_cell(mpiGrid[localCells[c]])) {
          local_propagated_cells.push_back(localCells[c]);
       }
    }

   // Figure out target spatial cells, result
   // independent of particle species.
   for (size_t c=0; c<localCells.size(); ++c) {
      if (mpiGrid[localCells[c]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         local_target_cells.push_back(localCells[c]);
      }
   }
   phiprof::stop("compute_cell_lists");

   // Translate all particle species
   for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      string profName = "translate "+getObjectWrapper().particleSpecies[popID].name;
      phiprof::start(profName);
      SpatialCell::setCommunicatedSpecies(popID);
      calculateSpatialTranslation(mpiGrid,localCells,local_propagated_cells,
                                  local_target_cells,remoteTargetCellsx,remoteTargetCellsy,
                                  remoteTargetCellsz,dt,popID);
      phiprof::stop(profName);
   }

   // Mapping complete, update moments and maximum dt limits //
momentCalculation:
   calculateMoments_R_maxdt(mpiGrid,localCells,true);

   Real minDT = 1e300;
   for (size_t c=0; c<localCells.size(); ++c) {
      if (mpiGrid[localCells[c]]->parameters[CellParams::MAXRDT] < minDT) 
         minDT = mpiGrid[localCells[c]]->parameters[CellParams::MAXRDT];
   }
   phiprof::stop("semilag-trans");
}

/*
  --------------------------------------------------
  Acceleration (velocity space propagation)
  --------------------------------------------------
*/

int getAccerelationSubcycles(SpatialCell* sc,Real dt,const int& popID) {
   return max( convert<int>(ceil(dt / sc->get_max_v_dt(popID))), 1);
}

/** Accelerate the given population to new time t+dt.
 * This function is AMR safe.
 * @param popID Particle population ID.
 * @param globalMaxSubcycles Number of times acceleration is subcycled.
 * @param step The current subcycle step.
 * @param mpiGrid Parallel grid library.
 * @param propagatedCells List of cells in which the population is accelerated.
 * @param dt Timestep.*/
void calculateAcceleration(const int& popID,const int& globalMaxSubcycles,const uint& step,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& propagatedCells,
                           const Real& dt) {
   // Set active population
   SpatialCell::setCommunicatedSpecies(popID);

   // Semi-Lagrangian acceleration for those cells which are subcycled
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t c=0; c<propagatedCells.size(); ++c) {
      const CellID cellID = propagatedCells[c];
      const Real maxVdt = mpiGrid[cellID]->get_max_v_dt(popID);

      //compute subcycle dt. The length is maxVdt on all steps
      //except the last one. This is to keep the neighboring
      //spatial cells in sync, so that two neighboring cells with
      //different number of subcycles have similar timesteps,
      //except that one takes an additional short step. This keeps
      //spatial block neighbors as much in sync as possible for
      //adjust blocks.
      Real subcycleDt;
      if( (step + 1) * maxVdt > dt) {
         subcycleDt = dt - step * maxVdt;
      } else{
         subcycleDt = maxVdt;
      }

      //generate pseudo-random order which is always the same irrespective of parallelization, restarts, etc
      char rngStateBuffer[256];
      random_data rngDataBuffer;

      // set seed, initialise generator and get value
      memset(&(rngDataBuffer), 0, sizeof(rngDataBuffer));
      #ifdef _AIX
         initstate_r(P::tstep + cellID, &(rngStateBuffer[0]), 256, NULL, &(rngDataBuffer));
         int64_t rndInt;
         random_r(&rndInt, &rngDataBuffer);
      #else
         initstate_r(P::tstep + cellID, &(rngStateBuffer[0]), 256, &(rngDataBuffer));
         int32_t rndInt;
         random_r(&rngDataBuffer, &rndInt);
      #endif
         
      uint map_order=rndInt%3;
      phiprof::start("cell-semilag-acc");
      cpu_accelerate_cell(mpiGrid[cellID],popID,map_order,subcycleDt);
      phiprof::stop("cell-semilag-acc");
   }

   //global adjust after each subcycle to keep number of blocks managable. Even the ones not
   //accelerating anyore participate. It is important to keep
   //the spatial dimension to make sure that we do not loose
   //stuff streaming in from other cells, perhaps not connected
   //to the existing distribution function in the cell.
   //- All cells update and communicate their lists of content blocks
   //- Only cells which were accerelated on this step need to be adjusted (blocks removed or added).
   //- Not done here on last step (done after loop)
   if(step < (globalMaxSubcycles - 1)) adjustVelocityBlocks(mpiGrid, propagatedCells, false, popID);
}

/** Accelerate all particle populations to new time t+dt. 
 * This function is AMR safe.
 * @param mpiGrid Parallel grid library.
 * @param dt Time step.*/
void calculateAcceleration(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           Real dt
                          ) {    
   typedef Parameters P;
   const vector<CellID>& cells = getLocalCells();

   if (dt == 0.0 && P::tstep > 0) {
      // Even if acceleration is turned off we need to adjust velocity blocks 
      // because the boundary conditions may have altered the velocity space, 
      // and to update changes in no-content blocks during translation.
      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
        adjustVelocityBlocks(mpiGrid, cells, true, popID);
      goto momentCalculation;
   }
   phiprof::start("semilag-acc");
    
    // Calculate first velocity moments, these are needed to 
    // calculate the transforms used in the accelerations.
    // Calculated moments are stored in the "_V" variables.
    calculateMoments_V(mpiGrid,cells,false);
    
    // Accelerate all particle species
    for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
       // Set active population
       SpatialCell::setCommunicatedSpecies(popID);

       // Iterate through all local cells and collect cells to propagate.
       // Ghost cells (spatial cells at the boundary of the simulation 
       // volume) do not need to be propagated:
       vector<CellID> propagatedCells;
       for (size_t c=0; c<cells.size(); ++c) {
          SpatialCell* SC = mpiGrid[cells[c]];
          const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = SC->get_velocity_mesh(popID);
          // disregard boundary cells and do not propagate spatial 
          // cells with no blocks (well, do not computes in practice)
          if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY && vmesh.size() != 0) {
             propagatedCells.push_back(cells[c]);
          }
       }

       // Compute global maximum for number of subcycles (collective operation).
       int maxSubcycles=0;
       int globalMaxSubcycles;
       for (size_t c=0; c<propagatedCells.size(); ++c) {
          const CellID cellID = propagatedCells[c];
          int subcycles = getAccerelationSubcycles(mpiGrid[cellID],dt,popID);
          mpiGrid[cellID]->parameters[CellParams::ACCSUBCYCLES] = subcycles;
          maxSubcycles=maxSubcycles < subcycles ? subcycles:maxSubcycles;
       }
       MPI_Allreduce(&maxSubcycles, &globalMaxSubcycles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

       // substep global max times
       for(uint step=0; step<globalMaxSubcycles; ++step) {
          // prune list of cells to propagate to only contained those which are now subcycled
          vector<CellID> temp;
          for (size_t c=0; c<propagatedCells.size(); ++c) {
             if (step < getAccerelationSubcycles(mpiGrid[propagatedCells[c]],dt,popID)) {
                temp.push_back(propagatedCells[c]);
             }
          }
          propagatedCells.swap(temp);
          calculateAcceleration(popID,globalMaxSubcycles,step,mpiGrid,propagatedCells,dt);
       } // for-loop over acceleration substeps
       
       // final adjust for all cells, also fixing remote cells.
       adjustVelocityBlocks(mpiGrid, cells, true, popID);
    } // for-loop over particle species

    phiprof::stop("semilag-acc");

   // Recalculate "_V" velocity moments
momentCalculation:
   calculateMoments_V(mpiGrid,cells,true);

   // Set CellParams::MAXVDT to be the minimum dt of all per-species values
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* cell = mpiGrid[cells[c]];
      cell->parameters[CellParams::MAXVDT] = numeric_limits<Real>::max();
      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         cell->parameters[CellParams::MAXVDT]
           = min(cell->get_max_v_dt(popID), cell->parameters[CellParams::MAXVDT]);
      }
   }
}

/*--------------------------------------------------
  Functions for computing moments
  --------------------------------------------------*/

void calculateInterpolatedVelocityMoments(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const int cp_rho,
   const int cp_rhovx,
   const int cp_rhovy,
   const int cp_rhovz,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33
) {
   const vector<CellID>& cells = getLocalCells();
   
   //Iterate through all local cells (excl. system boundary cells):
    #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      if(SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SC->parameters[cp_rho  ] = 0.5* ( SC->parameters[CellParams::RHO_R] + SC->parameters[CellParams::RHO_V] );
         SC->parameters[cp_rhovx] = 0.5* ( SC->parameters[CellParams::RHOVX_R] + SC->parameters[CellParams::RHOVX_V] );
         SC->parameters[cp_rhovy] = 0.5* ( SC->parameters[CellParams::RHOVY_R] + SC->parameters[CellParams::RHOVY_V] );
         SC->parameters[cp_rhovz] = 0.5* ( SC->parameters[CellParams::RHOVZ_R] + SC->parameters[CellParams::RHOVZ_V] );
         SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R] + SC->parameters[CellParams::P_11_V] );
         SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R] + SC->parameters[CellParams::P_22_V] );
         SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R] + SC->parameters[CellParams::P_33_V] );
      }
   }
}

void calculateInitialVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const vector<CellID>& cells = getLocalCells();
   phiprof::start("Calculate moments");

   // Iterate through all local cells (incl. system boundary cells):
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      calculateCellMoments(SC,true);

      // WARNING the following is sane as this function is only called by initializeGrid.
      // We need initialized _DT2 values for the dt=0 field propagation done in the beginning.
      // Later these will be set properly.
      SC->parameters[CellParams::RHO_DT2] = SC->parameters[CellParams::RHO];
      SC->parameters[CellParams::RHOVX_DT2] = SC->parameters[CellParams::RHOVX];
      SC->parameters[CellParams::RHOVY_DT2] = SC->parameters[CellParams::RHOVY];
      SC->parameters[CellParams::RHOVZ_DT2] = SC->parameters[CellParams::RHOVZ];
      SC->parameters[CellParams::P_11_DT2] = SC->parameters[CellParams::P_11];
      SC->parameters[CellParams::P_22_DT2] = SC->parameters[CellParams::P_22];
      SC->parameters[CellParams::P_33_DT2] = SC->parameters[CellParams::P_33];
   } // for-loop over spatial cells
   phiprof::stop("Calculate moments"); 
}
