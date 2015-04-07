/*
  This file is part of Vlasiator.

  Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include "omp.h"
#endif
#include <zoltan.h>

#include "../vlasovmover.h"
#include "phiprof.hpp"

#include "cpu_moments.h"
#include "cpu_acc_semilag.hpp"
#include "cpu_trans_map.hpp"

#include <stdint.h>
#include <dccrg.hpp>

#include "spatial_cell.hpp"
#include "../grid.h"
#include "../definitions.h"

using namespace std;
using namespace spatial_cell;

creal ZERO    = 0.0;
creal HALF    = 0.5;
creal FOURTH  = 1.0/4.0;
creal SIXTH   = 1.0/6.0;
creal ONE     = 1.0;
creal TWO     = 2.0;
creal EPSILON = 1.0e-25;


/*!
  
  Propagates the distribution function in spatial space. 
  
  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

*/

void calculateSpatialTranslation(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   creal dt
) {
   const size_t popID = 0;
   typedef Parameters P;
   int trans_timer;
   bool localTargetGridGenerated = false;
   const vector<CellID>& localCells = getLocalCells();
   vector<CellID> remoteTargetCellsx;
   vector<CellID> remoteTargetCellsy;
   vector<CellID> remoteTargetCellsz;
   vector<CellID> local_propagated_cells;
   vector<CellID> local_target_cells;

   phiprof::start("semilag-trans");

   // If dt=0 we are either initializing or distribution functions are not propagated. 
   // In both cases go to the end of this function and calculate the moments.
   if (dt == 0.0) goto momentCalculation;

   phiprof::start("compute_cell_lists");
   remoteTargetCellsx = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID);
   remoteTargetCellsy = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID);
   remoteTargetCellsz = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID);

   for (size_t c=0; c<localCells.size(); ++c) {
      if(do_translate_cell(mpiGrid[localCells[c]])){
         local_propagated_cells.push_back(localCells[c]);
      }
   }
   for (size_t c=0; c<localCells.size(); ++c) {
      if(mpiGrid[localCells[c]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         local_target_cells.push_back(localCells[c]);
      }
   }
   phiprof::stop("compute_cell_lists");

   // ------------- SLICE - map dist function in Z --------------- //
   if(P::zcells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-z","MPI");
      phiprof::start(trans_timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      mpiGrid.start_remote_neighbor_copy_updates(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      /*generate target grid in the temporary arrays, same size as
       *   original one. We only need to create these in target cells*/
      createTargetGrid(mpiGrid,remoteTargetCellsz);

      if(!localTargetGridGenerated){ 
         createTargetGrid(mpiGrid,local_target_cells);
         localTargetGridGenerated=true;
      }
      
      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      phiprof::start("compute-mapping-z");
      #pragma omp parallel
      {
         no_subnormals();
         for (size_t c=0; c<local_propagated_cells.size(); ++c) {
            trans_map_1d(mpiGrid,local_propagated_cells[c], 2, dt); // map along z//
         }
      }
      phiprof::stop("compute-mapping-z");

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(trans_timer);
      
      trans_timer=phiprof::initializeTimer("update_remote-z","MPI");
      phiprof::start("update_remote-z");
      update_remote_mapping_contribution(mpiGrid, 2, 1);
      update_remote_mapping_contribution(mpiGrid, 2, -1);
      phiprof::stop("update_remote-z");

      clearTargetGrid(mpiGrid,remoteTargetCellsz);
      swapTargetSourceGrid(mpiGrid, local_target_cells);
      zeroTargetGrid(mpiGrid, local_target_cells);
   }

   // ------------- SLICE - map dist function in X --------------- //
   if(P::xcells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-x","MPI");
      phiprof::start(trans_timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      mpiGrid.start_remote_neighbor_copy_updates(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      createTargetGrid(mpiGrid,remoteTargetCellsx);
       if(!localTargetGridGenerated){ 
         createTargetGrid(mpiGrid,local_target_cells);
         localTargetGridGenerated=true;
      }
       
      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);

      phiprof::start("compute-mapping-x");
      #pragma omp parallel
      {
         no_subnormals();
         for (size_t c=0; c<local_propagated_cells.size(); ++c) {
            trans_map_1d(mpiGrid,local_propagated_cells[c], 0, dt); // map along x//
         }
      }
      phiprof::stop("compute-mapping-x");

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(trans_timer);
      
      trans_timer=phiprof::initializeTimer("update_remote-x","MPI");
      phiprof::start("update_remote-x");
      update_remote_mapping_contribution(mpiGrid, 0, 1);
      update_remote_mapping_contribution(mpiGrid, 0, -1);
      phiprof::stop("update_remote-x");
      clearTargetGrid(mpiGrid,remoteTargetCellsx);
      swapTargetSourceGrid(mpiGrid, local_target_cells);
      zeroTargetGrid(mpiGrid, local_target_cells);

   }
   
   // ------------- SLICE - map dist function in Y --------------- //
   if(P::ycells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-y","MPI");
      phiprof::start(trans_timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      mpiGrid.start_remote_neighbor_copy_updates(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
      
      createTargetGrid(mpiGrid,remoteTargetCellsy);
      if(!localTargetGridGenerated){ 
         createTargetGrid(mpiGrid,local_target_cells);
         localTargetGridGenerated=true;
      }
      
      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);

      phiprof::start("compute-mapping-y");
      #pragma omp parallel
      {
         no_subnormals();
         for (size_t c=0; c<local_propagated_cells.size(); ++c) {
            trans_map_1d(mpiGrid,local_propagated_cells[c], 1, dt); // map along y//
         }
      }
      phiprof::stop("compute-mapping-y");

      phiprof::start(trans_timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(trans_timer);
      
      trans_timer=phiprof::initializeTimer("update_remote-y","MPI");
      phiprof::start("update_remote-y");
      update_remote_mapping_contribution(mpiGrid, 1, 1);
      update_remote_mapping_contribution(mpiGrid, 1, -1);
      phiprof::stop("update_remote-y");
      clearTargetGrid(mpiGrid,remoteTargetCellsy);
      swapTargetSourceGrid(mpiGrid, local_target_cells);
   }
   
   clearTargetGrid(mpiGrid,local_target_cells);

momentCalculation:
   
   // Mapping complete, update moments //
   phiprof::start("compute-moments-n-maxdt");
   // Note: Parallelization over blocks is not thread-safe
   #pragma omp parallel for
   for (size_t c=0; c<localCells.size(); ++c) {
      SpatialCell* SC=mpiGrid[localCells[c]];

      const Real dx=SC->parameters[CellParams::DX];
      const Real dy=SC->parameters[CellParams::DY];
      const Real dz=SC->parameters[CellParams::DZ];
      SC->parameters[CellParams::RHO_R  ] = 0.0;
      SC->parameters[CellParams::RHOVX_R] = 0.0;
      SC->parameters[CellParams::RHOVY_R] = 0.0;
      SC->parameters[CellParams::RHOVZ_R] = 0.0;
      SC->parameters[CellParams::P_11_R ] = 0.0;
      SC->parameters[CellParams::P_22_R ] = 0.0;
      SC->parameters[CellParams::P_33_R ] = 0.0;

      //Reset spatial max DT
      SC->parameters[CellParams::MAXRDT]=numeric_limits<Real>::max();
      for (vmesh::LocalID block_i=0; block_i<SC->get_number_of_velocity_blocks(); ++block_i) {
         const Real* const blockParams = SC->get_block_parameters(block_i);

         //compute maximum dt. Algorithm has a CFL condition, since it
         //is written only for the case where we have a stencil
         //supporting max translation of one cell
         for (unsigned int i=0; i<WID;i+=WID-1) {
            const Real Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
            const Real Vy = blockParams[BlockParams::VYCRD] + (i+HALF)*blockParams[BlockParams::DVY];
            const Real Vz = blockParams[BlockParams::VZCRD] + (i+HALF)*blockParams[BlockParams::DVZ];
            
            if(fabs(Vx)!=ZERO) SC->parameters[CellParams::MAXRDT]=min(dx/fabs(Vx),SC->parameters[CellParams::MAXRDT]);
            if(fabs(Vy)!=ZERO) SC->parameters[CellParams::MAXRDT]=min(dy/fabs(Vy),SC->parameters[CellParams::MAXRDT]);
            if(fabs(Vz)!=ZERO) SC->parameters[CellParams::MAXRDT]=min(dz/fabs(Vz),SC->parameters[CellParams::MAXRDT]);
         }
         
         //compute first moments for this block
         if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
            cpu_calcVelocityFirstMoments(
               SC,
               block_i,
               CellParams::RHO_R,
               CellParams::RHOVX_R,
               CellParams::RHOVY_R,
               CellParams::RHOVZ_R
            );   //set first moments after translation
      }
      // Second iteration needed as rho has to be already computed when computing pressure
      for (vmesh::LocalID block_i=0; block_i< SC->get_number_of_velocity_blocks(); ++block_i){
         //compute second moments for this block
         if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
            cpu_calcVelocitySecondMoments(
               SC,
               block_i,
               CellParams::RHO_R,
               CellParams::RHOVX_R,
               CellParams::RHOVY_R,
               CellParams::RHOVZ_R,
               CellParams::P_11_R,
               CellParams::P_22_R,
               CellParams::P_33_R
            );   //set second moments after translation
      }
   }
   phiprof::stop("compute-moments-n-maxdt");
   phiprof::stop("semilag-trans");
}

/*
  --------------------------------------------------
  Acceleration (velocity space propagation)
  --------------------------------------------------
*/

int getAccerelationSubcycles(SpatialCell* sc, Real dt){
   return max(convert<int>(ceil(dt / sc->parameters[CellParams::MAXVDT])),1);
}

void calculateAcceleration(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   Real dt
) {
   typedef Parameters P;
   const vector<CellID> cells = getLocalCells();
   vector<CellID> propagatedCells;
   int maxSubcycles=0;
   int globalMaxSubcycles;

   if (dt == 0.0 && P::tstep > 0) goto momentCalculation;
   
//    if(dt > 0)  // FIXME this has to be deactivated to support regular projects but it breaks test_trans support most likely, with this on dt stays 0
   phiprof::start("semilag-acc");

   // Iterate through all local cells and collect cells to propagate.
   // Ghost cells (spatial cells at the boundary of the simulation 
   // volume) do not need to be propagated:
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* SC = mpiGrid[cells[c]];
      //disregard boundary cells
      //do not integrate cells with no blocks  (well, do not computes in practice)
      if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY &&
          SC->get_number_of_velocity_blocks() != 0) {
         propagatedCells.push_back(cells[c]);
      }
   }

   //Compute global maximum for number of subcycles (collective operation)
   for (size_t c=0; c<propagatedCells.size(); ++c) {
      const CellID cellID = propagatedCells[c];
      int subcycles = getAccerelationSubcycles(mpiGrid[cellID], dt);
      mpiGrid[cellID]->parameters[CellParams::ACCSUBCYCLES] = subcycles;
      maxSubcycles=maxSubcycles < subcycles ? subcycles:maxSubcycles;
   }
   MPI_Allreduce(&maxSubcycles, &globalMaxSubcycles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

   //substep global max times
   for(uint step = 0;step < globalMaxSubcycles; step ++) {

      //prune list of cells to propagate to only contained those which are now subcycled
      vector<CellID> temp;
      for (size_t c=0; c<propagatedCells.size(); ++c) {
         if(step < getAccerelationSubcycles(mpiGrid[propagatedCells[c]], dt)) {
            temp.push_back(propagatedCells[c]);
         }
      }
      propagatedCells.swap(temp);

      //Semilagrangian acceleration for those cells which are subcycled
#pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<propagatedCells.size(); ++c) {
         no_subnormals();
         const CellID cellID = propagatedCells[c];
         const Real maxVdt = mpiGrid[cellID]->parameters[CellParams::MAXVDT]; 

         //compute subcycle dt. The length is maVdt on all steps
         //except the last one. This is to keep the neighboring
         //spatial cells in sync, so that two neighboring cells with
         //different number of subcycles have similar timesteps,
         //except that one takes an additional short step. This keeps
         //spatial block neighbors as much in sync as possible for
         //adjust blocks.
         Real subcycleDt;
         if( (step + 1) * maxVdt > dt) {
            subcycleDt = dt - step * maxVdt;
         }
         else{
            subcycleDt = maxVdt;
         }
         //generate pseudo-random order which is always the same irrespectiive of parallelization, restarts, etc
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
         cpu_accelerate_cell(mpiGrid[cellID],map_order,subcycleDt);
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
      if(step < (globalMaxSubcycles - 1))
         adjustVelocityBlocks(mpiGrid, propagatedCells, false);
   }
   //final adjust for all cells, also fixing remote cells.
   adjustVelocityBlocks(mpiGrid, cells, true);
   phiprof::stop("semilag-acc");   

momentCalculation:

   phiprof::start("Compute moments");
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      //compute moments after acceleration
      mpiGrid[cellID]->parameters[CellParams::RHO_V  ] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::RHOVX_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::RHOVY_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::RHOVZ_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::P_11_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::P_22_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::P_33_V] = 0.0;
      
      for (vmesh::LocalID block_i=0; block_i<mpiGrid[cellID]->get_number_of_velocity_blocks(); ++block_i) {
         cpu_calcVelocityFirstMoments(
            mpiGrid[cellID],
            block_i,
            CellParams::RHO_V,
            CellParams::RHOVX_V,
            CellParams::RHOVY_V,
            CellParams::RHOVZ_V
                                      );   //set first moments after acceleration
      }

      // Second iteration needed as rho has to be already computed when computing pressure
      for (vmesh::LocalID block_i=0; block_i<mpiGrid[cellID]->get_number_of_velocity_blocks(); ++block_i) {
         cpu_calcVelocitySecondMoments(
            mpiGrid[cellID],
            block_i,
            CellParams::RHO_V,
            CellParams::RHOVX_V,
            CellParams::RHOVY_V,
            CellParams::RHOVZ_V,
            CellParams::P_11_V,
            CellParams::P_22_V,
            CellParams::P_33_V);   //set second moments after acceleration
      }
   }
   phiprof::stop("Compute moments");
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




void calculateCellVelocityMoments(
   SpatialCell* SC,
   bool doNotSkip // default: false
) {
   // if doNotSkip == true then the first clause is false and we will never return, i.e. always compute
   // otherwise we skip DO_NOT_COMPUTE cells
   // or boundary cells of layer larger than 1
   if (!doNotSkip &&
       (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       (SC->sysBoundaryLayer != 1  &&
       SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
      ) return;

   SC->parameters[CellParams::RHO  ] = 0.0;
   SC->parameters[CellParams::RHOVX] = 0.0;
   SC->parameters[CellParams::RHOVY] = 0.0;
   SC->parameters[CellParams::RHOVZ] = 0.0;
   SC->parameters[CellParams::P_11 ] = 0.0;
   SC->parameters[CellParams::P_22 ] = 0.0;
   SC->parameters[CellParams::P_33 ] = 0.0;

   // Iterate through all velocity blocks in this spatial cell
   // and calculate velocity moments:
   for (vmesh::LocalID block_i=0; block_i<SC->get_number_of_velocity_blocks(); ++block_i) {
      cpu_calcVelocityFirstMoments(SC,
         block_i,
         CellParams::RHO,
         CellParams::RHOVX,
         CellParams::RHOVY,
         CellParams::RHOVZ
      );
   }

   // Second iteration needed as rho has to be already computed when computing pressure
   for (vmesh::LocalID block_i=0; block_i<SC->get_number_of_velocity_blocks(); ++block_i) {
      cpu_calcVelocitySecondMoments(
         SC,
         block_i,
         CellParams::RHO,
         CellParams::RHOVX,
         CellParams::RHOVY,
         CellParams::RHOVZ,
         CellParams::P_11,
         CellParams::P_22,
         CellParams::P_33
      );
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
      calculateCellVelocityMoments(SC);
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

   }
   phiprof::stop("Calculate moments"); 

}
