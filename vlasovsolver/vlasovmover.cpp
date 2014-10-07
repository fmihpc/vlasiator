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
   typedef Parameters P;
   int trans_timer;
   
   
   phiprof::start("semilag-trans");
   phiprof::start("compute_cell_lists");
   const vector<CellID> local_cells = mpiGrid.get_cells();
   const vector<CellID> remote_stencil_cells_x = 
      mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
   const vector<CellID> remote_stencil_cells_y = 
      mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
   const vector<CellID> remote_stencil_cells_z = 
      mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
   phiprof::stop("compute_cell_lists");

   // ------------- SLICE - map dist function in Z --------------- //
   if(P::zcells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-z","MPI");
      phiprof::start(trans_timer);
      //start by doing all transfers in a blocking fashion (communication stage can be optimized separately) //
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA_TO_FLUXES);
      mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      phiprof::stop(trans_timer);
#pragma omp parallel
      {
         phiprof::start("prepare-block-data-z");
         for (size_t c=0; c<local_cells.size(); ++c)
            trans_prepare_block_data(mpiGrid,local_cells[c]);
         for (size_t c=0; c<remote_stencil_cells_z.size(); ++c)
            trans_prepare_block_data(mpiGrid,remote_stencil_cells_z[c]);
#pragma omp barrier
         phiprof::stop("prepare-block-data-z");
         phiprof::start("compute-mapping-z");
         for (size_t c=0; c<local_cells.size(); ++c) {
            if(do_translate_cell(mpiGrid[local_cells[c]]))
               trans_map_1d(mpiGrid,local_cells[c], 2, dt); // map along z//
         }
         phiprof::stop("compute-mapping-z");
      }
      
      trans_timer=phiprof::initializeTimer("update_remote-z","MPI");
      phiprof::start("update_remote-z");
      update_remote_mapping_contribution(mpiGrid, 2, 1);
      update_remote_mapping_contribution(mpiGrid, 2, -1);
      phiprof::stop("update_remote-z");
   }   
   
// ------------- SLICE - map dist function in X --------------- //
   if(P::xcells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-x","MPI");
      phiprof::start(trans_timer);
      //start by doing all transfers in a blocking fashion (communication stage can be optimized separately) //
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA_TO_FLUXES);
      mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);  
      phiprof::stop(trans_timer);
      
      #pragma omp parallel
      {
         phiprof::start("prepare-block-data-x");
         for (size_t c=0; c<local_cells.size(); ++c)
            trans_prepare_block_data(mpiGrid,local_cells[c]);
         for (size_t c=0; c<remote_stencil_cells_x.size(); ++c)
            trans_prepare_block_data(mpiGrid,remote_stencil_cells_x[c]);
         #pragma omp barrier
         phiprof::stop("prepare-block-data-x");

	 phiprof::start("compute-mapping-x");
         for (size_t c=0; c<local_cells.size(); ++c) {
            if (do_translate_cell(mpiGrid[local_cells[c]]))
               trans_map_1d(mpiGrid,local_cells[c], 0, dt); // map along x//
         }
         phiprof::stop("compute-mapping-x");
      }

      trans_timer=phiprof::initializeTimer("update_remote-x","MPI");
      phiprof::start("update_remote-x");
      update_remote_mapping_contribution(mpiGrid, 0, 1);
      update_remote_mapping_contribution(mpiGrid, 0, -1);
      phiprof::stop("update_remote-x");
   }
   
// ------------- SLICE - map dist function in Y --------------- //
   if(P::ycells_ini > 1 ){
      trans_timer=phiprof::initializeTimer("transfer-stencil-data-y","MPI");
      phiprof::start(trans_timer);
      //start by doing all transfers in a blocking fashion (communication stage can be optimized separately) //
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA_TO_FLUXES);
      mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);  
      phiprof::stop(trans_timer);
#pragma omp parallel
      {
         phiprof::start("prepare-block-data-y");
         for (size_t c=0; c<local_cells.size(); ++c)
            trans_prepare_block_data(mpiGrid,local_cells[c]);
         for (size_t c=0; c<remote_stencil_cells_y.size(); ++c)
            trans_prepare_block_data(mpiGrid,remote_stencil_cells_y[c]);
#pragma omp barrier
         phiprof::stop("prepare-block-data-y");         
         phiprof::start("compute-mapping-y");
         for (size_t c=0; c<local_cells.size(); ++c) {
            if(do_translate_cell(mpiGrid[local_cells[c]]))
               trans_map_1d(mpiGrid,local_cells[c], 1, dt); // map along y//
         }
         phiprof::stop("compute-mapping-y");
      }
      
      trans_timer=phiprof::initializeTimer("update_remote-y","MPI");
      phiprof::start("update_remote-y");
      update_remote_mapping_contribution(mpiGrid, 1, 1);
      update_remote_mapping_contribution(mpiGrid, 1, -1);
      phiprof::stop("update_remote-y");

   }
   
   // Mapping complete, update moments //
   phiprof::start("compute-moments-n-maxdt");
   // Note: Parallelization over blocks is not thread-safe
#pragma omp  parallel for
   for (size_t c=0; c<local_cells.size(); ++c) {
      SpatialCell* SC=mpiGrid[local_cells[c]];
      
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


void calculateAcceleration(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   Real dt
) {
   typedef Parameters P;
   const vector<CellID> cells = mpiGrid.get_cells();
   vector<CellID> propagatedCells;
   // Iterate through all local cells and propagate distribution functions 
   // in velocity space. Ghost cells (spatial cells at the boundary of the simulation 
   // volume) do not need to be propagated:
   
   
   //set initial cells to propagate
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* SC = mpiGrid[cells[c]];
      //disregard boundary cells
      //do not integrate cells with no blocks  (well, do not computes in practice)
      if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY &&
         SC->get_number_of_velocity_blocks() != 0) {
         propagatedCells.push_back(cells[c]);
      }
   }

   //Semilagrangian acceleration
   phiprof::start("semilag-acc");
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t c=0; c<propagatedCells.size(); ++c) {
      const CellID cellID = propagatedCells[c];
      phiprof::start("cell-semilag-acc");
      cpu_accelerate_cell(mpiGrid[cellID],dt);
      phiprof::stop("cell-semilag-acc");
   }
   phiprof::stop("semilag-acc");   
   
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
				       CellParams::P_33_V
				      );   //set second moments after acceleration
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
   vector<CellID> cells;
   cells=mpiGrid.get_cells();
   
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




void calculateCellVelocityMoments(SpatialCell* SC,
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
   vector<CellID> cells;
   cells=mpiGrid.get_cells();
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
