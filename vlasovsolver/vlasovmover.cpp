/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute

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


/*--------------------------------------------------
Translation function (real space propagation)
--------------------------------------------------*/

/*!  

Propagates the distribution function in spatial space. 
 
Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
three‐dimensional monotone and conservative semi‐Lagrangian scheme
(SLICE‐3D) for transport problems." Quarterly Journal of the Royal
Meteorological Society 138.667 (2012): 1640-1651.

*/


void calculateSpatialTranslation(dccrg::Dccrg<SpatialCell>& mpiGrid,
				 creal dt) {
  typedef Parameters P;
  int trans_timer;
  

  phiprof::start("semilag-trans");
  phiprof::start("compute_cell_lists");
  const vector<CellID> local_cells = mpiGrid.get_cells();
  const vector<CellID> remote_translated_cells = 
    mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_SOURCE_NEIGHBORHOOD_ID);
  vector<CellID> propagated_cells; /*cells has all cells that we want to translate*/
  propagated_cells.reserve(local_cells.size() + remote_translated_cells.size() ); //reserve space
  /*add local cells */
  for(uint cell_i = 0;cell_i< local_cells.size();cell_i++) {
    const CellID cellID = local_cells[cell_i];
    SpatialCell* SC=mpiGrid[cellID];
    if(SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       (SC->sysBoundaryLayer != 1  &&
	SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)
       ) continue;
    propagated_cells.push_back(cellID);
  }
  /*add remote cells*/
  for(uint cell_i = 0;cell_i< remote_translated_cells.size();cell_i++) {
    const CellID cellID = remote_translated_cells[cell_i];
    SpatialCell* SC=mpiGrid[cellID];
    if(SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       (SC->sysBoundaryLayer != 1  &&
	SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)
       ) continue;
    propagated_cells.push_back(cellID);
  }
  phiprof::stop("compute_cell_lists");  


  /*propagate*/
  trans_timer=phiprof::initializeTimer("transfer-stencil-data-z","MPI");
  phiprof::start(trans_timer);
  /*start by doing all transfers in a blocking fashion (communication stage can be optimized separately) */
  SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
  mpiGrid.update_remote_neighbor_data(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);  
  phiprof::stop(trans_timer);
  phiprof::start("compute-mapping-z");
  for (size_t c=0; c<propagated_cells.size(); ++c) {
     map_trans_1d(mpiGrid,propagated_cells[c], 2, dt); /*< map along z*/
  }
  phiprof::stop("compute-mapping-z");

  trans_timer=phiprof::initializeTimer("transfer-stencil-data-x","MPI");
  phiprof::start(trans_timer);
  /*start by doing all transfers in a blocking fashion (communication stage can be optimized separately) */
  SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
  mpiGrid.update_remote_neighbor_data(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);  
  phiprof::stop(trans_timer);
  phiprof::start("compute-mapping-x");
  for (size_t c=0; c<propagated_cells.size(); ++c) {
     map_trans_1d(mpiGrid,propagated_cells[c], 0, dt); /*< map along x*/
  }
  phiprof::stop("compute-mapping-x");
    
  trans_timer=phiprof::initializeTimer("transfer-stencil-data-y","MPI");
  phiprof::start(trans_timer);
  /*start by doing all transfers in a blocking fashion (communication stage can be optimized separately) */
  SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
  mpiGrid.update_remote_neighbor_data(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);  
  phiprof::stop(trans_timer);
  phiprof::start("compute-mapping-y");
  for (size_t c=0; c<propagated_cells.size(); ++c) {
     map_trans_1d(mpiGrid,propagated_cells[c], 1, dt); /*< map along y*/
  }
  phiprof::stop("compute-mapping-y");
  

  phiprof::start("compute-moments-n-maxdt");
  //move data from the fx arrays where the new results are stored, to the actual data array
#pragma omp  parallel for
  for (size_t c=0; c<local_cells.size(); ++c) {
    SpatialCell* SC=mpiGrid[local_cells[c]];
    const Real dx=SC->parameters[CellParams::DX];
    const Real dy=SC->parameters[CellParams::DY];
    const Real dz=SC->parameters[CellParams::DZ];
    SC->parameters[CellParams::RHO_V  ] = 0.0;
    SC->parameters[CellParams::RHOVX_V] = 0.0;
    SC->parameters[CellParams::RHOVY_V] = 0.0;
    SC->parameters[CellParams::RHOVZ_V] = 0.0;

    //Reset spatial max DT
    SC->parameters[CellParams::MAXRDT]=numeric_limits<Real>::max();
    for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
      unsigned int block = SC->velocity_block_list[block_i];
      Velocity_Block* block_ptr = SC->at(block);
      const Real* const blockParams = block_ptr->parameters;
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

      //compute moments for this block
      cpu_calcVelocityMoments(SC,block,CellParams::RHO_R,CellParams::RHOVX_R,CellParams::RHOVY_R,CellParams::RHOVZ_R);   //set moments after translation
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
   dccrg::Dccrg<SpatialCell>& mpiGrid,
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
      if(SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY &&
         SC->number_of_blocks != 0) {
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

      for(unsigned int block_i=0; block_i< mpiGrid[cellID]->number_of_blocks;block_i++){
         unsigned int block = mpiGrid[cellID]->velocity_block_list[block_i];         
         cpu_calcVelocityMoments(mpiGrid[cellID],block,CellParams::RHO_V,CellParams::RHOVX_V,CellParams::RHOVY_V,CellParams::RHOVZ_V);   //set moments after acceleration
      }
   }
   phiprof::stop("Compute moments");
}




/*--------------------------------------------------
Functions for computing moments
--------------------------------------------------*/




void calculateInterpolatedVelocityMoments(dccrg::Dccrg<SpatialCell>& mpiGrid,
                                          const int cp_rho, const int cp_rhovx, const int cp_rhovy, const int cp_rhovz) {
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
      
      }
   }
}




void calculateCellVelocityMoments(
   SpatialCell* SC,
   bool doNotSkip
){
   if(!doNotSkip &&
         (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
            (SC->sysBoundaryLayer != 1  &&
             SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
   ) return;
   SC->parameters[CellParams::RHO  ] = 0.0;
   SC->parameters[CellParams::RHOVX] = 0.0;
   SC->parameters[CellParams::RHOVY] = 0.0;
   SC->parameters[CellParams::RHOVZ] = 0.0;
   //Iterate through all velocity blocks in this spatial cell 
   // and calculate velocity moments:
   for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
      unsigned int block = SC->velocity_block_list[block_i];
      cpu_calcVelocityMoments(SC,block,CellParams::RHO,CellParams::RHOVX,CellParams::RHOVY,CellParams::RHOVZ);
   }
}
 

void calculateVelocityMoments(dccrg::Dccrg<SpatialCell>& mpiGrid){
   vector<CellID> cells;
   cells=mpiGrid.get_cells();
   phiprof::start("Calculate moments"); 
   // Iterate through all local cells (incl. system boundary cells):
#pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      calculateCellVelocityMoments(SC);
   }
   phiprof::stop("Calculate moments"); 
}
 

