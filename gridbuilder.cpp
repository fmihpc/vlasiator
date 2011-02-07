#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <map>
#include <mpi.h>

#include "definitions.h"
#include "mpiconversion.h"
#include "common.h"
#include "parameters.h"
#include "cell_spatial.h"
#include "project.h"
#include "gridbuilder.h"

using namespace std;

// ******************************
// ***** GRIDBUILDERFACTORY *****
// ******************************

GridBuilderCreator GridBuilderFactory::gbc = NULL;

GridBuilderFactory::GridBuilderFactory() { }

GridBuilderFactory::~GridBuilderFactory() { }

GridBuilder* GridBuilderFactory::createBuilder() {
   if (gbc == NULL) return NULL;
   return (*gbc)();
}

bool GridBuilderFactory::deleteBuilder(GridBuilder*& gb) {
   delete gb;
   gb = NULL;
   return true;
}

bool GridBuilderFactory::registerBuilder(GridBuilderCreator gbc) {
   //cerr << "registerBuilder called" << endl;
   if (gbc == NULL) return false;
   GridBuilderFactory::gbc = gbc;
   return true;
}

// **********************************
// ***** GRIDBUILDER BASE CLASS *****
// **********************************

GridBuilder::GridBuilder() { }

GridBuilder::~GridBuilder() { }




inline uint velblock(cuint& iv,cuint& jv,cuint& kv) {
   typedef Parameters P;
   return kv*P::vyblocks_ini*P::vxblocks_ini + jv*P::vxblocks_ini + iv;
}

uint cellIndex(cuint& i,cuint& j,cuint& k) {
   typedef Parameters P;
   return k*P::ycells_ini*P::xcells_ini + j*P::xcells_ini + i;
}

/** Set up a spatial cell.
 * @param cell The spatial cell which is to be initialized.
 * @param xmin x-coordinate of the lower left corner of the cell.
 * @param ymin y-coordinate of the lower left corner of the cell.
 * @param zmin z-coordinate of the lower left corner of the cell.
 * @param dx Size of the cell in x-direction.
 * @param dy Size of the cell in y-direction.
 * @param dz Size of the cell in z-direction.
 * @param isRemote If true, the given cell is a remote cell (resides on another process) 
 * and its initial state need not be calculated.
 * @return If true, the cell was initialized successfully. Otherwise an error has 
 * occurred and the simulation should be aborted.
 */
bool buildSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,
		      creal& zmin,creal& dx,creal& dy,creal& dz,
		     const bool& isRemote) {
   typedef Parameters P;

   cuint VELBLOCKS = P::vxblocks_ini*P::vyblocks_ini*P::vzblocks_ini;
   // Set up cell parameters:
   cell.cpu_cellParams[CellParams::XCRD] = xmin;
   cell.cpu_cellParams[CellParams::YCRD] = ymin;
   cell.cpu_cellParams[CellParams::ZCRD] = zmin;
   cell.cpu_cellParams[CellParams::DX  ] = dx;
   cell.cpu_cellParams[CellParams::DY  ] = dy;
   cell.cpu_cellParams[CellParams::DZ  ] = dz;
   calcCellParameters(cell.cpu_cellParams,0.0);
   cell.cpu_cellParams[CellParams::RHO  ] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVX] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVY] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVZ] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = (P::vxmax-P::vxmin)/P::vxblocks_ini; // Size of a block in vx-direction
   creal dvy_block = (P::vymax-P::vymin)/P::vyblocks_ini; //                    vy
   creal dvz_block = (P::vzmax-P::vzmin)/P::vzblocks_ini; //                    vz
   creal dvx_blockCell = dvx_block / WID;                 // Size of one cell in a block in vx-direction
   creal dvy_blockCell = dvy_block / WID;                 //                                vy
   creal dvz_blockCell = dvz_block / WID;                 //                                vz
   creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;  // Volume of one cell in a block
   Real* const blockParams = cell.cpu_blockParams;
   Real* const avgs = cell.cpu_avgs;
   uint* const nbrsVel = cell.cpu_nbrsVel;

   for (uint kv=0; kv<P::vzblocks_ini; ++kv) for (uint jv=0; jv<P::vyblocks_ini; ++jv) for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
      cuint velIndex = velblock(iv, jv, kv);
      
      creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
      creal vy_block = P::vymin + jv*dvy_block; // vy-
      creal vz_block = P::vzmin + kv*dvz_block; // vz-
      
      calcBlockParameters(blockParams + velIndex*SIZE_BLOCKPARAMS);
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VXCRD] = vx_block;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VYCRD] = vy_block;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VZCRD] = vz_block;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVX  ] = dvx_blockCell;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVY  ] = dvy_blockCell;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVZ  ] = dvz_blockCell;

      if (isRemote == true) continue;
      // Calculate volume average of distrib. function for each cell in the block.
      for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	 creal vx_cell = vx_block + ic*dvx_blockCell;
	 creal vy_cell = vy_block + jc*dvy_blockCell;
	 creal vz_cell = vz_block + kc*dvz_blockCell;
	 avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic] =
	   calcPhaseSpaceDensity(xmin,ymin,zmin,dx,dy,dz,vx_cell,vy_cell,vz_cell,dvx_blockCell,dvy_blockCell,dvz_blockCell);
	 
	 // Add contributions to spatial cell velocity moments:
	 cell.cpu_cellParams[CellParams::RHO  ] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*dV;
	 cell.cpu_cellParams[CellParams::RHOVX] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vx_block + (ic+convert<Real>(0.5))*dvx_blockCell)*dV;
	 cell.cpu_cellParams[CellParams::RHOVY] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vy_block + (jc+convert<Real>(0.5))*dvy_blockCell)*dV;
	 cell.cpu_cellParams[CellParams::RHOVZ] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vz_block + (kc+convert<Real>(0.5))*dvz_blockCell)*dV;
      }
      
      // Since the whole velocity space is inside the spatial cell and the initial 
      // velocity grid contains just unrefined blocks,
      // the velocity neighbour lists can be constructed now:
      uint vxneg_nbr = iv-1;
      uint vxpos_nbr = iv+1;
      uint vyneg_nbr = jv-1;
      uint vypos_nbr = jv+1;
      uint vzneg_nbr = kv-1;
      uint vzpos_nbr = kv+1;
      uint state = 0;
      if (iv == 0) {              // -vx boundary
	 state = state | NbrsVel::VX_NEG_BND;
	 vxneg_nbr = P::vxblocks_ini-1;
      }
      if (iv == P::vxblocks_ini-1) { // +vx boundary
	 state = state | NbrsVel::VX_POS_BND;
	 vxpos_nbr = 0;
      }
      if (jv == 0) {              // -vy boundary
	 state = state | NbrsVel::VY_NEG_BND;
	 vyneg_nbr = P::vyblocks_ini-1;
      }
      if (jv == P::vyblocks_ini-1) { // +vy boundary
	 state = state | NbrsVel::VY_POS_BND;
	 vypos_nbr = 0;
      }
      if (kv == 0) {              // -vz boundary
	 state = state | NbrsVel::VZ_NEG_BND;
	 vzneg_nbr = P::vzblocks_ini-1;
      }
      if (kv == P::vzblocks_ini-1) { // +vz boundary
	 state = state | NbrsVel::VZ_POS_BND;
	 vzpos_nbr = 0;
      }
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] = state;
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::MYIND] = velIndex;
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VXNEG] = velblock(vxneg_nbr,jv,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VXPOS] = velblock(vxpos_nbr,jv,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VYNEG] = velblock(iv,vyneg_nbr,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VYPOS] = velblock(iv,vypos_nbr,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VZNEG] = velblock(iv,jv,vzneg_nbr);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VZPOS] = velblock(iv,jv,vzpos_nbr);
   }
   return true;
}

bool buildGrid(ParGrid<SpatialCell>& mpiGrid) {
   bool globalSuccess = 0;
   bool success = true;
   const int MASTER_RANK = 0;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   int N_processes;
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   
   // Create a GridBuilder:
   GridBuilder* builder = GridBuilderFactory::createBuilder();
   if (builder == NULL) {
      mpilogger << "(BUILDGRID) ERROR: Received NULL GridBuilder!" << endl << write;
      success = false;
   }
   if (success == true && builder->initialize() == false) {
      mpilogger << "(BUILDGRID) ERROR: GridBuilder failed to initialize!" << endl << write;
      success = false;
   }

   // Check that everything is ok:
   if (MPI_Allreduce(&success,&globalSuccess,1,MPI_Type<uchar>(),MPI_MIN,MPI_COMM_WORLD) != MPI_SUCCESS) {
      cerr << "(BUILDGRID) ERROR: Allreduce failed!" << std::endl; exit(1);
   }
   if (globalSuccess == false) {
      GridBuilderFactory::deleteBuilder(builder);
      return false;
   }
   
   // Get total number of cells to create:
   lluint N_cells;
   if (builder->getTotalNumberOfCells(N_cells) == false) {
      mpilogger << "(BUILDGRID) ERROR: Failed to get number of cells to create from GridBuilder!" << endl << write;
      success = false;
   }
   
   // ********************************************
   // ***** START TO PASS CELLS TO PROCESSES *****
   // ********************************************
   
   cuint MAX_NBRS = 49; // Cells can have max. 48 spatial neighbours
   lluint cellID;
   lluint nbrs[MAX_NBRS];
   uchar nbrTypes[MAX_NBRS];
   Real coords[3];
   Real sizes[3];
   MPI_Status mpiStatus;
   for (uint i=0; i<49; ++i) nbrs[i] = numeric_limits<lluint>::max();
   if (myrank == MASTER_RANK) {
      lluint cellsCreated = 0;
      for (int proc=0; proc<N_processes; ++proc) {
	 // Count the number of cells process proc receives:
	 lluint cellsToProc = N_cells/N_processes;
	 if (proc < N_cells%N_processes) ++cellsToProc;
	 
	 lluint cellCounter = 0;
	 while (cellCounter < cellsToProc) {
	    for (uint i=0; i<MAX_NBRS; ++i) nbrs[i] = numeric_limits<lluint>::max();
	    if (builder->getNextCell(48,cellID,coords,sizes,nbrs,nbrTypes) == false) {
	       // Check that we received the correct number of cells:
	       if (cellsCreated != N_cells) {
		  mpilogger << "(BUILDGRID) ERROR: GridBuilder gave too few cells!" << endl << write; 
		  success = false; 
	       }
	       break; // while
	    }
	    
	    if (proc == myrank) {
	       // Cells are for me. Add them to grid:
	       if (mpiGrid.addCell(cellID,coords,sizes,nbrs,nbrTypes) == false) {
		  mpilogger << "(BUILDGRID) ERROR: Failed to add cell to grid!" << endl << write;
		  success = false; 
		  break; // while
	       }
	    } else {
	       // Send cellID, nbrs, nbrTypes, coords, sizes to process proc (this is 
	       // definitely not the most efficient way to send data, but suffices for now):
	       if (MPI_Send(&cellID,1,MPI_Type<lluint>(),proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
		  mpilogger << "(BUILDGRID) ERROR occurred while sending cell data!" << endl << write; 
		  success = false;
	       }
	       if (MPI_Send(coords,3,MPI_Type<Real>(),proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
		  mpilogger << "(BUILDGRID) ERROR occurred while sending cell data!" << endl << write; 
		  success = false;
	       }
	       if (MPI_Send(sizes,3,MPI_Type<Real>(),proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
		  mpilogger << "(BUILDGRID) ERROR occurred while sending cell data!" << endl << write; 
		  success = false;
	       }
	       if (MPI_Send(nbrs,MAX_NBRS,MPI_Type<lluint>(),proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
		  mpilogger << "(BUILDGRID) ERROR occurred while sending cell data!" << endl << write; 
		  success = false;
	       }
	       if (MPI_Send(nbrTypes,MAX_NBRS,MPI_Type<uchar>(),proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
		  mpilogger << "(BUILDGRID) ERROR occurred while sending cell data!" << endl << write; 
		  success = false;
	       }
	    }
	    ++cellCounter;
	    ++cellsCreated;
	 }
	 if (success == false) break; // for
      }
      // Tell everyone that we are done:
      cellID = numeric_limits<lluint>::max();
      for (int proc=0; proc<N_processes; ++proc) {
	 if (MPI_Send(&cellID,1,MPI_Type<lluint>(),proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
	    cerr << "(BUILDGRID) ERROR occurred while sending exit signal!" << endl;
	 }
      }
   } else {
      do {
	 // Keep receiving cell IDs until master sends an invalid one:
	 if (MPI_Recv(&cellID,1,MPI_Type<lluint>(),MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	    mpilogger << "(BUILDGRID) ERROR occurred while receiving cells!" << endl << write; success = false;
	 }
	 if (cellID == numeric_limits<lluint>::max()) break;
	 // Received a valid cellID. Now receive the coordinates and connectivity data:
	 if (MPI_Recv(coords,3,MPI_Type<Real>(),MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	    mpilogger << "(BUILDGRID) ERROR occurred while receiving cells!" << endl << write; success = false;
	 }
	 if (MPI_Recv(sizes,3,MPI_Type<Real>(),MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	    mpilogger << "(BUILDGRID) ERROR occurred while receiving cells!" << endl << write; success = false;
	 }
	 if (MPI_Recv(nbrs,MAX_NBRS,MPI_Type<lluint>(),MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	    mpilogger << "(BUILDGRID) ERROR occurred while receiving cells!" << endl << write; success = false;
	 }
	 if (MPI_Recv(nbrTypes,MAX_NBRS,MPI_Type<uchar>(),MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	    mpilogger << "(BUILDGRID) ERROR occurred while receiving cells!" << endl << write; success = false;
	 }
	 // Add cell to grid:
	 if (mpiGrid.addCell(cellID,coords,sizes,nbrs,nbrTypes) == false) {
	    mpilogger << "(BUILDGRID) ERROR: Failed to add cell to grid!" << endl << write;
	    success = false;
	 }
      } while (true);
   }
   
   // Check that everything is ok:
   if (MPI_Allreduce(&success,&globalSuccess,1,MPI_Type<uchar>(),MPI_MIN,MPI_COMM_WORLD) != MPI_SUCCESS) {
      cerr << "(BUILDGRID) ERROR: Allreduce failed!" << std::endl; exit(1);
   }
   if (globalSuccess == false) return false;

   // At this point we know the cells and their connectivity. 
   // We can do an initial load balance.

   
   // **************************************************
   // ***** START TO PASS CELL PHYSICAL GEOMETRIES *****
   // **************************************************
   
   if (myrank == MASTER_RANK) {
      
   } else {
      
   }

   builder->finalize();
   GridBuilderFactory::deleteBuilder(builder);
   return success;
}

