/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <cmath>
#include <limits>
#include <map>
#include <mpi.h>
#include <iostream>

#include "definitions.h"
#include "mpiconversion.h"
#include "common.h"
#include "parameters.h"
#include "spatial_cell.hpp"
#include "project.h"
#include "gridbuilder.h"

using namespace std;
namespace VC = VirtualCell;

// ******************************
// ***** GRIDBUILDERFACTORY *****
// ******************************

// Init static member functions:
GridBuilderCreator GridBuilderFactory::gbc = NULL;

/** Constructor for GridBuilderFactory. The constructor is empty.
 */
GridBuilderFactory::GridBuilderFactory() { }

/** Destructor for GridBuilderFactory. The destructor is empty.
 */
GridBuilderFactory::~GridBuilderFactory() { }

/** Request a GridBuilder from GridBuilderFactory. The returned 
 * GridBuilder must be deallocated via GridBuilderFactory::deleteBuilder.
 * @return A pointer to a new GridBuilder. If a GridBuilder has not 
 * been registered, a NULL pointer is returned.
 * @see GridBuilderFactory::deleteBuilder
 * @see GridBuilderFactory::registerBuilder
 */
GridBuilder* GridBuilderFactory::createBuilder() {
   if (gbc == NULL) return NULL;
   return (*gbc)();
}

/** Deallocate a GridBuilder.
 * @param gb Pointer to GridBuilder which has been created with a call to GridBuilderFactory::createBuilder.
 * @return If true, the GridBuilder was deallocated successfully.
 * @see GridBuilderFactory::createBuilder
 */
bool GridBuilderFactory::deleteBuilder(GridBuilder*& gb) {
   delete gb;
   gb = NULL;
   return true;
}

/** Register a GridBuilder to GridBuilderFactory. If this is called multiple 
 * times, newer GridBuilderCreator functions overwrites the old one.
 * @param gbc A function which creates a new GridBuilder.
 * @return If true, the GridBuilder was registered successfully.
 */
bool GridBuilderFactory::registerBuilder(GridBuilderCreator gbc) {
   if (gbc == NULL) return false;
   GridBuilderFactory::gbc = gbc;
   return true;
}

// **********************************
// ***** GRIDBUILDER BASE CLASS *****
// **********************************

/** GridBuilder base class constructor. The constructor is empty.
 */
GridBuilder::GridBuilder() { }

/** GridBuilder base class virtual destructor. The destructor is empty.
 */
GridBuilder::~GridBuilder() { }

bool GridBuilder::doInitialLoadBalance() {return true;}


inline uint velblock(cuint& iv,cuint& jv,cuint& kv) {
   typedef Parameters P;
   return kv*P::vyblocks_ini*P::vxblocks_ini + jv*P::vxblocks_ini + iv;
}

uint cellIndex(cuint& i,cuint& j,cuint& k) {
   typedef Parameters P;
   return k*P::ycells_ini*P::xcells_ini + j*P::xcells_ini + i;
}

/** Set up a spatial cell.
 * DEPRECATED: TO BE REMOVED IN THE FUTURE!
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

   // Set up cell parameters:
   cell.parameters[CellParams::XCRD] = xmin;
   cell.parameters[CellParams::YCRD] = ymin;
   cell.parameters[CellParams::ZCRD] = zmin;
   cell.parameters[CellParams::DX  ] = dx;
   cell.parameters[CellParams::DY  ] = dy;
   cell.parameters[CellParams::DZ  ] = dz;
   calcCellParameters(cell.parameters,0.0);
   cell.parameters[CellParams::RHO  ] = 0.0;
   cell.parameters[CellParams::RHOVX] = 0.0;
   cell.parameters[CellParams::RHOVY] = 0.0;
   cell.parameters[CellParams::RHOVZ] = 0.0;
   
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
	 cell.parameters[CellParams::RHO  ] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*dV;
	 cell.parameters[CellParams::RHOVX] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vx_block + (ic+convert<Real>(0.5))*dvx_blockCell)*dV;
	 cell.parameters[CellParams::RHOVY] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vy_block + (jc+convert<Real>(0.5))*dvy_blockCell)*dV;
	 cell.parameters[CellParams::RHOVZ] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vz_block + (kc+convert<Real>(0.5))*dvz_blockCell)*dV;
      }
      // Create velocity neighbour list entry:
      uint nbrFlags = 0;
      for (uint kkv=0; kkv<3; ++kkv) for (uint jjv=0; jjv<3; ++jjv) for (uint iiv=0; iiv<3; ++iiv) {
	 // By default set the block as non-existing:
	 nbrsVel[velIndex*SIZE_NBRS_VEL + kkv*9+jjv*3+iiv] = numeric_limits<uint>::max();
	 
	 // Then check if the block actually exists:
	 cuint iindex = iv + iiv - 1;
	 cuint jindex = jv + jjv - 1;
	 cuint kindex = kv + kkv - 1;
	 if (iindex >= P::vxblocks_ini) continue;
	 if (jindex >= P::vyblocks_ini) continue;
	 if (kindex >= P::vzblocks_ini) continue;
	 
	 // Block exists, override values set above:
	 nbrsVel[velIndex*SIZE_NBRS_VEL + kkv*9+jjv*3+iiv] = velblock(iindex,jindex,kindex);
	 nbrFlags = (nbrFlags | (1 << (kkv*9+jjv*3+iiv)));
      }
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::NBRFLAGS] = nbrFlags;	 
   }
   creal spatialVolume = cell.parameters[CellParams::DX]*cell.parameters[CellParams::DY]*cell.parameters[CellParams::DZ];
   cell.parameters[CellParams::RHO  ] /= spatialVolume;
   cell.parameters[CellParams::RHOVX] /= spatialVolume;
   cell.parameters[CellParams::RHOVY] /= spatialVolume;
   cell.parameters[CellParams::RHOVZ] /= spatialVolume;
   return true;
}

bool broadcastMandatoryParameters(const int& myrank,const int& MASTER_RANK,MPI_Comm comm,GridBuilder* const builder) {
   bool success = true;
   Real q,m,dt,t_min;
   luint timestep,max_timesteps;
   if (myrank == MASTER_RANK) {
      q             = builder->getSpeciesCharge();
      m             = builder->getSpeciesMass();
      dt            = builder->getDt();
      t_min         = builder->getTmin();
      timestep      = builder->getTimestep();
      max_timesteps = builder->getMaxTimesteps();
   }
   if (MPI_Bcast(&q            ,1,MPI_Type<Real>(),MASTER_RANK,comm)  != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&m            ,1,MPI_Type<Real>(),MASTER_RANK,comm)  != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&dt           ,1,MPI_Type<Real>(),MASTER_RANK,comm)  != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&t_min        ,1,MPI_Type<Real>(),MASTER_RANK,comm)  != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&timestep     ,1,MPI_Type<luint>(),MASTER_RANK,comm) != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&max_timesteps,1,MPI_Type<luint>(),MASTER_RANK,comm) != MPI_SUCCESS) success = false;
   Parameters::q = q;
   Parameters::m = m;
   Parameters::q_per_m = q/m;
   Parameters::dt = dt;
   Parameters::t = t_min + timestep*dt;
   Parameters::tstep = timestep;
   Parameters::tstep_min = timestep;
   Parameters::tsteps = max_timesteps;
   return success;
}

/** Build the simulation grid. This is done so that master process first requests new cells 
 * from a user-defined GridBuilder, and then distributes them to other processes.
 * @param mpiGrid The parallel grid.
 * @return If true, the grid was created and distributed successfully. It should be 
 * load balanced.
 */ 
#include "mpilogger.h"
extern MPILogger mpilogger;

bool buildGrid(MPI_Comm comm,const int& MASTER_RANK) {
   bool success = true;
   int myrank,N_processes;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   
   GridBuilder* builder = GridBuilderFactory::createBuilder();
   if (builder == NULL) {
      mpilogger << "(BUILDGRID) ERROR: Received NULL GridBuilder!" << endl << write;
      success = false;
   }
   if (success == true && builder->initialize(comm,MASTER_RANK) == false) {
      mpilogger << "(BUILDGRID) ERROR: GridBuilder failed to initialize!" << endl << write;
      success = false;
   }
   if (broadcastMandatoryParameters(myrank,MASTER_RANK,comm,builder) == false) success = false;
   
   
   // Deallocate memory and exit:
   builder->finalize();
   GridBuilderFactory::deleteBuilder(builder);
   return success;
}
