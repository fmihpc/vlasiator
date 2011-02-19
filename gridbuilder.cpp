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
   creal spatialVolume = cell.cpu_cellParams[CellParams::DX]*cell.cpu_cellParams[CellParams::DY]*cell.cpu_cellParams[CellParams::DZ];
   cell.cpu_cellParams[CellParams::RHO  ] /= spatialVolume;
   cell.cpu_cellParams[CellParams::RHOVX] /= spatialVolume;
   cell.cpu_cellParams[CellParams::RHOVY] /= spatialVolume;
   cell.cpu_cellParams[CellParams::RHOVZ] /= spatialVolume;
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
#ifndef PARGRID
#include "mpilogger.h"
extern MPILogger mpilogger;

bool buildGrid(MPI_Comm comm,const int& MASTER_RANK) {
   bool success = true;
   GridBuilder* builder = GridBuilderFactory::createBuilder();
   if (builder == NULL) {
      mpilogger << "(BUILDGRID) ERROR: Received NULL GridBuilder!" << endl << write;
      success = false;
   }
   if (success == true && builder->initialize(comm,MASTER_RANK) == false) {
      mpilogger << "(BUILDGRID) ERROR: GridBuilder failed to initialize!" << endl << write;
      success = false;
   }
   builder->finalize();
   GridBuilderFactory::deleteBuilder(builder);
   return success;
}
#else
bool buildGrid(ParGrid<SpatialCell>& mpiGrid,MPI_Comm comm,const int& MASTER_RANK) {
   bool globalSuccess = 0;
   bool success = true;
   int myrank;
   MPI_Comm_rank(comm,&myrank);
   int N_processes;
   MPI_Comm_size(comm,&N_processes);
   VC::ID INVALID_CELL_ID = numeric_limits<VC::ID>::max();
   
   // Create a GridBuilder:
   GridBuilder* builder = GridBuilderFactory::createBuilder();
   if (builder == NULL) {
      mpilogger << "(BUILDGRID) ERROR: Received NULL GridBuilder!" << endl << write;
      success = false;
   }
   if (success == true && builder->initialize(comm,MASTER_RANK) == false) {
      mpilogger << "(BUILDGRID) ERROR: GridBuilder failed to initialize!" << endl << write;
      success = false;
   }
   if (success == false) cerr << "builder init failed on proc #" << myrank << endl;
   
   // Check that everything is ok:
   if (MPI_Allreduce(&success,&globalSuccess,1,MPI_Type<uchar>(),MPI_MIN,comm) != MPI_SUCCESS) {
      cerr << "(BUILDGRID) ERROR: Allreduce failed!" << endl; exit(1);
   }
   if (globalSuccess == false) {
      builder->finalize();
      GridBuilderFactory::deleteBuilder(builder);
      return false;
   }
   
   // Get values of mandatory parameters and broadcast them to all processes
   if (broadcastMandatoryParameters(myrank,MASTER_RANK,comm,builder) == false) success = false;
   
   // Get total number of cells to create:
   VC::ID N_cells;
   if (myrank == MASTER_RANK) if (builder->getTotalNumberOfCells(N_cells) == false) {
      mpilogger << "(BUILDGRID) ERROR: Failed to get number of cells to create from GridBuilder!" << endl << write;
      success = false;
   }
   
   // ********************************************
   // ***** START TO PASS CELLS TO PROCESSES *****
   // ********************************************
   vector<VC::ID> cellIDs;
   vector<uchar>  spatNbrsPerCell;
   // Request cells from GridBuilder:
   if (myrank == MASTER_RANK) if (builder->getCellIDs(cellIDs,spatNbrsPerCell) == false) success = false;
   // Check that master received cell list successfully:
   if (MPI_Allreduce(&success,&globalSuccess,1,MPI_Type<uchar>(),MPI_MIN,comm) != MPI_SUCCESS) {
      cerr << "(BUILDGRID) ERROR: Allreduce failed!" << endl; exit(1);
   }
   if (globalSuccess == false) {
      builder->finalize();
      GridBuilderFactory::deleteBuilder(builder);
      return false;
   }
   
   // Count the number of cells allocated to each process and scatter to all processes:
   VC::ID* cellsPerProcess;
   if (myrank == MASTER_RANK) {
      cellsPerProcess = new VC::ID[N_processes];
      for (int proc=0; proc<N_processes; ++proc) {
	 cellsPerProcess[proc] = N_cells/N_processes;
	 if (proc < N_cells%N_processes) ++cellsPerProcess[proc];
      }
   }
   VC::ID N_myCells = 0;
   if (MPI_Scatter(cellsPerProcess,1,MPI_Type<VC::ID>(),&N_myCells,1,MPI_Type<VC::ID>(),MASTER_RANK,comm) != MPI_SUCCESS) success = false;
   
   // Send allocated cells, and number of spatial neighbours per each cell, to all processes:
   VC::ID* myCellIDs  = new VC::ID[N_myCells];
   uchar*  myCellNbrs = new uchar[N_myCells];
   MPI_Request* mpiRequests = NULL;
   if (myrank == MASTER_RANK) {
      VC::ID counter = 0;
      mpiRequests = new MPI_Request[2*N_processes];
      for (uint i=0; i<2*N_processes; ++i) mpiRequests[i] = MPI_REQUEST_NULL;
      for (int proc=0; proc<N_processes; ++proc) {
	 if (myrank != proc) {
	    MPI_Isend(&(cellIDs[counter])        ,cellsPerProcess[proc],MPI_Type<VC::ID>(),proc,proc,comm,&(mpiRequests[2*proc+0]));
	    MPI_Isend(&(spatNbrsPerCell[counter]),cellsPerProcess[proc],MPI_Type<uchar>() ,proc,proc,comm,&(mpiRequests[2*proc+1]));
	 }
	 counter += cellsPerProcess[proc];
      }
      // Master copies its own cells & neighbours manually (assume MASTER_RANK=0):
      for (VC::ID i=0; i<N_myCells; ++i) {
	 myCellIDs[i]  = cellIDs[i];
	 myCellNbrs[i] = spatNbrsPerCell[i];
      }
      if (MPI_Waitall(2*N_processes,mpiRequests,MPI_STATUSES_IGNORE) != MPI_SUCCESS) success = false;
      delete cellsPerProcess;
      cellIDs.resize(0);
      spatNbrsPerCell.resize(0);
   } else {
      // Slave processes post receives and wait:
      mpiRequests = new MPI_Request[2];
      MPI_Irecv(myCellIDs ,N_myCells,MPI_Type<VC::ID>(),MASTER_RANK,myrank,comm,&(mpiRequests[0]));
      MPI_Irecv(myCellNbrs,N_myCells,MPI_Type<uchar>() ,MASTER_RANK,myrank,comm,&(mpiRequests[1]));
      if (MPI_Waitall(2,mpiRequests,MPI_STATUSES_IGNORE) != MPI_SUCCESS) success = false;
   }
   delete mpiRequests;
   
   // Now every process knows how many cells it will initially have, as well as 
   // the number of spatial neighbours and velocity blocks for every cell. Very likely 
   // the obtained partitioning of cells is a poor one, thus we should do an initial 
   // load balance with minimum amount of (cell) data.
   // 
   // Before a load balance we need to get the cell coordinates and neighbour data 
   // from GridBuilder.
   // 
   // Cell coordinates are only requested at this point because of Zoltan's coordinate-based 
   // load balancing methods (RCB, RIB, HSFC, ...).
   
   // Count the sum of neighbours of each cell (needed for mem allocation):
   VC::ID neighbourCount = 0;
   for (VC::ID i=0; i<N_myCells; ++i) neighbourCount += myCellNbrs[i];
   VC::ID N_myNbrCount = neighbourCount;

   // Allocate buffers:
   Real* coordsBuffer    = new Real[6*N_myCells];
   VC::ID* nbrIDBuffer   = new VC::ID[N_myNbrCount];
   uchar* nbrTypesBuffer = new uchar[N_myNbrCount];

   // Send requests for neighbour data (all), process requests (master), 
   // and wait for data to arrive (all):
   neighbourCount = 0;

   if (builder->addCellNbrRequests(N_myCells,N_myNbrCount,myCellIDs,myCellNbrs,coordsBuffer,nbrIDBuffer,nbrTypesBuffer) == false) return false;
   if (builder->processCellNbrRequests() == false) success = false;
   if (builder->waitCellNbrRequests() == false) success = false;

   if (success == false) {
      cerr << "Failed to get cell nbr data!" << endl;
   }

   // Now we have the cell IDs, their global IDs of their neighbours and cell coordinates. 
   // We can add cells to parallel grid:
   neighbourCount = 0;
   for (VC::ID i=0; i<N_myCells; ++i) {
      vector<VC::ID> tmpNbrIDs;
      vector<uchar>  tmpNbrTypes;
      tmpNbrIDs.resize(myCellNbrs[i]);
      tmpNbrTypes.resize(myCellNbrs[i]);
      for (uchar j=0; j<myCellNbrs[i]; ++j) tmpNbrIDs[j]   = nbrIDBuffer[neighbourCount + j];
      for (uchar j=0; j<myCellNbrs[i]; ++j) tmpNbrTypes[j] = nbrTypesBuffer[neighbourCount + j];
      if (mpiGrid.addCell(myCellIDs[i],&(coordsBuffer[i*6]),&(coordsBuffer[i*6+3]),tmpNbrIDs,tmpNbrTypes) == false) success = false;
      neighbourCount += myCellNbrs[i];
   }
   if (success == true) if (mpiGrid.initialize() == false) success = false;
   
   // Check that everything is ok:
   if (MPI_Allreduce(&success,&globalSuccess,1,MPI_Type<uchar>(),MPI_MIN,comm) != MPI_SUCCESS) {
      cerr << "(BUILDGRID) ERROR: Allreduce failed!" << std::endl; exit(1);
   }
   if (globalSuccess == false) success = false;
   delete myCellIDs;
   delete myCellNbrs;
   delete coordsBuffer;
   delete nbrIDBuffer;
   delete nbrTypesBuffer;
   nbrTypesBuffer = NULL;
   nbrIDBuffer = NULL;
   coordsBuffer = NULL;
   myCellIDs = NULL;
   myCellNbrs = NULL;
   
   // Parallel grid has possibly done a load balance, so we need to request them 
   // from the grid (instead of using the ones we got from GridBuilder).
   vector<ID::type> localCells;
   mpiGrid.getCells(localCells);
   myCellIDs = new VC::ID[localCells.size()];
   N_myCells = localCells.size();
   for (size_t i=0; i<localCells.size(); ++i) {
      // Copy cell coordinates to SpatialCell
      SpatialCell* SC = mpiGrid[localCells[i]];
      SC->cpu_cellParams[CellParams::XCRD] = mpiGrid.get_cell_x(localCells[i]);
      SC->cpu_cellParams[CellParams::YCRD] = mpiGrid.get_cell_y(localCells[i]);
      SC->cpu_cellParams[CellParams::ZCRD] = mpiGrid.get_cell_z(localCells[i]);
      SC->cpu_cellParams[CellParams::DX  ] = mpiGrid.get_cell_x_size(localCells[i]);
      SC->cpu_cellParams[CellParams::DY  ] = mpiGrid.get_cell_y_size(localCells[i]);
      SC->cpu_cellParams[CellParams::DZ  ] = mpiGrid.get_cell_z_size(localCells[i]);
      
      // Request cell parameters for each spatial cell (need to copy cell IDs from 
      // mpiGrid to a temp array, because pargrid has uint IDs and we need VirtualCell::IDs
      // here):
      myCellIDs[i] = localCells[i];
   }
   Real* cellParamsBuffer = new Real[SIZE_CELLPARAMS*N_myCells];
   if (builder->addCellParamsRequests(N_myCells,myCellIDs,cellParamsBuffer) == false) success = false;
   if (builder->processCellParamsRequests() == false) success = false;
   if (builder->waitCellParamsRequests() == false) success = false;

   for (VC::ID i=0; i<N_myCells; ++i) {
      SpatialCell* SC = mpiGrid[localCells[i]];
      for (uint j=0; j<SIZE_CELLPARAMS; ++j) SC->cpu_cellParams[j] = cellParamsBuffer[i*SIZE_CELLPARAMS+j];
   }
   delete cellParamsBuffer;
   cellParamsBuffer = NULL;
   
   // Get the number of velocity blocks per cell so that we can allocate memory for 
   // blocks and block neighbour lists:
   uint* blocksPerCell = new uint[N_myCells];
   builder->addCellBlockNumberRequests(N_myCells,myCellIDs,blocksPerCell);
   builder->processCellBlockNumberRequests();
   builder->waitCellBlockNumberRequests();
   for (size_t i=0; i<localCells.size(); ++i) {
      mpiGrid[localCells[i]]->initialize(blocksPerCell[i]);
   }

   // Now we can finally load values of distribution function, block parameters 
   // and neighbour lists for each block:
   Real** avgsBuffer        = new Real*[N_myCells];
   Real** blockParamsBuffer = new Real*[N_myCells];
   uint** nbrsVelBuffer     = new uint*[N_myCells];
   for (size_t i=0; i<localCells.size(); ++i) {
      SpatialCell* SC      = mpiGrid[localCells[i]];
      avgsBuffer[i]        = SC->cpu_avgs;
      blockParamsBuffer[i] = SC->cpu_blockParams;
      nbrsVelBuffer[i]     = SC->cpu_nbrsVel;
      
      for (uint j=0; j<SIZE_VELBLOCK*blocksPerCell[i]; ++j) avgsBuffer[i][j] = 1.0;
   }
   if (builder->addCellBlockDataRequests(N_myCells,myCellIDs,blocksPerCell,avgsBuffer,blockParamsBuffer,nbrsVelBuffer) == false) success = false;
   if (builder->processCellBlockDataRequests() == false) success = false;
   if (builder->waitCellBlockDataRequests() == false) success = false;

   // Deallocate memory:
   for (VC::ID i=0; i<localCells.size(); ++i) {
      avgsBuffer[i] = NULL;
      blockParamsBuffer[i] = NULL;
      nbrsVelBuffer[i] = NULL;
   }
   delete avgsBuffer;
   delete blockParamsBuffer;
   delete nbrsVelBuffer;
   delete blocksPerCell;
   blocksPerCell = NULL;   
   delete myCellIDs;
   myCellIDs = NULL;
   
   builder->finalize();
   GridBuilderFactory::deleteBuilder(builder);   
   return success;
}
#endif
