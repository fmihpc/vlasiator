// GENERAL NOTES
// 
// "Task queues" seem to be ~50% faster on MPI waiting times on meteo than the "simple" 
// version (-DSIMPLE) in which inner cells are calculated first, then arrival 
// of all neighbour data is waited, and finally boundary cells calculated.

#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/mpi.hpp>
#include <zoltan.h>

#include "../vlasovmover.h"
#include "../profile.h"
#include "../memalloc.h"
#include "cpu_acc_leveque.h"
#include "cpu_trans_leveque.h"
#include "../priorityqueue.h"
#include "../mpilogger.h"
using namespace std;

#include <stdint.h>
#define DCCRG_SEND_SINGLE_CELLS
#define DCCRG_CELL_DATA_SIZE_FROM_USER
#define DCCRG_USER_MPI_DATA_TYPE
#include <dccrg.hpp>
typedef uint64_t CellID;

#include "../transferstencil.h"
static TransferStencil<CellID> stencilAverages(INVALID_CELLID);
static TransferStencil<CellID> stencilUpdates(INVALID_CELLID);

static set<CellID> ghostCells;


static map<pair<CellID,int>,Real*> updateBuffers; /**< For each local cell receiving one or more remote df/dt updates,
						   * MPI rank of remote process sending an update and address to the 
						   * allocated buffer. */
static map<CellID,set<Real*> > remoteUpdates;     /**< For each local cell receiving one or more remote df/dt updates, 
						   * a set containing addresses of all allocated buffers. Note that map 
						   * remoteUpdates is only used to iterate over all df/dt buffers, which 
						   * is inconvenient to do with updateBuffers. updateBuffers is in convenient 
						   * form to post MPI receives, remoteUpdates is convenient to iterate 
						   * all cell's remote updates.*/


//??
namespace ID {
   typedef unsigned int type;
}

extern MPILogger mpilogger;

CellID getNeighbourID(
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid,
	#else
	const dccrg<SpatialCell>& mpiGrid,
	#endif
	const CellID& cellID,
	const uchar& i,
	const uchar& j,
	const uchar& k
) {
   #ifdef PARGRID
   const uchar nbrTypeID = calcNbrTypeID(i,j,k);
   return mpiGrid.getNeighbour(cellID,nbrTypeID);
   #else
   // TODO: merge this with the one in lond...anna.cpp
   const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of(cellID, int(i) - 2, int(j) - 2, int(k) - 2);
   if (neighbors.size() == 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << " No neighbor for cell " << cellID
         << " at offsets " << int(i) - 2 << ", " << int(j) - 2 << ", " << int(k) - 2
         << std::endl;
      abort();
   }
   // TODO support spatial refinement
   return neighbors[0];
   #endif
}

bool initializeMover(dccrg<SpatialCell>& mpiGrid) { 
   
   // Populate spatial neighbour list:
   vector<CellID> cells,remoteCells;
   std::vector<CellID> nbrs; //temporary vector for neighbors at certain offset
/* REPLACED  mpiGrid.getAllCells(cells); */
   
   cells=mpiGrid.get_cells();
   remoteCells=mpiGrid.get_list_of_remote_cells_with_local_neighbours();
   cells.insert( cells.end(), remoteCells.begin(), remoteCells.end() );
   
   for (size_t cell=0; cell<cells.size(); ++cell) {
       cuint cellID = cells[cell];
      uint* const nbrsSpa = mpiGrid[cellID]->cpu_nbrsSpa;
      bool isGhost = false;
      
      // Get spatial neighbour IDs and store them into a vector:
      uint counter = 0;
      vector<CellID> nbrIDs;
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
          // in dccrg cells are neighbors of themselves only with periodic boundaries
          if (i == 0 && j == 0 && k == 0) {
             nbrIDs.push_back(cellID);
          } else {
              //REPLACED nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+i,2+j,2+k)));
              nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2+i, 2+j, 2+k));
              if (nbrIDs.back() == INVALID_CELLID) {
                  isGhost = true;
              }
          }
          ++counter;
      }
      /*REPLACED
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(0,2,2))); // i-2,j,k nbr, goes to index 27
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,0,2))); // i,j-2,k nbr, goes to index 28
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,0))); // i,j,k-2 nbr, goes to index 29
      */
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 0, 2, 2));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2, 0, 2));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2, 2, 0));

      
      // Store neighbour offsets into a vector:
      vector<uint> cellOffsets(nbrIDs.size());
      for (size_t i=0; i<nbrIDs.size(); ++i) {
	 if (nbrIDs[i] == INVALID_CELLID) cellOffsets[i] = numeric_limits<uint>::max();
	 else cellOffsets[i] = mpiGrid[nbrIDs[i]]->cpuIndex;
      }
      
      // Create spatial neighbour list entry for each block. Offsets to missing 
      // neighbours are replaced by offsets to this cell so that cells on the 
      // boundary of the simulation domain (="ghost cells") can be propagated.
      // This allows boundary values (set by user) to propagate into simulation 
      // domain.
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 uint boundaryFlag = 0;
	 // Store offsets to each spatial neighbour of this block. Note that
	 // the offset to this block is stored to index 13:
	 for (size_t i=0; i<nbrIDs.size(); ++i) {
	    if (cellOffsets[i] == numeric_limits<uint>::max()) {
	       //nbrsSpa[block*SIZE_NBRS_SPA + i] = numeric_limits<CellID>::max();
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = cellOffsets[13] + block;
	    } else {
	       boundaryFlag = boundaryFlag | (1 << i);
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = cellOffsets[i] + block;
	    }
	 }
	 // Boundary flags are stored to the last position in nbrsSpa array:
	 nbrsSpa[block*SIZE_NBRS_SPA + 30] = boundaryFlag;
      }
      
      if (isGhost == true) {
	 ghostCells.insert(cellID);
      }
   }
   
   // ***** Calculate MPI send/receive stencils *****

   // Send/receive stencils for avgs:
   vector<Offset> nbrOffsets;   
   nbrOffsets.push_back(Offset(-1, 0, 0));
   nbrOffsets.push_back(Offset( 1, 0, 0));
   nbrOffsets.push_back(Offset( 0,-1, 0));
   nbrOffsets.push_back(Offset( 0, 1, 0));
   nbrOffsets.push_back(Offset( 0, 0,-1));
   nbrOffsets.push_back(Offset( 0, 0, 1));
   nbrOffsets.push_back(Offset(-2, 0, 0));
   nbrOffsets.push_back(Offset( 0,-2, 0));
   nbrOffsets.push_back(Offset( 0, 0,-2));
   stencilAverages.addReceives(mpiGrid,nbrOffsets);
   nbrOffsets.clear();

   nbrOffsets.push_back(Offset(-1, 0, 0));
   nbrOffsets.push_back(Offset( 1, 0, 0));
   nbrOffsets.push_back(Offset( 0,-1, 0));
   nbrOffsets.push_back(Offset( 0, 1, 0));
   nbrOffsets.push_back(Offset( 0, 0,-1));
   nbrOffsets.push_back(Offset( 0, 0, 1));
   nbrOffsets.push_back(Offset( 2, 0, 0 ));
   nbrOffsets.push_back(Offset( 0, 2, 0));
   nbrOffsets.push_back(Offset( 0, 0, 2));
   stencilAverages.addSends(mpiGrid,nbrOffsets);
   nbrOffsets.clear();

      // Send/receive stencils for df/dt updates:

   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
       if (i == 0 && j == 0 && k == 0) {
           continue;
       }
       nbrOffsets.push_back(Offset(i,j,k));
   }
   stencilUpdates.addRemoteUpdateReceives(mpiGrid,nbrOffsets);
   stencilUpdates.addRemoteUpdateSends(mpiGrid,nbrOffsets);
   /*
     REPLACED
   nbrTypeIDs.push_back(calcNbrTypeID(2-1,2  ,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2+1,2  ,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2-1,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2+1,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2  ,2-1));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2  ,2+1));
   nbrTypeIDs.push_back(calcNbrTypeID(2-2,2  ,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2-2,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2  ,2-2));
   stencilAverages.addReceives(mpiGrid,nbrTypeIDs);
   
   nbrTypeIDs.clear();
   nbrTypeIDs.push_back(calcNbrTypeID(2-1,2  ,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2+1,2  ,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2-1,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2+1,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2  ,2-1));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2  ,2+1));
   nbrTypeIDs.push_back(calcNbrTypeID(2+2,2  ,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2+2,2  ));
   nbrTypeIDs.push_back(calcNbrTypeID(2  ,2  ,2+2));
   stencilAverages.addSends(mpiGrid,nbrTypeIDs);
   nbrTypeIDs.clear();
   // Send/receive stencils for df/dt updates:

   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
   nbrTypeIDs.push_back(calcNbrTypeID(2+i,2+j,2+k));
   }
   stencilUpdates.addRemoteUpdateReceives(mpiGrid,nbrTypeIDs);
   stencilUpdates.addRemoteUpdateSends(mpiGrid,nbrTypeIDs);
   
   */

   // Allocate receive buffers for all local cells that 
   // have at least one remote neighbour. For GPUs the first 
   // buffer must be allocated using page-locked memory:

   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
       cint host            = it->first.first;
       //cint tag             = it->first.second;
       const CellID localID = it->second;
       Real* buffer = NULL;
       const size_t elements = mpiGrid[localID]->N_blocks*SIZE_VELBLOCK;
       allocateArray(&buffer,elements);
       updateBuffers.insert(make_pair(make_pair(localID,host),buffer));
       remoteUpdates[localID].insert(buffer);
   }
   
   return true;
}

bool finalizeMover() {
   // Free allocated buffers:
   for (map<pair<CellID,int>,Real*>::iterator it=updateBuffers.begin(); it!=updateBuffers.end(); ++it) {
      freeArray(it->second);
   }
   
   return true;
}

void calculateSimParameters(dccrg<SpatialCell>& mpiGrid,creal& t,Real& dt) { }
void calculateCellParameters(dccrg<SpatialCell>& mpiGrid,creal& t,ID::type cell) { }


void calculateAcceleration(dccrg<SpatialCell>& mpiGrid) { 
   typedef Parameters P;
   
   const vector<CellID> cells = mpiGrid.get_cells();

   // Iterate through all local cells and propagate distribution functions 
   // in velocity space. Ghost cells (spatial cells at the boundary of the simulation 
   // volume) do not need to be propagated:
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      if (ghostCells.find(cellID) != ghostCells.end()) continue;
      SpatialCell* SC = mpiGrid[cellID];
      
      // Clear df/dt contributions:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_clearVelFluxes<Real>(*SC,block);
      }
      
      // Calculate df/dt contributions of all blocks in the cell:
      profile::start("df/dt updates in velocity space");
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcVelFluxes<Real>(*SC,block,P::dt,NULL);
      }
      profile::stop("df/dt updates in velocity space");
      
      // Propagate distribution functions in velocity space:
      profile::start("velocity acceleration");
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_propagateVel<Real>(*SC,block,P::dt);
      }
      profile::stop("velocity acceleration");

   }
}

void calculateSpatialDerivatives(dccrg<SpatialCell>& mpiGrid) { }

#ifdef SIMPLE
void calculateSpatialFluxes(dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
   std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/
   
   // TEMPORARY SOLUTION
   vector<CellID> cells = mpiGrid.get_cells();
   cuint avgsByteSize = mpiGrid[cells[0]]->N_blocks * SIZE_VELBLOCK * sizeof(Real);
   // END TEMPORARY SOLUTION

   // Post receives for avgs: 
// REPLACE    mpiGrid.startSingleMode();
   MPIsendRequests.clear();                                                                                                  
   MPIrecvRequests.clear();
   
   
   for (map<pair<int,int>,CellID>::iterator it=stencilAverages.recvs.begin(); it!=stencilAverages.recvs.end(); ++it) {
      cint host           = it->first.first;
      cint tag            = it->first.second;
      const CellID cellID = it->second;
      char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
      cuint byteSize      = avgsByteSize;
      
      //REPLACED mpiGrid.singleReceive(host,tag,byteSize,buffer,cellID);
      MPIrecvRequests.push_back(MPI_Request());
      MPI_Irecv(buffer,byteSize,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back()));
      
   }

   // Apply boundary conditions. This must be done before sending avgs!
   for (set<CellID>::const_iterator it=ghostCells.begin(); it!=ghostCells.end(); ++it) {
      const CellID cellID = *it;
      cuint* const nbrsSpa   = mpiGrid[cellID]->cpu_nbrsSpa;
      cuint existingCells    = nbrsSpa[30];
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
      vlasovBoundaryCondition(cellID,existingCells,nonExistingCells,mpiGrid);
   }
   
   // Post sends for avgs:
   for (multimap<CellID,pair<int,int> >::iterator it=stencilAverages.sends.begin(); it!=stencilAverages.sends.end(); ++it) {
       const CellID cellID = it->first;
       cint host           = it->second.first;
      cint tag            = it->second.second;
      char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
      cuint byteSize      = avgsByteSize;

      //REPLACED       mpiGrid.singleSend(host,tag,byteSize,buffer,cellID);
      MPIsendRequests.push_back(MPI_Request());
      //std::cerr << "ParGrid proc #" << myrank << " MPIsendRequests.size() = " << MPIsendRequests.size() << std::endl;
      if (MPI_Isend(buffer,byteSize,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
          std::cerr << "calculateSpatialFlux failed to send data!" << std::endl;
      }
   }
   
   // Clear spatial fluxes to zero value. Remote neighbour df/dt arrays 
   // need to be cleared as well:
   profile::start("df/dt updates in spatial space");
   cells=mpiGrid.get_cells();
   vector<CellID> remoteCells=mpiGrid.get_list_of_remote_cells_with_local_neighbours();
   cells.insert( cells.end(), remoteCells.begin(), remoteCells.end() );
   
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      Real*  const dfdt   = mpiGrid[cellID]->cpu_fx;      
      for (uint i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) dfdt[i] = 0.0;
   }

   // Iterate through all local cells and calculate their contributions to 
   // time derivatives of distribution functions df/dt in spatial space. Ghost cell 
   // (spatial cells at the boundary of the simulation volume) contributions need to be
   // calculated as well:
   for (set<CellID>::iterator cell=stencilAverages.innerCells.begin(); cell!=stencilAverages.innerCells.end(); ++cell) {
      const CellID cellID      = *cell;
      creal* const avgs        = grid.getAvgs();
      creal* const cellParams  = mpiGrid[cellID]->cpu_cellParams;
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real*  const dfdt        = grid.getFx();
      cuint* const nbrsSpa     = mpiGrid[cellID]->cpu_nbrsSpa;
      
      // Iterate through all velocity blocks in the spatial cell and calculate 
      // contributions to df/dt:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcSpatDfdt(avgs,cellParams,blockParams,dfdt,nbrsSpa,block,P::dt);
      }
   }
   profile::stop("df/dt updates in spatial space");

   // Wait for remote avgs:
   profile::start("(MPI) receive remote averages");
   //REPLACED mpiGrid.waitAllReceives();
   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE);
   // Free memory:
   MPIrecvRequests.clear();
   profile::stop("(MPI) receive remote averages");
   
   // Iterate through the rest of local cells:
   profile::start("df/dt updates in spatial space");
   for (set<CellID>::iterator cell=stencilAverages.boundaryCells.begin(); cell!=stencilAverages.boundaryCells.end(); ++cell) {
      const CellID cellID      = *cell;
      creal* const avgs        = grid.getAvgs();
      creal* const cellParams  = mpiGrid[cellID]->cpu_cellParams;
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real*  const dfdt        = grid.getFx();
      cuint* const nbrsSpa     = mpiGrid[cellID]->cpu_nbrsSpa;
      
      // Iterate through all velocity blocks in the spatial cell and calculate
      // contributions to df/dt:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcSpatDfdt(avgs,cellParams,blockParams,dfdt,nbrsSpa,block,P::dt);
      }
   }
   profile::stop("df/dt updates in spatial space");
   
   // Wait for sends to complete:
   profile::start("(MPI) send averages");
   //REPLACED mpiGrid.waitAllSends();
#ifdef NDEBUG
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),MPI_STATUSES_IGNORE);
#else
   std::vector<MPI_Status> MPIstatuses;
   MPIstatuses.resize(MPIsendRequests.size());
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),&(MPIstatuses[0]));   
   for (uint i=0; i<MPIsendRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS){
       mpilogger <<" Some sends failed for leveque solver "<<endl<<write;
   }
#endif
   MPIsendRequests.clear();
   
   profile::stop("(MPI) send averages");
}

#else // #ifdef SIMPLE
//COMMENT away non-simple version to ease porting work, first the SIMPLE version will be ported to DCCRG
/*
void calculateSpatialFluxes(dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
      
   vector<CellID> cells;
         
   // TEMPORARY SOLUTION
   mpiGrid.getCells(cells);
   cuint avgsByteSize = mpiGrid[cells[0]]->N_blocks * SIZE_VELBLOCK * sizeof(Real);
   const size_t SIZE_DFDT = mpiGrid[cells[0]]->N_blocks * SIZE_VELBLOCK*sizeof(Real);
   // END TEMPORARY SOLUTION

   // Post receives for avgs:
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::iterator it=stencilAverages.recvs.begin(); it!=stencilAverages.recvs.end(); ++it) {
      cint host           = it->first.first;
      cint tag            = it->first.second;
      const CellID cellID = it->second;
      char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
      cuint byteSize      = avgsByteSize;
      mpiGrid.singleReceive(host,tag,byteSize,buffer,cellID);
   }

   // Apply boundary conditions. This must be done before sending avgs!
   for (set<CellID>::const_iterator it=ghostCells.begin(); it!=ghostCells.end(); ++it) {
      const CellID cellID = *it;
      cuint* const nbrsSpa   = mpiGrid[cellID]->cpu_nbrsSpa;
      cuint existingCells    = nbrsSpa[30];
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
      vlasovBoundaryCondition(cellID,existingCells,nonExistingCells,mpiGrid);
   }
   
   // Post sends for avgs:
   for (multimap<CellID,pair<int,int> >::iterator it=stencilAverages.sends.begin(); it!=stencilAverages.sends.end(); ++it) {
      const CellID cellID = it->first;
      cint host           = it->second.first;
      cint tag            = it->second.second;
      char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
      cuint byteSize      = avgsByteSize;
      mpiGrid.singleSend(host,tag,byteSize,buffer,cellID);
   }

   // Post receives for df/dt updates:
   mpiGrid.startSingleMode2();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      const CellID localID  = it->second;
      cint host             = it->first.first;
      cint tag              = it->first.second;
      
      map<pair<CellID,int>,Real*>::iterator it = updateBuffers.find(make_pair(localID,host));
      if (it == updateBuffers.end()) {cerr << "FATAL ERROR: Could not find update buffer!" << endl; exit(1);}
      char* const buffer    = reinterpret_cast<char*>(it->second);
      
      mpiGrid.singleReceive2(host,tag,SIZE_DFDT,buffer,localID);
   }

   // Clear spatial fluxes to zero value. Remote neighbour df/dt arrays need to be cleared as well:
   mpiGrid.getAllCells(cells);
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      Real*  const dfdt   = mpiGrid[cellID]->cpu_fx;
      for (uint i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) dfdt[i] = 0.0;
   }
   
   // Clear the number of received avgs & updates:
   for (map<CellID,pair<uint,uint> >::iterator it=stencilAverages.neighbours.begin(); it!=stencilAverages.neighbours.end(); ++it) {
      it->second.second = 0;
   }
   for (map<CellID,pair<uint,uint> >::iterator it=stencilUpdates.neighbours.begin(); it!=stencilUpdates.neighbours.end(); ++it) {
      it->second.second = 0;
   }
   // Clear the number of computed local updates to remote df/dt:
   for (map<CellID,pair<uint,uint> >::iterator it=stencilUpdates.updates.begin(); it!=stencilUpdates.updates.end(); ++it) {
      it->second.second = 0;
   }

   // Push all inner cells in stencilAverages to readyCells:
   for (set<CellID>::const_iterator it=stencilAverages.innerCells.begin(); it!=stencilAverages.innerCells.end(); ++it) {
      const CellID cellID = *it;
      cuint priority      = stencilUpdates.sends.count(cellID);
      readyCells.insert(cellID,priority);
   }
   
   CellID cellID;
   bool allTasksCompleted;
   size_t calculatedCells = 0;
   uint priority;

   // *****************************************
   // ***** CALCULATE DF/DT CONTRIBUTIONS *****
   // *****************************************
   do {
      allTasksCompleted = true;
      
      // Check if avgs have been received from remote neighbours:
      profile::start("(MPI) receive remote averages");
      mpiGrid.singleModeWaitSome();
      while (mpiGrid.getReadyCell(cellID) == true) {
	 // Increase counter on all local cells requiring the received neighbour data.
	 // If all required data has been received, insert the cell into readyCells:
	 for (multimap<CellID,CellID>::const_iterator it=stencilAverages.remoteToLocalMap.lower_bound(cellID); it!=stencilAverages.remoteToLocalMap.upper_bound(cellID); ++it) {
	    map<CellID,pair<uint,uint> >::iterator localCell = stencilAverages.neighbours.find(it->second);
	    ++(localCell->second.second);
	    if (localCell->second.second == localCell->second.first) {
	       readyCells.insert(localCell->first,stencilUpdates.sends.count(localCell->first));
	    }
	 }
      }
      profile::stop("(MPI) receive remote averages");
      
      // Calculate a cell if possible:
      profile::start("df/dt updates in spatial space");
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 creal* const avgs        = grid.getAvgs();
	 creal* const cellParams  = mpiGrid[cellID]->cpu_cellParams;
	 creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
	 Real*  const dfdt        = grid.getFx();
	 cuint* const nbrsSpa     = mpiGrid[cellID]->cpu_nbrsSpa;


	 // Iterate through all velocity blocks in the spatial cell and calculate
	 // contributions to df/dt:
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_calcSpatDfdt(avgs,cellParams,blockParams,dfdt,nbrsSpa,block,P::dt);
	 }

	 // increase counter on all df/dt updates to remote cells. If all local 
	 // modifications have been calculated, send update to remote neighbour:
	 for (multimap<CellID,CellID>::iterator it=stencilUpdates.remoteToLocalMap.lower_bound(cellID); it!=stencilUpdates.remoteToLocalMap.upper_bound(cellID); ++it) {
	    const CellID nbrID = it->second;
	    ++(stencilUpdates.updates[nbrID].second);
	    if (stencilUpdates.updates[nbrID].second == stencilUpdates.updates[nbrID].first) {
	       multimap<CellID,pair<int,int> >::const_iterator jt=stencilUpdates.sends.find(nbrID);
	       cint host             = jt->second.first;
	       cint tag              = jt->second.second;
	       char* buffer          = reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_fx);
	       mpiGrid.singleSend2(host,tag,SIZE_DFDT,buffer,tag);
	    }
	 }
	 ++calculatedCells;
      }
      profile::stop("df/dt updates in spatial space");
      
      if (calculatedCells != mpiGrid.getNumberOfLocalCells()) allTasksCompleted = false;
   } while (allTasksCompleted == false);
   
   profile::start("(MPI) send averages");
   mpiGrid.singleModeWaitAllSends();
   profile::stop("(MPI) send averages");
   
   // **************************************
   // ***** PROPAGATE IN SPATIAL SPACE *****
   // **************************************
   calculatedCells = 0;
   readyCells.clear();
   for (set<CellID>::const_iterator it=stencilUpdates.innerCells.begin(); it!=stencilUpdates.innerCells.end(); ++it) {
      cellID = *it;
      readyCells.insert(cellID,0);
   }
   
   do {
      allTasksCompleted = true;

      // Check if remote contributions to df/dt have been received. If a local cell has 
      // all required neighbour updates add it to readyCells:
      profile::start("(MPI) receive remote updates");
      mpiGrid.singleModeWaitSome2();
      while (mpiGrid.getReadyCell2(cellID) == true) {
	 ++(stencilUpdates.neighbours[cellID].second);
	 
	 if (stencilUpdates.neighbours[cellID].second == stencilUpdates.neighbours[cellID].first) {
	    // Sum df/dt updates to the first buffer, if necessary:
	    map<CellID,set<Real*> >::iterator it=remoteUpdates.find(cellID);
	    set<Real*>::iterator buffer = it->second.begin();
	    Real* const sumBuffer = *buffer;
	    ++buffer;
	    while (buffer != it->second.end()) {
	       for (uint i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) sumBuffer[i] += (*buffer)[i];
	       ++buffer;
	    }
	    // Insert cell into readyCells:
	    readyCells.insert(cellID,0);
	 }
      }
      profile::stop("(MPI) receive remote updates");

      // Propagate a local cell if possible:
      profile::start("spatial translation");
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 
	 Real* const avgs         = mpiGrid[cellID]->cpu_avgs;
	 creal* const dfdt        = mpiGrid[cellID]->cpu_fx;
	 creal* nbr_dfdt          = NULL;
	 map<CellID,set<Real*> >::const_iterator it=remoteUpdates.find(cellID);
	 if (it != remoteUpdates.end()) nbr_dfdt = *(it->second.begin());
	 creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
	 Real* const cellParams   = mpiGrid[cellID]->cpu_cellParams;
	 
	 // Clear velocity moments that have been calculated during the previous time step:
	 cellParams[CellParams::RHO  ] = 0.0;
	 cellParams[CellParams::RHOVX] = 0.0;
	 cellParams[CellParams::RHOVY] = 0.0;
	 cellParams[CellParams::RHOVZ] = 0.0;


	 if (ghostCells.find(cellID) == ghostCells.end()) {
	    for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	       cpu_propagateSpatWithMoments(avgs,dfdt,nbr_dfdt,blockParams,cellParams,block);
	    }
	 } else {
	    for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	       cpu_calcVelocityMoments(avgs,blockParams,cellParams,block);
	    }
	 }	 

	 ++calculatedCells;
      }
      profile::stop("spatial translation");
      
      if (calculatedCells != mpiGrid.getNumberOfLocalCells()) allTasksCompleted = false;
   } while (allTasksCompleted == false);

   
   profile::start("(MPI) send updates");
   mpiGrid.singleModeWaitAllSends2();
   profile::stop("(MPI) send updates");

}
*/

#endif // #ifdef SIMPLE
   
#ifdef SIMPLE

void calculateSpatialPropagation(dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { 
    std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
    std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/


    vector<CellID> cells;
   // Post receives for remote updates:
   
   // TEMPORARY SOLUTION
   //REPLACED mpiGrid.getCells(cells);
   cells=mpiGrid.get_cells();
   const size_t SIZE_DFDT = mpiGrid[cells[0]]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
   // END TEMPORARY SOLUTION
   
   //REPLACED mpiGrid.startSingleMode();
   MPIsendRequests.clear(); 
   MPIrecvRequests.clear();
   
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      const CellID localID  = it->second;
      cint host             = it->first.first;
      cint tag              = it->first.second;
      map<pair<CellID,int>,Real*>::iterator it2 = updateBuffers.find(make_pair(localID,host));
      if (it2 == updateBuffers.end()) {cerr << "FATAL ERROR: Could not find update buffer!" << endl; exit(1);}
      char* const buffer    = reinterpret_cast<char*>(it2->second);
      
//REPLACED      mpiGrid.singleReceive(host,tag,SIZE_DFDT,buffer,localID);
      MPIrecvRequests.push_back(MPI_Request());
      MPI_Irecv(buffer,SIZE_DFDT,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back()));

   }
   // Post sends for remote updates:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilUpdates.sends.begin(); it!=stencilUpdates.sends.end(); ++it) {
      const CellID nbrID    = it->first;
      cint host             = it->second.first;
      cint tag              = it->second.second;
      char* buffer          = reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_fx);

// REPLACED    mpiGrid.singleSend(host,tag,SIZE_DFDT,buffer,tag);
      MPIsendRequests.push_back(MPI_Request());
      //std::cerr << "ParGrid proc #" << myrank << " MPIsendRequests.size() = " << MPIsendRequests.size() << std::endl;
      if (MPI_Isend(buffer,SIZE_DFDT,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
          std::cerr << "calculateSpatialPropagation failed to send data!" << std::endl;
      }
      
      
   }
   
   profile::start("spatial translation");
   for (set<CellID>::iterator c=stencilUpdates.innerCells.begin(); c!=stencilUpdates.innerCells.end(); ++c) {
      const CellID cellID = *c;
      Real* const avgs         = mpiGrid[cellID]->cpu_avgs;
      creal* const dfdt        = mpiGrid[cellID]->cpu_fx;
      creal* const nbr_dfdt    = NULL;
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real* const cellParams   = mpiGrid[cellID]->cpu_cellParams;
      
      // Clear velocity moments that have been calculated during the previous time step:
      cellParams[CellParams::RHO  ] = 0.0;
      cellParams[CellParams::RHOVX] = 0.0;
      cellParams[CellParams::RHOVY] = 0.0;
      cellParams[CellParams::RHOVZ] = 0.0;
      
      // If the spatial cell is classified as a ghost cell, apply
      // boundary condition before calculating its df/dt contributions:
      if (ghostCells.find(cellID) == ghostCells.end()) {
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_propagateSpatWithMoments(avgs,dfdt,nbr_dfdt,blockParams,cellParams,block);
	 }
      } else {
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_calcVelocityMoments(avgs,blockParams,cellParams,block);
	 }
      }
   }
   profile::stop("spatial translation");
   
   // Wait for remote neighbour updates to arrive:
   profile::start("(MPI) receive remote updates");
    // REPLACED   mpiGrid.waitAllReceives();
   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE);
   // Free memory:
   MPIrecvRequests.clear();
   profile::stop("(MPI) receive remote updates");
   
   // Sum remote neighbour updates to the first buffer of each 
   // local cell, if necessary:
   profile::start("spatial translation");
   for (map<CellID,set<Real*> >::iterator it=remoteUpdates.begin(); it!=remoteUpdates.end(); ++it) {
      const CellID cellID = it->first;
      set<Real*>::iterator buffer = it->second.begin();
      Real* sumBuffer = *buffer;
      ++buffer;
      while (buffer != it->second.end()) {
	 for (uint i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) sumBuffer[i] += (*buffer)[i];
	 ++buffer;
      }
   }
   
   for (set<CellID>::iterator c=stencilUpdates.boundaryCells.begin(); c!=stencilUpdates.boundaryCells.end(); ++c) {
      const CellID cellID = *c;
      Real* const avgs         = mpiGrid[cellID]->cpu_avgs;
      creal* const dfdt        = mpiGrid[cellID]->cpu_fx;
      creal* const nbr_dfdt    = *(remoteUpdates[cellID].begin());
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real* const cellParams   = mpiGrid[cellID]->cpu_cellParams;
      
      cellParams[CellParams::RHO  ] = 0.0;
      cellParams[CellParams::RHOVX] = 0.0;
      cellParams[CellParams::RHOVY] = 0.0;
      cellParams[CellParams::RHOVZ] = 0.0;

      if (ghostCells.find(cellID) == ghostCells.end()) {
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_propagateSpatWithMoments(avgs,dfdt,nbr_dfdt,blockParams,cellParams,block);
	 }
      } else {
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_calcVelocityMoments(avgs,blockParams,cellParams,block);
	 }
      }
   }
   profile::stop("spatial translation");

   // Wait for neighbour update sends:
   profile::start("(MPI) send updates");
// REPLACE   mpiGrid.singleModeWaitAllSends();
#ifdef NDEBUG
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),MPI_STATUSES_IGNORE);
#else
   std::vector<MPI_Status> MPIstatuses;
   MPIstatuses.resize(MPIsendRequests.size());
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),&(MPIstatuses[0]));   
   for (uint i=0; i<MPIsendRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS){
       mpilogger <<" Some sends failed for leveque solver "<<endl<<write;
   }
#endif
   // Free memory:
   MPIsendRequests.clear();
   profile::stop("(MPI) send updates");
}

#else //ifdef SIMPLE
/*
void calculateSpatialPropagation(dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { }
*/

#endif

void initialLoadBalance(dccrg<SpatialCell>& mpiGrid) { }

void calculateVelocityMoments(dccrg<SpatialCell>& mpiGrid) { 
   vector<CellID> cells;
   //REPLCED mpiGrid.getCells(cells);
   cells=mpiGrid.get_cells();
   
   // Iterate through all local cells (incl. ghosts):
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      creal* const avgs        = mpiGrid[cellID]->cpu_avgs;
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real* const cellParams   = mpiGrid[cellID]->cpu_cellParams;
      
      // Clear velocity moments:
      cellParams[CellParams::RHO  ] = 0.0;
      cellParams[CellParams::RHOVX] = 0.0;
      cellParams[CellParams::RHOVY] = 0.0;
      cellParams[CellParams::RHOVZ] = 0.0;
      
      // Iterate through all velocity blocks in this spatial cell 
      // and calculate velocity moments:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcVelocityMoments(avgs,blockParams,cellParams,block);
      }
   }
}


