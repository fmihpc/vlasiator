// GENERAL NOTES
// 
// "Task queues" seem to be ~50% faster on MPI waiting times on meteo than the "simple" 
// version (-DSIMPLE) in which inner cells are calculated first, then arrival 
// of all neighbour data is waited, and finally boundary cells calculated.

#include <cstdlib>
#include <iostream>
#include <vector>

#include "../vlasovmover.h"
#include "../timer.h"
#include "../memalloc.h"
#include "cpu_acc_leveque.h"
#include "cpu_trans_leveque.h"
#include "../priorityqueue.h"

using namespace std;

#ifdef PARGRID
   #include <limits>
   #include "../pargrid.h"
   #include "../transferstencil.h"

   typedef uint CellID;
   const CellID INVALID_CELLID = numeric_limits<CellID>::max();
   static TransferStencil<CellID> stencilAverages(INVALID_CELLID);
   static TransferStencil<CellID> stencilUpdates(INVALID_CELLID);
#else
   #include <stdint.h>
   #include <dccrg.hpp>
   typedef uint64_t CellID;
#endif

#ifdef PAT_PROFILE
   #include "pat_api.h"
#endif

#ifdef PARGRID

namespace timer {
   uint calcVelFluxes;
   uint calcVelPropag;
   uint calcUpdates;
   uint calcPropag;
   uint sendAverages;
   uint recvAverages;
   uint sendUpdates;
   uint recvUpdates;
   
   uint memoryCopies;
}

static map<pair<CellID,int>,Real*> updateBuffers;
static map<CellID,set<Real*> > remoteUpdates;

static PriorityQueue<CellID> readyCells;

inline uchar calcNbrTypeID(cuchar& i,cuchar& j,cuchar& k) {
   return k*25 + j*5 + i;
}

bool initializeMover(ParGrid<SpatialCell>& mpiGrid) { 
   timer::calcVelFluxes = Timer::create("(COMPUTATION) df/dt updates in velocity space : ");
   timer::calcVelPropag = Timer::create("(COMPUTATION) velocity acceleration           : ");
   timer::calcUpdates   = Timer::create("(COMPUTATION) df/dt updates in spatial space  : ");
   timer::calcPropag    = Timer::create("(COMPUTATION) spatial translation             : ");
   timer::recvAverages  = Timer::create("(    MPI    ) receive remote averages         : ");
   timer::sendAverages  = Timer::create("(    MPI    ) send averages                   : ");
   timer::recvUpdates   = Timer::create("(    MPI    ) receive remote updates          : ");
   timer::sendUpdates   = Timer::create("(    MPI    ) send updates                    : ");
   
   timer::memoryCopies  = Timer::create("( MEM COPY  ) memory copies                   : ");
   
   // Populate spatial neighbour list:
   vector<CellID> cells;
   mpiGrid.getAllCells(cells);
   for (size_t cell=0; cell<cells.size(); ++cell) {
      cuint cellID = cells[cell];
      uint* const nbrsSpa = mpiGrid[cellID]->cpu_nbrsSpa;
      
      // Get spatial neighbour IDs and store them into a vector:
      uint counter = 0;
      vector<uint> nbrIDs;
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
	 if (i == 0 & (j == 0 & k == 0)) nbrIDs.push_back(cellID); // in ParGrid cells do not consider themselves as their own neighbours
	 else nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+i,2+j,2+k)));
	 ++counter;
      }
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(0,2,2))); // i-1,j,k nbr, goes to index 27
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,0,2))); // i,j-2,k nbr, goes to index 28
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,0))); // i,j,k-2 nbr, goes to index 29
      
      // Store neighbour offsets into a vector:
      vector<uint> cellOffsets(nbrIDs.size());
      for (size_t i=0; i<nbrIDs.size(); ++i) {
	 if (nbrIDs[i] == numeric_limits<uint>::max()) cellOffsets[i] = numeric_limits<uint>::max();
	 else cellOffsets[i] = mpiGrid[nbrIDs[i]]->cpuIndex;
      }
      
      // Create spatial neighbour list entry for each block:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 uint boundaryFlag = 0;
	 // Store offsets to each spatial neighbour of this block. Note that
	 // the offset to this block is stored to index 13:
	 for (size_t i=0; i<nbrIDs.size(); ++i) {
	    if (cellOffsets[i] == numeric_limits<uint>::max()) {
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = numeric_limits<uint>::max();
	    } else {
	       boundaryFlag = boundaryFlag | (1 << i);
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = cellOffsets[i] + block;
	    }
	 }
	 // Boundary flags are stored to the last position in nbrsSpa array:
	 nbrsSpa[block*SIZE_NBRS_SPA + 30] = boundaryFlag;
      }
   }
   
   // ***** Calculate MPI send/receive stencils *****

   // Send/receive stencils for avgs:
   vector<uchar> nbrTypeIDs;
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
   
   // Send/receive stencils for df/dt updates:
   nbrTypeIDs.clear();
   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
      nbrTypeIDs.push_back(calcNbrTypeID(2+i,2+j,2+k));
   }
   stencilUpdates.addRemoteUpdateReceives(mpiGrid,nbrTypeIDs);
   stencilUpdates.addRemoteUpdateSends(mpiGrid,nbrTypeIDs);

   // Allocate receive buffers for all local cells that 
   // have at least one remote neighbour. For GPUs the first 
   // buffer must be allocated using page-locked memory:
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      cint host            = it->first.first;
      cint tag             = it->first.second;
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
   Timer::print();
   
   // Free allocated buffers:
   for (map<pair<CellID,int>,Real*>::iterator it=updateBuffers.begin(); it!=updateBuffers.end(); ++it) {
      freeArray(it->second);
   }
   
   return true;
}

void calculateSimParameters(ParGrid<SpatialCell>& mpiGrid,creal& t,Real& dt) { }

void calculateCellParameters(ParGrid<SpatialCell>& mpiGrid,creal& t,ID::type cell) { }

void calculateAcceleration(ParGrid<SpatialCell>& mpiGrid) { 
   typedef Parameters P;
   
   vector<CellID> cells;
   mpiGrid.getCells(cells);

   // Iterate through all local cells and propagate distribution functions 
   // in velocity space. Ghost cells (spatial cells at the boundary of the simulation 
   // volume) do not need to be propagated:
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      
      // Clear df/dt contributions:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_clearVelFluxes<Real>(*SC,block);
      }
      
      #ifdef PAT_PROFILE
         PAT_region_begin(3,"calcVelFluxes");
      #endif
      
      // Calculate df/dt contributions of all blocks in the cell:
      Timer::start(timer::calcVelFluxes);
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcVelFluxes<Real>(*SC,block,P::dt,NULL);
      }
      Timer::stop(timer::calcVelFluxes);
      
      #ifdef PAT_PROFILE
         PAT_region_end(3);
         PAT_region_begin(4,"calcVelPropag");
      #endif
      
      // Propagate distribution functions in velocity space:
      Timer::start(timer::calcVelPropag);
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_propagateVel<Real>(*SC,block,P::dt);
      }
      Timer::stop(timer::calcVelPropag);
      #ifdef PAT_PROFILE
         PAT_region_end(4);
      #endif
   }
}

void calculateSpatialDerivatives(ParGrid<SpatialCell>& mpiGrid) { }

#ifdef SIMPLE
void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   
   vector<CellID> cells;
   
   // TEMPORARY SOLUTION
   mpiGrid.getCells(cells);
   cuint avgsByteSize = mpiGrid[cells[0]]->N_blocks * SIZE_VELBLOCK * sizeof(Real);
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
   
   // Post sends for avgs:
   for (multimap<CellID,pair<int,int> >::iterator it=stencilAverages.sends.begin(); it!=stencilAverages.sends.end(); ++it) {
      const CellID cellID = it->first;
      cint host           = it->second.first;
      cint tag            = it->second.second;
      char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
      cuint byteSize      = avgsByteSize;
      mpiGrid.singleSend(host,tag,byteSize,buffer,cellID);
   }
   
   // Clear spatial fluxes to zero value. Remote neighbour df/dt arrays 
   // need to be cleared as well:
   Timer::start(timer::calcUpdates);
   mpiGrid.getAllCells(cells);
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
   Timer::stop(timer::calcUpdates);

   // Wait for remote avgs:
   Timer::start(timer::recvAverages);
   mpiGrid.waitAllReceives();
   Timer::stop(timer::recvAverages);
   
   // Iterate through the rest of local cells:
   Timer::start(timer::calcUpdates);
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
   Timer::stop(timer::calcUpdates);
   
   // Wait for sends to complete:
   Timer::start(timer::sendAverages);
   mpiGrid.waitAllSends();
   Timer::stop(timer::sendAverages);
   /*
   int procs;
   int myrank;
   MPI_Comm_size(MPI_COMM_WORLD,&procs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
   // Post receives for remote updates:
   
   // TEMPORARY SOLUTION
   const size_t SIZE_DFDT = mpiGrid[cells[0]]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
   // END TEMPORARY SOLUTION
   
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      const CellID localID  = it->second;
      cint host             = it->first.first;
      cint tag              = it->first.second;
      //const size_t byteSize = mpiGrid[localID]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
      //char* buffer          = reinterpret_cast<char*>(updateBuffers[make_pair(localID,make_pair(host,tag))]);
      
      map<pair<CellID,int>,Real*>::iterator it = updateBuffers.find(make_pair(localID,host));
      if (it == updateBuffers.end()) {cerr << "FATAL ERROR: Could not find update buffer!" << endl; exit(1);}
      char* const buffer    = reinterpret_cast<char*>(it->second);
      
      mpiGrid.singleReceive(host,tag,SIZE_DFDT,buffer,localID);
   }

   // Post sends for remote updates:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilUpdates.sends.begin(); it!=stencilUpdates.sends.end(); ++it) {
      const CellID nbrID    = it->first;
      cint host             = it->second.first;
      cint tag              = it->second.second;
      //const size_t byteSize = mpiGrid[cells[0]]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
      //const size_t byteSize = mpiGrid[nbrID]->N_blocks*SIZE_VELBLOCK*sizeof(Real); // TARKISTA
      char* buffer          = reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_fx);
      mpiGrid.singleSend(host,tag,SIZE_DFDT,buffer,tag);
   }
    */
}

#else // #ifdef SIMPLE

void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid) {
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
      Timer::start(timer::recvAverages);
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
      Timer::stop(timer::recvAverages);
      
      // Calculate a cell if possible:
      Timer::start(timer::calcUpdates);
      if (readyCells.empty() == false) {
	 readyCells.pop(cellID,priority);
	 creal* const avgs        = grid.getAvgs();
	 creal* const cellParams  = mpiGrid[cellID]->cpu_cellParams;
	 creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
	 Real*  const dfdt        = grid.getFx();
	 cuint* const nbrsSpa     = mpiGrid[cellID]->cpu_nbrsSpa;

	 #ifdef PAT_PROFILE
	    PAT_region_begin(1,"calcSpatDfdt");
	 #endif
	 
	 // Iterate through all velocity blocks in the spatial cell and calculate
	 // contributions to df/dt:
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_calcSpatDfdt(avgs,cellParams,blockParams,dfdt,nbrsSpa,block,P::dt);
	 }

	 #ifdef PAT_PROFILE
	    PAT_region_end(1);
	 #endif
	 
	 // Increase counter on all df/dt updates to remote cells. If all local 
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
      Timer::stop(timer::calcUpdates);
      
      if (calculatedCells != mpiGrid.getNumberOfLocalCells()) allTasksCompleted = false;
   } while (allTasksCompleted == false);
   
   Timer::start(timer::sendAverages);
   mpiGrid.singleModeWaitAllSends();
   Timer::stop(timer::sendAverages);
   
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
      Timer::start(timer::recvUpdates);
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
      Timer::stop(timer::recvUpdates);

      // Propagate a local cell if possible:
      Timer::start(timer::calcPropag);
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

	 #ifdef PAT_PROFILE
	    PAT_region_begin(2,"propagateWithMoments");
	 #endif
	 
	 for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	    cpu_propagateSpatWithMoments(avgs,dfdt,nbr_dfdt,blockParams,cellParams,block);
	 }
	 
	 #ifdef PAT_PROFILE
	    PAT_region_end(2);
	 #endif
	 
	 // Divide calculated velocity moments by spatial cell volume to get SI units:
	 creal VOLUME = cellParams[CellParams::DX]*cellParams[CellParams::DY]*cellParams[CellParams::DZ];
	 cellParams[CellParams::RHO  ] /= VOLUME;
	 cellParams[CellParams::RHOVX] /= VOLUME;
	 cellParams[CellParams::RHOVY] /= VOLUME;
	 cellParams[CellParams::RHOVZ] /= VOLUME;
	 
	 ++calculatedCells;
      }
      Timer::stop(timer::calcPropag);
      
      if (calculatedCells != mpiGrid.getNumberOfLocalCells()) allTasksCompleted = false;
   } while (allTasksCompleted == false);

   
   Timer::start(timer::sendUpdates);
   mpiGrid.singleModeWaitAllSends2();
   Timer::stop(timer::recvUpdates);

}

#endif // #ifdef SIMPLE
   
#ifdef SIMPLE

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { 
   vector<CellID> cells;
   // Post receives for remote updates:
   
   // TEMPORARY SOLUTION
   mpiGrid.getCells(cells);
   const size_t SIZE_DFDT = mpiGrid[cells[0]]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
   // END TEMPORARY SOLUTION
   
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      const CellID localID  = it->second;
      cint host             = it->first.first;
      cint tag              = it->first.second;
      map<pair<CellID,int>,Real*>::iterator it = updateBuffers.find(make_pair(localID,host));
      if (it == updateBuffers.end()) {cerr << "FATAL ERROR: Could not find update buffer!" << endl; exit(1);}
      char* const buffer    = reinterpret_cast<char*>(it->second);
      
      mpiGrid.singleReceive(host,tag,SIZE_DFDT,buffer,localID);
   }
   // Post sends for remote updates:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilUpdates.sends.begin(); it!=stencilUpdates.sends.end(); ++it) {
      const CellID nbrID    = it->first;
      cint host             = it->second.first;
      cint tag              = it->second.second;
      char* buffer          = reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_fx);
      mpiGrid.singleSend(host,tag,SIZE_DFDT,buffer,tag);
   }
      
   Timer::start(timer::calcPropag);
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
      
      
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_propagateSpatWithMoments(avgs,dfdt,nbr_dfdt,blockParams,cellParams,block);
      }
      
      // Divide calculated velocity moments by spatial cell volume 
      // in order to get SI units:
      creal VOLUME = cellParams[CellParams::DX]*cellParams[CellParams::DY]*cellParams[CellParams::DZ];
      cellParams[CellParams::RHO  ] /= VOLUME;
      cellParams[CellParams::RHOVX] /= VOLUME;
      cellParams[CellParams::RHOVY] /= VOLUME;
      cellParams[CellParams::RHOVZ] /= VOLUME;
   }
   Timer::stop(timer::calcPropag);

   // Wait for remote neighbour updates to arrive:
   Timer::start(timer::recvUpdates);
   mpiGrid.waitAllReceives();
   Timer::stop(timer::recvUpdates);
   
   // Sum remote neighbour updates to the first buffer of each 
   // local cell, if necessary:
   Timer::start(timer::calcPropag);
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
      
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_propagateSpatWithMoments(avgs,dfdt,nbr_dfdt,blockParams,cellParams,block);
      }
      
      creal VOLUME = cellParams[CellParams::DX]*cellParams[CellParams::DY]*cellParams[CellParams::DZ];
      cellParams[CellParams::RHO  ] /= VOLUME;
      cellParams[CellParams::RHOVX] /= VOLUME;
      cellParams[CellParams::RHOVY] /= VOLUME;
      cellParams[CellParams::RHOVZ] /= VOLUME;
   }
   Timer::stop(timer::calcPropag);

   // Wait for neighbour update sends:
   Timer::start(timer::sendUpdates);
   mpiGrid.singleModeWaitAllSends();
   Timer::stop(timer::sendUpdates);
}

#else

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { }

#endif

void initialLoadBalance(ParGrid<SpatialCell>& mpiGrid) { }

void calculateVelocityMoments(ParGrid<SpatialCell>& mpiGrid) { 
   vector<CellID> cells;
   mpiGrid.getCells(cells);
   
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
      
      // Divide calculated velocity moments by spatial cell volume 
      // in order to get SI units:
      creal VOLUME = cellParams[CellParams::DX]*cellParams[CellParams::DY]*cellParams[CellParams::DZ];
      cellParams[CellParams::RHO  ] /= VOLUME;
      cellParams[CellParams::RHOVX] /= VOLUME;
      cellParams[CellParams::RHOVY] /= VOLUME;
      cellParams[CellParams::RHOVZ] /= VOLUME;
   }
}

#endif
