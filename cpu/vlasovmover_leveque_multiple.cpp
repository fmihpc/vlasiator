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
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include "omp.h"
#endif
#include <boost/mpi.hpp>
#include <zoltan.h>

#include "../vlasovmover.h"
#include "profile.hpp"
#include "../memalloc.h"
#include "cpu_acc_leveque.h"
#include "cpu_trans_leveque.h"
#include "../priorityqueue.h"
#include "../mpilogger.h"
#ifdef CRAY_TOPOLOGY_OPTIMIZATION
#include "graph.h"
#include "mapping.h"
#include "crayxttorus.h"
#endif


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
	const dccrg::Dccrg<SpatialCell>& mpiGrid,
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

bool initializeMover(dccrg::Dccrg<SpatialCell>& mpiGrid) { 
   
   // Populate spatial neighbour list:
   vector<CellID> cells,remoteCells;
   std::vector<CellID> nbrs; //temporary vector for neighbors at certain offset
   
   cells=mpiGrid.get_cells();
   //remoteCells=mpiGrid.get_list_of_remote_cells_with_local_neighbours();
   //cells.insert( cells.end(), remoteCells.begin(), remoteCells.end() );
   
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
              nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2+i, 2+j, 2+k));
              if (nbrIDs.back() == INVALID_CELLID) {
                  isGhost = true;
              }
          }
          ++counter;
      }
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

      if (isGhost == true) mpiGrid[cellID]->isGhostCell = true;
      else mpiGrid[cellID]->isGhostCell = false;
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

   // Exchange ghostFlags between neighbouring processes, 
   // so that boundary condition functions are correctly called 
   // for remote ghost cells:
   SpatialCell::base_address_identifier = 7;
   mpiGrid.start_remote_neighbour_data_update();
   mpiGrid.wait_neighbour_data_update_receives();
   
   // Now iterate through all cells (local + remote), and insert the cells 
   // with isGhostCell flag turned on into ghostCells list. Boundary condition 
   // functions are called for every cell in ghostCells:
   remoteCells=mpiGrid.get_list_of_remote_cells_with_local_neighbours();
   cells.insert( cells.end(), remoteCells.begin(), remoteCells.end() );
   for (uint c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      if (mpiGrid[cellID]->isGhostCell == true) ghostCells.insert(cellID);
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

void calculateSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid,creal& t,Real& dt) { }
void calculateCellParameters(dccrg::Dccrg<SpatialCell>& mpiGrid,creal& t,ID::type cell) { }


void calculateAcceleration(dccrg::Dccrg<SpatialCell>& mpiGrid) { 
   typedef Parameters P;
   
   const vector<CellID> cells = mpiGrid.get_cells();
   vector<CellID> nonGhostCells;

   // Iterate through all local cells and propagate distribution functions 
   // in velocity space. Ghost cells (spatial cells at the boundary of the simulation 
   // volume) do not need to be propagated:
  
   for (size_t c=0; c<cells.size(); ++c) {
       const CellID cellID = cells[c];
       if (ghostCells.find(cellID) != ghostCells.end()) continue;
       nonGhostCells.push_back(cellID);
   }
   //Operations for each cell is local, thus a threaded loop should be safe
#pragma omp parallel for
   for (size_t c=0; c<nonGhostCells.size(); ++c) {
      const CellID cellID = nonGhostCells[c];
      SpatialCell* SC = mpiGrid[cellID];
      
      // Clear df/dt contributions:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
         cpu_clearVelFluxes<Real>(*SC,block);
      }
      
      // Calculate df/dt contributions of all blocks in the cell:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
         cpu_calcVelFluxes<Real>(*SC,block,P::dt,NULL);
      }
      
      // Propagate distribution functions in velocity space:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
         cpu_propagateVel<Real>(*SC,block,P::dt);
      }
   }
   

}

void calculateSpatialDerivatives(dccrg::Dccrg<SpatialCell>& mpiGrid) { }

void calculateSpatialFluxes(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   /*
   // TEMPORARY SOLUTION
   vector<CellID> cells = mpiGrid.get_cells();
   cuint avgsByteSize = mpiGrid[cells[0]]->N_blocks * SIZE_VELBLOCK * sizeof(Real);
   // END TEMPORARY SOLUTION
    */
   profile::start("calculateSpatialFluxes");

   profile::start("Preparation");
   // Prepare vectors with inner and outer cells, ghost cells which are inner and boundary cells, and compute maximum n blocks.
   // These are needed for efficient openMP parallelization, we cannot use iterators
   int maxNblocks=0;
   vector<CellID> ghostInnerCellIds;
   for (set<CellID>::iterator cell=stencilAverages.innerCells.begin(); cell!=stencilAverages.innerCells.end(); ++cell) {
       if(mpiGrid[*cell]->N_blocks>maxNblocks)
           maxNblocks=mpiGrid[*cell]->N_blocks;
       if (ghostCells.find(*cell) != ghostCells.end())
           ghostInnerCellIds.push_back(*cell);
   }

   vector<CellID> ghostBoundaryCellIds;
   for (set<CellID>::iterator cell=stencilAverages.boundaryCells.begin(); cell!=stencilAverages.boundaryCells.end(); ++cell) {
       if(mpiGrid[*cell]->N_blocks>maxNblocks)
           maxNblocks=mpiGrid[*cell]->N_blocks;
       if (ghostCells.find(*cell) != ghostCells.end()) 
           ghostBoundaryCellIds.push_back(*cell);
   }

   //store send/receive map keys to parallelize the MPI loops
   vector<pair<int,int> > receiveHostTag;
   for (map<pair<int,int>,CellID>::iterator it=stencilAverages.recvs.begin(); it!=stencilAverages.recvs.end(); ++it) {
      receiveHostTag.push_back(it->first);
   }

   vector<CellID> sendCellId;
   vector<pair<int,int> > sendHostTag;
   for (multimap<CellID,pair<int,int> >::iterator it=stencilAverages.sends.begin(); it!=stencilAverages.sends.end(); ++it) {
      sendCellId.push_back(it->first);
      sendHostTag.push_back(it->second);
   }       

   profile::stop("Preparation");
//threaded parallel region is entered
#pragma omp parallel    
 {
   int counter;
   //private MPI_Requst tables
   std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
   std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/
   MPIrecvRequests.clear();
   MPIsendRequests.clear();
   
#pragma omp master      
   profile::start("Boundary conditions (inner)");
   // Apply boundary conditions on local ghost cells, these need to be up-to-date for 
   // local cell propagation below:
   //For the openmp parallelization to be safe, only the local cellID should be updated by boundarycondition, and only non-ghost cells should be read.
   
#pragma omp  for  
   for(int i=0;i<ghostInnerCellIds.size();i++){
      const CellID cellID = ghostInnerCellIds[i];
      cuint* const nbrsSpa   = mpiGrid[cellID]->cpu_nbrsSpa;
      cuint existingCells    = nbrsSpa[30];
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
      vlasovBoundaryCondition(cellID,existingCells,nonExistingCells,mpiGrid);
   }
   
#pragma omp master           
   profile::stop("Boundary conditions (inner)",ghostInnerCellIds.size(),"Cells");
   
   // Post receives for avgs. This is done in parallel by all theads (MPI_THREAD_MULTIPLE threadsafety needed)  
#pragma omp master           
   profile::start("(MPI) Start receives");


#pragma omp parallel for
   for(int i=0;i<receiveHostTag.size();i++){
      map<pair<int,int>,CellID>::iterator it=stencilAverages.recvs.find(receiveHostTag[i]);
      cint host           = it->first.first;
      cint tag            = it->first.second;
      const CellID cellID = it->second;
      
      char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
      //cuint byteSize      = avgsByteSize; // NOTE: N_blocks should be ok in buffer cells
      cuint byteSize      = mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
      MPIrecvRequests.push_back(MPI_Request());
      MPI_Irecv(buffer,byteSize,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back()));      
   }
#pragma omp master              
   profile::stop("(MPI) Start receives");
   
   // Post sends for avgs:
   
#pragma omp master              
   profile::start("(MPI) Start sends");

#pragma omp parallel for
   for (int i=0;i<sendCellId.size();i++){
       const CellID cellID = sendCellId[i];
       cint host           = sendHostTag[i].first;
       cint tag            = sendHostTag[i].second;
       char* buffer        = reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs);
       //cuint byteSize      = avgsByteSize; // NOTE: N_blocks should be ok in buffer cells
       cuint byteSize      = mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK*sizeof(Real);
       MPIsendRequests.push_back(MPI_Request());
      
       if (MPI_Isend(buffer,byteSize,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
          std::cerr << "calculateSpatialFlux failed to send data!" << std::endl;
       }
   }

#pragma omp master              
   profile::stop("(MPI) Start sends");
   
   // Clear spatial fluxes to zero value. Remote neighbour df/dt arrays 
   // need to be cleared as well:       

//FIXME, this could be done outside of parallel region
   vector<CellID> cells=mpiGrid.get_cells();
   vector<CellID> remoteCells=mpiGrid.get_list_of_remote_cells_with_local_neighbours();
   cells.insert( cells.end(), remoteCells.begin(), remoteCells.end() );

#pragma omp master
   profile::start("Mark uninit fluxes");

#pragma omp  parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      Real*  const dfdt   = mpiGrid[cellID]->cpu_fx;      
      for (uint i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; i+=SIZE_VELBLOCK) dfdt[i] = numeric_limits<Real>::max() ;
   }
#pragma omp master      
   profile::stop("Mark uninit fluxes");
   
   // Iterate through all local cells and calculate their contributions to 
   // time derivatives of distribution functions df/dt in spatial space. Ghost cell 
   // (spatial cells at the boundary of the simulation volume) contributions need to be
   // calculated as well:
#pragma omp master   
   profile::start("df/dt in real space (inner)");  

   for (set<CellID>::iterator cell=stencilAverages.innerCells.begin(); cell!=stencilAverages.innerCells.end(); ++cell) {
      const CellID cellID      = *cell;
      creal* const avgs        = grid.getAvgs();
      creal* const cellParams  = mpiGrid[cellID]->cpu_cellParams;
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real*  const dfdt        = grid.getFx();
      cuint* const nbrsSpa     = mpiGrid[cellID]->cpu_nbrsSpa;

      // Iterate through all velocity blocks in the spatial cell and calculate 
      // contributions to df/dt:
      //fixme/check, for ghost cells the nbrsSpa can point to the cell itself for non-existing neighbours (see initialization in beginning of this file)
      //this means that for ghost cells fluxes to cell itself may have race condition (does it matter, is its flux even used?)
#pragma omp for
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcSpatDfdt(avgs,cellParams,blockParams,dfdt,nbrsSpa,block,HALF*P::dt);
      }
   }
   
#pragma omp master   
   profile::stop("df/dt in real space (inner)");

   // Wait for remote avgs:
#pragma omp master   
   profile::start("(MPI) Wait receives");
   // Wait for all receives to complete:

   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE);

// all threads sync their Waitall operations
#pragma omp barrier
   // Free memory:
   MPIrecvRequests.clear();
#pragma omp master         
   profile::stop("(MPI) Wait receives");
   
   // Apply boundary conditions on local ghost cells, these need to be up-to-date for 
   // boundary cell propagation below:
#pragma omp master         
   profile::start("Boundary conditions (boundary)");
   //For the openmp parallelization to be safe, only the local cellID should be updated by boundarycondition, and only non-ghost cells should be read.
   //Fixme/check for the copy boundary conditions these are not fulfilled as it does not check for ghost status
#pragma omp  for 
   for(int i=0;i<ghostBoundaryCellIds.size();i++){
       const CellID cellID = ghostBoundaryCellIds[i];
       if (ghostCells.find(cellID) == ghostCells.end()) continue;
       cuint* const nbrsSpa   = mpiGrid[cellID]->cpu_nbrsSpa;
       cuint existingCells    = nbrsSpa[30];
       cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
       vlasovBoundaryCondition(cellID,existingCells,nonExistingCells,mpiGrid);
   }
#pragma omp master            
   profile::stop("Boundary conditions (boundary)",ghostBoundaryCellIds.size(),"Cells");
   
   // Iterate through the rest of local cells:
#pragma omp master         
   profile::start("df/dt in real space (boundary)");

   for (set<CellID>::iterator cell=stencilAverages.boundaryCells.begin(); cell!=stencilAverages.boundaryCells.end(); ++cell) {
      const CellID cellID      = *cell;
      creal* const avgs        = grid.getAvgs();
      creal* const cellParams  = mpiGrid[cellID]->cpu_cellParams;
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real*  const dfdt        = grid.getFx();
      cuint* const nbrsSpa     = mpiGrid[cellID]->cpu_nbrsSpa;
      
      // Iterate through all velocity blocks in the spatial cell and calculate
      // contributions to df/dt:
      // Loop safe as each thread is working on the same cell(?)
      //fixme/check, for ghost cells the nbrsSpa can point to the cell itself for non-existing neighbours (see initialization in beginning of this file)
      //this means that for ghost cells fluxes to cell itself may have race condition. Flux not used in propagation so this should not matter(?)
#pragma omp for
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 cpu_calcSpatDfdt(avgs,cellParams,blockParams,dfdt,nbrsSpa,block,HALF*P::dt);
      }
   }
#pragma omp master
   profile::stop("df/dt in real space (boundary)");
   
   // Wait for sends to complete:
#pragma omp master
   profile::start("(MPI) Wait sends");

   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),MPI_STATUSES_IGNORE);
   MPIsendRequests.clear();
 }//end omp parallel
   profile::stop("(MPI) Wait sends");
   profile::stop("calculateSpatialFluxes");
}

void calculateSpatialPropagation(dccrg::Dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { 
    std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
    std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/

   profile::start("calculateSpatialPropagation");

   vector<CellID> cells;

   vector<CellID> innerCellIds;
   for (set<CellID>::iterator c=stencilUpdates.innerCells.begin(); c!=stencilUpdates.innerCells.end(); ++c) {
       innerCellIds.push_back(*c);
   }

   vector<CellID> boundaryCellIds;
   for (set<CellID>::iterator c=stencilUpdates.boundaryCells.begin(); c!=stencilUpdates.boundaryCells.end(); ++c) {
       boundaryCellIds.push_back(*c);
   }
   
// Post receives for remote updates:
   
   // FIXME: Support variable number of velocity blocks
   cells=mpiGrid.get_cells();
   const size_t SIZE_DFDT = cells.size() > 0 ? mpiGrid[cells[0]]->N_blocks*SIZE_VELBLOCK*sizeof(Real) : 0;
   
   MPIsendRequests.clear(); 
   MPIrecvRequests.clear();
   profile::start("(MPI) Start receives");
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      const CellID localID  = it->second;
      cint host             = it->first.first;
      cint tag              = it->first.second;
      map<pair<CellID,int>,Real*>::iterator it2 = updateBuffers.find(make_pair(localID,host));
      if (it2 == updateBuffers.end()) {cerr << "FATAL ERROR: Could not find update buffer!" << endl; exit(1);}
      char* const buffer    = reinterpret_cast<char*>(it2->second);
      
      MPIrecvRequests.push_back(MPI_Request());
      MPI_Irecv(buffer,SIZE_DFDT,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back()));

   }
   profile::stop("(MPI) Start receives");
   profile::start("(MPI) Start sends");
   // Post sends for remote updates:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilUpdates.sends.begin(); it!=stencilUpdates.sends.end(); ++it) {
      const CellID nbrID    = it->first;
      cint host             = it->second.first;
      cint tag              = it->second.second;
      char* buffer          = reinterpret_cast<char*>(mpiGrid[nbrID]->cpu_fx);

      MPIsendRequests.push_back(MPI_Request());
      //std::cerr << "ParGrid proc #" << myrank << " MPIsendRequests.size() = " << MPIsendRequests.size() << std::endl;
      if (MPI_Isend(buffer,SIZE_DFDT,MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
          std::cerr << "calculateSpatialPropagation failed to send data!" << std::endl;
      }
   }
   profile::stop("(MPI) Start sends");
   
   profile::start("Spatial translation (inner)");
//cpu_propagetSpatWithMoments only write to data in cell cellID, parallel for safe
#pragma parallel for 
   for (int i=0;i<innerCellIds.size();i++){
      const CellID cellID = innerCellIds[i];
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
      
      // Do not propagate ghost cells, only calculate their velocity moments:
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
   profile::stop("Spatial translation (inner)");
   
   // Wait for remote neighbour updates to arrive:
   profile::start("(MPI) Wait receives");
   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE);
   // Free memory:
   MPIrecvRequests.clear();
   profile::stop("(MPI) Wait receives");

      
   // Sum remote neighbour updates to the first receive buffer of each 
   // local cell (if necessary):
   profile::start("Sum remote updates");
   for (map<CellID,set<Real*> >::iterator it=remoteUpdates.begin(); it!=remoteUpdates.end(); ++it) {
      const CellID cellID = it->first;
      set<Real*>::iterator buffer = it->second.begin();
      Real* sumBuffer = *buffer;
      ++buffer;
      while (buffer != it->second.end()) {
#pragma omp parallel for        
	 for (uint i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
             sumBuffer[i] += (*buffer)[i];
	 ++buffer;
      }
   }
   profile::stop("Sum remote updates");
   profile::start("Spatial translation (boundary)");
   // Propagate boundary cells:
//cpu_propagetSpatWithMoments only write to data in cell cellID, parallel for safe
#pragma omp parallel for     
   for (int i=0;i<boundaryCellIds.size();i++){
      const CellID cellID = boundaryCellIds[i];
      Real* const avgs         = mpiGrid[cellID]->cpu_avgs;
      creal* const dfdt        = mpiGrid[cellID]->cpu_fx;
      creal* const nbr_dfdt    = *(remoteUpdates[cellID].begin());
      creal* const blockParams = mpiGrid[cellID]->cpu_blockParams;
      Real* const cellParams   = mpiGrid[cellID]->cpu_cellParams;
      
      cellParams[CellParams::RHO  ] = 0.0;
      cellParams[CellParams::RHOVX] = 0.0;
      cellParams[CellParams::RHOVY] = 0.0;
      cellParams[CellParams::RHOVZ] = 0.0;

      // Do not propagate boundary cells, only calculate their velocity moments:
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
   profile::stop("Spatial translation (boundary)");

   // Wait for neighbour update sends:
   profile::start("(MPI) Wait sends");
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
   profile::stop("(MPI) Wait sends");
   profile::stop("calculateSpatialPropagation");
}

void initialLoadBalance(dccrg::Dccrg<SpatialCell>& mpiGrid) {
  SpatialCell::base_address_identifier = 5;
  mpiGrid.balance_load();

#ifdef CRAY_TOPOLOGY_OPTIMIZATION

  profile::start("Optimize topology mapping");
  //compute first all edges in communication graph so that processes are the vertexes
  int rank,nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
  vector<uint64_t> localCells=mpiGrid.get_cells();
  map<uint64_t,int>  processConnections;
  for(vector<uint64_t>::const_iterator cell=localCells.begin();
      cell!=localCells.end();cell++){
      const vector<uint64_t> *neighbors=mpiGrid.get_neighbours(*cell);
      for(vector<uint64_t>::const_iterator neighbor=neighbors->begin();
          neighbor!=neighbors->end();neighbor++){
          if((*neighbor)>0 && mpiGrid.get_process(*neighbor)!=rank)
              processConnections[mpiGrid.get_process(*neighbor)]++;
      }
  }

  //create torus (hard-coded meteo)
  //torus 1x12x16
  //processors have 6 cores
  CrayXtTorus t(1,12,16,6,MPI_COMM_WORLD);

  //generate graph, add vertexes & edges in grid
  Graph g(MPI_COMM_WORLD);
  g.addVertex(rank);
  for(map<uint64_t,int>::const_iterator edge=processConnections.begin();edge!=processConnections.end();edge++){
      g.addEdge(rank,(int)edge->first,(double)edge->second);
  }
  g.commitChanges();

  //generate map, add local vertex
  std::map<int,int> idMap;
  idMap[rank]=rank; //rank is id for a crayXTtorus graphs
  //create map with seed from clock
  Mapping<Graph,CrayXtTorus> m(g,t,idMap,-1,MPI_COMM_WORLD);

  if(rank==0){
      cout << "Initial Weights " << m.getWeight() <<endl;
  }
  m.simulatedAnnealingOptimizer();
  if(rank==0){
      cout << "Final Weights " << m.getWeight()<<endl;
  }
  profile::stop("Optimize topology mapping");
  profile::start("Migrate cells");
  for(vector<uint64_t>::const_iterator cell=localCells.begin();
      cell!=localCells.end();cell++){
      //pin all local cells to the rank where they should be sen
      if(m.getNetworkVertex(rank)>-1 && m.getNetworkVertex(rank)<nProcs){
          mpiGrid.pin(*cell,m.getNetworkVertex(rank));
      }
      else{
          cout << "ERROR mapping "<<*cell << " " <<rank <<" -> " <<  m.getNetworkVertex(rank)<<endl;
      }
  }

  mpiGrid.migrate_cells();
  mpiGrid.unpin_all_cells();
  profile::stop("Migrate cells");
  
#endif // #ifdef CRAY_TOPOLOGY_OPTIMIZATION 

}

void calculateVelocityMoments(dccrg::Dccrg<SpatialCell>& mpiGrid) { 
   vector<CellID> cells;
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


