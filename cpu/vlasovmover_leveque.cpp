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

#include <stdint.h>
#define DCCRG_SEND_SINGLE_CELLS
#define DCCRG_CELL_DATA_SIZE_FROM_USER
#define DCCRG_USER_MPI_DATA_TYPE
#include <dccrg.hpp>
typedef uint64_t CellID;

#include "../transferstencil.h"
#include "spatial_cell.hpp"

using namespace std;
using namespace spatial_cell;

static TransferStencil<CellID> stencilAverages(INVALID_CELLID);
static TransferStencil<CellID> stencilUpdates(INVALID_CELLID);

static set<CellID> ghostCells;


static map<pair<CellID,int>,Real*> updateBuffers; /**< For each local cell receiving one or more remote df/dt updates,
						   * MPI rank of remote process sending an update and address to the 
						   * allocated buffer. */
static map<CellID,vector< vector<Real> > > remoteUpdates;     /**< For each local cell receiving one or more remote df/dt updates, 
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
	const dccrg::Dccrg<SpatialCell>& mpiGrid,
	const CellID& cellID,
	const uchar& i,
	const uchar& j,
	const uchar& k
) {
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
      SpatialCell *SC = mpiGrid[cellID];

      SC->neighbors.clear();
      SC->isGhostCell = false;
      SC->boundaryFlag=0;
      // Get spatial neighbour IDs and store them into a vector:
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
          // in dccrg cells are neighbors of themselves only with periodic boundaries
         if (i == 0 && j == 0 && k == 0) {
            SC->neighbors.push_back(cellID);
         }
         else {
            SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 2+i, 2+j, 2+k));
            if (SC->neighbors.back() == INVALID_CELLID) {
               //we only use a stencil of one to set if a cell is a ghost cell.
               SC->isGhostCell = true;
            }
         }
      }      

      if (getNeighbourID(mpiGrid, cellID, 0, 2, 2) != INVALID_CELLID) {
         SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 0, 2, 2));
      }
      else{
         //outside boundary, attempt to use one closer 
         SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 1, 2, 2));
      }

      if (getNeighbourID(mpiGrid, cellID, 2, 0, 2) != INVALID_CELLID) {
         SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 2, 0, 2));
      }
      else{
         //outside boundary, attempt to use one closer 
         SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 2, 1, 2));
      }

      if (getNeighbourID(mpiGrid, cellID, 2, 2, 0) != INVALID_CELLID) {
         SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 2, 2, 0));
      }
      else{
         //outside boundary, attempt to use one closer 
         SC->neighbors.push_back(getNeighbourID(mpiGrid, cellID, 2, 2, 1));
      }
                  
      
      for(int i=0;i<SC->neighbors.size();i++){
         if (SC->neighbors[i] == INVALID_CELLID) {
      // Offsets to missing 
      // neighbours are replaced by offsets to this cell so that cells on the 
      // boundary of the simulation domain (="ghost cells") can be propagated.
      // This allows boundary values (set by user) to propagate into simulation 
      // domain.
            SC->neighbors[i]=cellID;
         }
         else
            SC->boundaryFlag = SC->boundaryFlag | (1 << i);
      }
   }
   
   // ***** Calculate MPI send/receive stencils *****
   // Send/receive stencils for avgs:
   stencilAverages.clear();

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

   stencilUpdates.clear();
   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
      if (i == 0 && j == 0 && k == 0) {
           continue;
      }
      nbrOffsets.push_back(Offset(i,j,k));
   }
   stencilUpdates.addRemoteUpdateReceives(mpiGrid,nbrOffsets);
   stencilUpdates.addRemoteUpdateSends(mpiGrid,nbrOffsets);


   // Exchange ghostFlags between neighbouring processes, 
   // so that boundary condition functions are correctly called 
   // for remote ghost cells:
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_GHOSTFLAG);
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
      calculateCellAcceleration(mpiGrid,cellID);
   }
}

void calculateCellAcceleration(dccrg::Dccrg<SpatialCell>& mpiGrid,CellID cellID) {
   profile::start("Acceleration");
   SpatialCell* SC = mpiGrid[cellID];
   if (ghostCells.find(cellID) != ghostCells.end()){
      profile::stop("Acceleration",0,"Blocks");
      return;
   }
   profile::start("clearVelFluxes");
   // Clear df/dt contributions:
   for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
      unsigned int block = SC->velocity_block_list[block_i];         
      cpu_clearVelFluxes(SC,block);
   }
   profile::stop("clearVelFluxes");

   profile::start("calcVelFluxes");
   // Calculatedf/dt contributions of all blocks in the cell:
   for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
      unsigned int block = SC->velocity_block_list[block_i];         
      cpu_calcVelFluxes(SC,block,P::dt,NULL);
   }
   profile::stop("calcVelFluxes");

   profile::start("propagateVel");
      // Propagate distribution functions in velocity space:
   for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
      unsigned int block = SC->velocity_block_list[block_i];         
      cpu_propagateVel(SC,block,P::dt);
   }
   profile::stop("propagateVel");
   profile::stop("Acceleration",SC->number_of_blocks,"Blocks");
}

void calculateSpatialDerivatives(dccrg::Dccrg<SpatialCell>& mpiGrid) { }

void calculateSpatialFluxes(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   int counter;
   std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
   std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/
   std::vector<MPI_Datatype> MPIrecvTypes;               /**< Container for active datatypes due to receives    .*/
   std::vector<MPI_Datatype> MPIsendTypes;               /**< Container for active datatypes due to sends.*/
   
   /*
   // TEMPORARY SOLUTION
   vector<CellID> cells = mpiGrid.get_cells();
   cuint avgsByteSize = mpiGrid[cells[0]]->N_blocks * SIZE_VELBLOCK * sizeof(Real);
   // END TEMPORARY SOLUTION
    */
   profile::start("calculateSpatialFluxes");
   MPIsendRequests.clear(); // Note: unnecessary
   MPIrecvRequests.clear(); // Note: unnecessary

   // Prepare vectors with inner and outer cells, ghost cells which are inner and boundary cells, and compute maximum n blocks.
   // These are needed for efficient openMP parallelization, we cannot use iterators
   int maxNblocks=0;
   vector<CellID> ghostInnerCellIds;
   for (set<CellID>::iterator cell=stencilAverages.innerCells.begin(); cell!=stencilAverages.innerCells.end(); ++cell) {
       if(mpiGrid[*cell]->size()>maxNblocks)
           maxNblocks=mpiGrid[*cell]->size();
       if (ghostCells.find(*cell) != ghostCells.end())
           ghostInnerCellIds.push_back(*cell);
   }

   vector<CellID> ghostBoundaryCellIds;
   for (set<CellID>::iterator cell=stencilAverages.boundaryCells.begin(); cell!=stencilAverages.boundaryCells.end(); ++cell) {
       if(mpiGrid[*cell]->size()>maxNblocks)
           maxNblocks=mpiGrid[*cell]->size();
       if (ghostCells.find(*cell) != ghostCells.end()) 
           ghostBoundaryCellIds.push_back(*cell);
   }

   
   
   profile::start("Boundary conditions (inner)");
   // Apply boundary conditions on local ghost cells, these need to be up-to-date for 
   // local cell propagation below:
   //For the openmp parallelization to be safe, only the local cellID should be updated by boundarycondition, and only non-ghost cells should be read.
#pragma omp parallel for  
   for(int i=0;i<ghostInnerCellIds.size();i++){
      const CellID cellID = ghostInnerCellIds[i];
//      cuint* const nbrsSpa   = mpiGrid[cellID]->neighbors;
      cuint existingCells    = mpiGrid[cellID]->boundaryFlag;
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
      vlasovBoundaryCondition(cellID,existingCells,nonExistingCells,mpiGrid);
      counter++;
   }
   profile::stop("Boundary conditions (inner)",ghostInnerCellIds.size(),"Cells");
   
   // Post receives for avgs:
   profile::initializeTimer("Start receives","MPI");
   profile::start("Start receives");
   for (map<pair<int,int>,CellID>::iterator it=stencilAverages.recvs.begin(); it!=stencilAverages.recvs.end(); ++it) {
      cint host           = it->first.first;
      cint tag            = it->first.second;
      const CellID cellID = it->second;
      //cuint byteSize      = avgsByteSize; // NOTE: N_blocks should be ok in buffer cells
      mpiGrid[cellID]->set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      MPIrecvTypes.push_back(mpiGrid[cellID]->mpi_datatype());
      MPI_Type_commit(&(MPIrecvTypes.back()));
      MPIrecvRequests.push_back(MPI_Request());
      MPI_Irecv(mpiGrid[cellID],1,MPIrecvTypes.back(),host,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back()));      
   }
   profile::stop("Start receives");

   // Post sends for avgs:
   profile::initializeTimer("Start sends","MPI");
   profile::start("Start sends");
   for (multimap<CellID,pair<int,int> >::iterator it=stencilAverages.sends.begin(); it!=stencilAverages.sends.end(); ++it) {
      const CellID cellID = it->first;
      cint host           = it->second.first;
      cint tag            = it->second.second;
      mpiGrid[cellID]->set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      MPIsendTypes.push_back(mpiGrid[cellID]->mpi_datatype());
      MPI_Type_commit(&(MPIsendTypes.back()));

      MPIsendRequests.push_back(MPI_Request());      
      if (MPI_Isend(mpiGrid[cellID],1,MPIsendTypes.back(),host,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
         std::cerr << "calculateSpatialFlux failed to send data!" << std::endl;
      }
   }
   profile::stop("Start sends");
   
   // Clear spatial fluxes to zero value. Remote neighbour df/dt arrays 
   // need to be cleared as well:       
   vector<CellID> cells=mpiGrid.get_cells();
   vector<CellID> remoteCells=mpiGrid.get_list_of_remote_cells_with_local_neighbours();
   cells.insert( cells.end(), remoteCells.begin(), remoteCells.end() );
   
   profile::start("zero fluxes");
#pragma omp  parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];

      SpatialCell* SC=mpiGrid[cellID];
      for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
         unsigned int block = SC->velocity_block_list[block_i];         
         Velocity_Block* block_ptr = SC->at(block);     
         for (uint i=0; i<SIZE_VELBLOCK; i++) block_ptr->fx[i]=0.0;
      }
   }
   profile::stop("zero fluxes");
   
   // Iterate through all local cells and calculate their contributions to 
   // time derivatives of distribution functions df/dt in spatial space. Ghost cell 
   // (spatial cells at the boundary of the simulation volume) contributions need to be
   // calculated as well:
   profile::start("df/dt in real space (inner)");  

#pragma omp parallel
   
   for (set<CellID>::iterator cell=stencilAverages.innerCells.begin(); cell!=stencilAverages.innerCells.end(); ++cell) {
      SpatialCell* SC = mpiGrid[*cell];
      // Iterate through all velocity blocks in the spatial cell and calculate 
      // contributions to df/dt:
      //fixme/check, for ghost cells the nbrsSpa can point to the cell itself for non-existing neighbours (see initialization in beginning of this file)
      //this means that for ghost cells fluxes to cell itself may have race condition (does it matter, is its flux even used?)
//#pragma omp for
      for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
         unsigned int block = SC->velocity_block_list[block_i];         
	 cpu_calcSpatDfdt(mpiGrid,SC,block,HALF*P::dt);
      }
   }
   profile::stop("df/dt in real space (inner)");

   // Wait for remote avgs:
   profile::initializeTimer("Wait receives","MPI","Wait");
   profile::start("Wait receives");
   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE);
   // Free memory:
   MPIrecvRequests.clear();
   for(int i=0;i<MPIrecvTypes.size();i++){
      MPI_Type_free(&(MPIrecvTypes[i]));
   }
   MPIrecvTypes.clear();
   profile::stop("Wait receives");
   
   // Apply boundary conditions on local ghost cells, these need to be up-to-date for 
   // boundary cell propagation below:

   profile::start("Boundary conditions (boundary)");
   //For the openmp parallelization to be safe, only the local cellID should be updated by boundarycondition, and only non-ghost cells should be read.
   //Fixme/check for the copy boundary conditions these are not fulfilled as it does not check for ghost status
#pragma omp parallel for 
   for(int i=0;i<ghostBoundaryCellIds.size();i++){
       const CellID cellID = ghostBoundaryCellIds[i];
       if (ghostCells.find(cellID) == ghostCells.end()) continue;
       //cuint* const nbrsSpa   = mpiGrid[cellID]->cpu_nbrsSpa;
       cuint existingCells    = mpiGrid[cellID]->boundaryFlag;
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
//       cuint existingCells    = nbrsSpa[30];
      //   cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
       vlasovBoundaryCondition(cellID,existingCells,nonExistingCells,mpiGrid);
   }
   profile::stop("Boundary conditions (boundary)",ghostBoundaryCellIds.size(),"Cells");
   
   // Iterate through the rest of local cells:
   profile::start("df/dt in real space (boundary)");

#pragma omp parallel 
   for (set<CellID>::iterator cell=stencilAverages.boundaryCells.begin(); cell!=stencilAverages.boundaryCells.end(); ++cell) {
      const CellID cellID      = *cell;
      SpatialCell* SC = mpiGrid[cellID];
      // Iterate through all velocity blocks in the spatial cell and calculate 
      // contributions to df/dt:
      //fixme/check, for ghost cells the nbrsSpa can point to the cell itself for non-existing neighbours (see initialization in beginning of this file)
      //this means that for ghost cells fluxes to cell itself may have race condition (does it matter, is its flux even used?)
#pragma omp for
      for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
         unsigned int block = SC->velocity_block_list[block_i];         
	 cpu_calcSpatDfdt(mpiGrid,SC,block,HALF*P::dt);
      }


   }
   profile::stop("df/dt in real space (boundary)");
   
   // Wait for sends to complete:
   profile::initializeTimer("Wait sends","MPI","Wait");
   profile::start("Wait sends");

   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),MPI_STATUSES_IGNORE);

   MPIsendRequests.clear();
   for(int i=0;i<MPIsendTypes.size();i++){
      MPI_Type_free(&(MPIsendTypes[i]));
   }
   MPIsendTypes.clear();
   profile::stop("Wait sends");
   profile::stop("calculateSpatialFluxes");
}

void calculateSpatialPropagation(dccrg::Dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { 
   std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
   std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/
   std::vector<MPI_Datatype> MPIsendTypes;               /**< Container for active datatypes due to sends.*/
   
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
   

   cells=mpiGrid.get_cells();
   
   MPIsendRequests.clear(); 
   MPIrecvRequests.clear();
   MPIsendTypes.clear();
   profile::initializeTimer("Start receives","MPI");
   profile::start("Start receives");


   
   // Allocate receive buffers for all local cells that 
   // have at least one remote neighbour. For GPUs the first 
   // buffer must be allocated using page-locked memory:
   //FIXME, avoid excessiv deallocation/allocation
   remoteUpdates.clear();
   updateBuffers.clear();
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      cint host            = it->first.first;
      const CellID localID = it->second;
      const size_t elements = mpiGrid[localID]->size()*SIZE_VELBLOCK;
      remoteUpdates[localID].push_back(vector<Real>(elements));
      updateBuffers.insert(make_pair(make_pair(localID,host), &(remoteUpdates[localID].back()[0]) ));      
   }

   
   for (map<pair<int,int>,CellID>::const_iterator it=stencilUpdates.recvs.begin(); it!=stencilUpdates.recvs.end(); ++it) {
      const CellID localID  = it->second;
      cint host             = it->first.first;
      cint tag              = it->first.second;

      
      map<pair<CellID,int>,Real*>::iterator it2 = updateBuffers.find(make_pair(localID,host));
      if (it2 == updateBuffers.end()) {cerr << "FATAL ERROR: Could not find update buffer!" << endl; exit(1);}
      char* const buffer    = reinterpret_cast<char*>(it2->second);

      //receive as bytestream (convert to sparse format later on)       
      MPIrecvRequests.push_back(MPI_Request());
      MPI_Irecv(buffer, mpiGrid[localID]->size()*SIZE_VELBLOCK*sizeof(Real),MPI_BYTE,host,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back()));

   }
   profile::stop("Start receives");
   profile::initializeTimer("Start sends","MPI");
   profile::start("Start sends");
   // Post sends for remote updates:
   for (multimap<CellID,pair<int,int> >::const_iterator it=stencilUpdates.sends.begin(); it!=stencilUpdates.sends.end(); ++it) {
      const CellID nbrID    = it->first;
      cint host             = it->second.first;
      cint tag              = it->second.second;
      mpiGrid[nbrID]->set_mpi_transfer_type(Transfer::VEL_BLOCK_FLUXES);
      MPIsendTypes.push_back(mpiGrid[nbrID]->mpi_datatype());
      MPI_Type_commit(&(MPIsendTypes.back()));
      MPIsendRequests.push_back(MPI_Request());
      
      if (MPI_Isend(mpiGrid[nbrID],1,MPIsendTypes.back(),host,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
          std::cerr << "calculateSpatialPropagation failed to send data!" << std::endl;
      }
      
   }
   profile::stop("Start sends");
   
   profile::start("Spatial translation (inner)");
//cpu_propagetSpatWithMoments only write to data        in cell cellID, parallel for safe
#pragma omp parallel for        
   for (int i=0;i<innerCellIds.size();i++){
      const CellID cellID = innerCellIds[i];

      creal* const nbr_dfdt    = NULL;
      SpatialCell* SC = mpiGrid[cellID];
      // Clear velocity moments that have been calculated during the previous time step:
      SC->parameters[CellParams::RHO  ] = 0.0;
      SC->parameters[CellParams::RHOVX] = 0.0;
      SC->parameters[CellParams::RHOVY] = 0.0;
      SC->parameters[CellParams::RHOVZ] = 0.0;

      
      // Do not propagate boundary cells, only calculate their velocity moments:
      if (ghostCells.find(cellID) == ghostCells.end()) {
         for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
            unsigned int block = SC->velocity_block_list[block_i];         
            cpu_propagateSpatWithMoments(nbr_dfdt,SC,block,block_i);
         }  
      } else {
         for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
            unsigned int block = SC->velocity_block_list[block_i];         
	    cpu_calcVelocityMoments(SC,block);
	 }
      }
      if(!secondStep)
         calculateCellAcceleration(mpiGrid,cellID);
   }
   profile::stop("Spatial translation (inner)");
   
   // Wait for remote neighbour updates to arrive:
   profile::initializeTimer("Wait receives","MPI","Wait");
   profile::start("Wait receives");
   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE);
   // Free memory:
   MPIrecvRequests.clear();
   profile::stop("Wait receives");

      
   // Sum remote neighbour updates to the first receive buffer of each 
   // local cell (if necessary):
   //FIXME, is this old comment taken care of?:  NEED TO SUM THE DENSE TEMPORARY ARRAY FLUXES INTO THE SPARSE ONE. ORDER AS IN THE CREATION OF THE DATATYPE (....)
   profile::start("Sum remote updates");
   for (map<CellID,vector<vector<Real> > >::iterator it=remoteUpdates.begin(); it!=remoteUpdates.end(); ++it) {
      const CellID cellID = it->first;
      vector<vector<Real> >::iterator buffer = it->second.begin();
      //sum results into first receive buffer if is not empty
      if(buffer != it->second.end()) {
         Real* sumBuffer = &((*buffer)[0]);
         ++buffer;
         while (buffer != it->second.end()) {
            for (uint i=0; i< (*buffer).size();i++)
               sumBuffer[i] += (*buffer)[i];
            ++buffer;
         }
      }
   }
   profile::stop("Sum remote updates");
   profile::start("Spatial translation (boundary)");
   // Propagate boundary cells:
//cpu_propagetSpatWithMoments only write to data in cell cellID, parallel for safe

   for (int i=0;i<boundaryCellIds.size();i++){
      const CellID cellID = boundaryCellIds[i];
      creal* const nbr_dfdt    = &(remoteUpdates[cellID][0][0]);
      SpatialCell* SC = mpiGrid[cellID];
      SC->parameters[CellParams::RHO  ] = 0.0;
      SC->parameters[CellParams::RHOVX] = 0.0;
      SC->parameters[CellParams::RHOVY] = 0.0;
      SC->parameters[CellParams::RHOVZ] = 0.0;

      
      // Do not propagate boundary cells, only calculate their velocity moments:
      if (ghostCells.find(cellID) == ghostCells.end()) {
         for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
            unsigned int block = SC->velocity_block_list[block_i];         
            cpu_propagateSpatWithMoments(nbr_dfdt,SC,block,block_i);
         }  
      } else {
         for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
            unsigned int block = SC->velocity_block_list[block_i];         
	    cpu_calcVelocityMoments(SC,block);
	 }
      }
      if(!secondStep)
         calculateCellAcceleration(mpiGrid,cellID);
   }
   
   profile::stop("Spatial translation (boundary)");

   // Wait for neighbour update sends:
   profile::initializeTimer("Wait sends","MPI","Wait");
   profile::start("Wait sends");

   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),MPI_STATUSES_IGNORE);

   // Free memory:
   MPIsendRequests.clear();
   for(int i=0;i<MPIsendTypes.size();i++){
      MPI_Type_free(&(MPIsendTypes[i]));
   }
   MPIsendTypes.clear();
   
   profile::stop("Wait sends");
   profile::stop("calculateSpatialPropagation");
}








void calculateVelocityMoments(dccrg::Dccrg<SpatialCell>& mpiGrid) { 
   vector<CellID> cells;
   cells=mpiGrid.get_cells();
   
   // Iterate through all local cells (incl. ghosts):
#pragma omp parallel for        
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      // Clear velocity moments:
      SC->parameters[CellParams::RHO  ] = 0.0;
      SC->parameters[CellParams::RHOVX] = 0.0;
      SC->parameters[CellParams::RHOVY] = 0.0;
      SC->parameters[CellParams::RHOVZ] = 0.0;
      
      for(unsigned int block_i=0; block_i< SC->number_of_blocks;block_i++){
         unsigned int block = SC->velocity_block_list[block_i];         
      
      // Iterate through all velocity blocks in this spatial cell 
      // and calculate velocity moments:        
         cpu_calcVelocityMoments(SC,block);
      }
   }
}


      


