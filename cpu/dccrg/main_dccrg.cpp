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

#ifndef PARGRID // Do not use these functions if ParGrid is used:

#include <cstdlib>
#include <iostream>
#include <vector>
#include <boost/mpi.hpp>
#include <zoltan.h>

#define DCCRG_SEND_SINGLE_CELLS
#define DCCRG_CELL_DATA_SIZE_FROM_USER
#define DCCRG_USER_MPI_DATA_TYPE
#include <dccrg.hpp>

#include <mpilogger.h>
#include <definitions.h>
#include <parameters.h>
#include <cell_spatial.h>
#include <project.h>
#include "vlasovmover.h"
#include <profile.h>

extern MPILogger mpilogger;

#ifdef PARGRID
typedef uint CellID;
#else
typedef uint64_t CellID;
#endif

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_calcVelocityMoments(SpatialCell& cell);

namespace Main {
   std::vector<CellID> cells;
  std::vector<const SpatialCell*> nbrPtrs(6,(SpatialCell *)NULL);
   SpatialCell* cellPtr;
   
   uint calcAcc;
   uint calcSpatDerivs;
   uint spatDerivsMPIRecv;
   uint spatDerivsMPISend;
   uint calcSpatFluxes;
   uint spatFluxesMPIRecv;
   uint spatFluxesMPISend;
   uint calcSpatProp;
   uint spatPropMPIRecv;
   uint spatPropMPISend;
}

bool finalizeMover() {
   return true;
}

// TODO these are almost identical to the ones in ldz.cpp, merge
inline uchar calcNbrNumber(const uchar& i,const uchar& j,const uchar& k) {return k*9+j*3+i;}

inline uchar calcNbrTypeID(const uchar& i,const uchar& j,const uchar& k) {return k*25+j*5+i;}

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
   // TODO: merge this with the one in lond...anna.cpp and probably in ...mover_leveque.cpp
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
   Main::cells = mpiGrid.get_cells();
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      CellID cellID = Main::cells[cell];
      uint* const nbrsSpa = mpiGrid[cellID]->cpu_nbrsSpa;

      // Get spatial neighbour IDs and store them into a vector:
      uint counter = 0;
      std::vector<uint> nbrIDs;
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2 - 1, 2    , 2    ));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2 + 1, 2    , 2    ));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2    , 2 - 1, 2    ));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2    , 2 + 1, 2    ));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2    , 2    , 2 - 1));
      nbrIDs.push_back(getNeighbourID(mpiGrid, cellID, 2    , 2    , 2 + 1));

      // Store neighbour offsets into a vector:
      std::vector<uint> cellOffsets(nbrIDs.size());
      for (size_t i=0; i<nbrIDs.size(); ++i) {
	 if (nbrIDs[i] == INVALID_CELLID) cellOffsets[i] = INVALID_CELLID;
	 else cellOffsets[i] = mpiGrid[nbrIDs[i]]->cpuIndex * SIZE_VELBLOCK;
      }
      
      // Create spatial neighbour list entry for each block:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 uint boundaryFlag = 0;
	 // Store offsets to each spatial neighbour of this block. Note that 
	 // the offset to this block is stored to index 13:
	 for (size_t i=0; i<nbrIDs.size(); ++i) {
	    if (cellOffsets[i] == INVALID_CELLID) {
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = INVALID_CELLID;
	    } else {
	       boundaryFlag = boundaryFlag | (1 << i);
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = cellOffsets[i] + block*SIZE_VELBLOCK;
	    }
	 }
	 // Boundary flags are stored to the last position in nbrsSpa array:
	 nbrsSpa[block*SIZE_NBRS_SPA + 30] = boundaryFlag;
      }
   }

   return true;
}

void initialLoadBalance(dccrg<SpatialCell>& mpiGrid) {
   SpatialCell::base_address_identifier = 5;
   mpiGrid.balance_load();
}

/*!
Sets spatial neighbor pointers of given cell to given list.

Non-existing neighbor pointers are set to NULL.
*/
bool findNeighbours(
	std::vector<const SpatialCell*>& nbrPtr,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid,
	#else
	const dccrg<SpatialCell>& mpiGrid,
	#endif
	const CellID& CELLID
) {

   CellID nbrID;
   for (int i=0; i<6; ++i) nbrPtr[i] = NULL;
   nbrID = getNeighbourID(mpiGrid, CELLID, 2-1, 2  , 2  );
   if (nbrID != INVALID_CELLID) nbrPtr[0] = mpiGrid[nbrID];
   nbrID = getNeighbourID(mpiGrid, CELLID, 2+1, 2  , 2  );
   if (nbrID != INVALID_CELLID) nbrPtr[1] = mpiGrid[nbrID];
   nbrID = getNeighbourID(mpiGrid, CELLID, 2  , 2-1, 2  );
   if (nbrID != INVALID_CELLID) nbrPtr[2] = mpiGrid[nbrID];
   nbrID = getNeighbourID(mpiGrid, CELLID, 2  , 2+1, 2  );
   if (nbrID != INVALID_CELLID) nbrPtr[3] = mpiGrid[nbrID];
   nbrID = getNeighbourID(mpiGrid, CELLID, 2  , 2  , 2-1);
   if (nbrID != INVALID_CELLID) nbrPtr[4] = mpiGrid[nbrID];
   nbrID = getNeighbourID(mpiGrid, CELLID, 2  , 2  , 2+1);
   if (nbrID != INVALID_CELLID) nbrPtr[5] = mpiGrid[nbrID];
   
   return true;
}

void calculateVelocityMoments(dccrg<SpatialCell>& mpiGrid) {
   Main::cells = mpiGrid.get_cells();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (Main::cellPtr != NULL) cpu_calcVelocityMoments(*Main::cellPtr);
   }
}

void calculateSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t, Real& dt) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcSimParameters(mpiGrid, t, dt);
}

void calculateCellParameters(dccrg<SpatialCell>& mpiGrid,creal& t, uint64_t cell) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcCellParameters(mpiGrid[cell]->cpu_cellParams,t);
}

void calculateAcceleration(dccrg<SpatialCell>& mpiGrid) {
   profile::start("calcAcceleration");
   
   // Calculate acceleration for all cells (inner + boundary):
   Main::cells = mpiGrid.get_cells();
   for (size_t i=0; i<Main::cells.size(); ++i) {
      Main::cellPtr = mpiGrid[Main::cells[i]];
      if (Main::cellPtr != NULL) cpu_acceleration(*Main::cellPtr);
   }
   
   profile::stop("calcAcceleration",Main::cells.size(),"SpatialCells");
}

void calculateSpatialDerivatives(dccrg<SpatialCell>& mpiGrid) {
   profile::start("calcSpatDerivatives");
   profile::start("Start data exchange");
   unsigned int computedCells;
   /* TODO: update N_blocks first if needed?
   SpatialCell::base_address_identifier = 5;
   mpiGrid.update_remote_neighbour_data();*/

   // Start neighbour data exchange:
   SpatialCell::base_address_identifier = 0;
   mpiGrid.start_remote_neighbour_data_update();
   profile::stop("Start data exchange");
   profile::start("Compute inner cells");
   // Calculate derivatives for inner cells:
   Main::cells = mpiGrid.get_cells_with_local_neighbours();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
   profile::stop("Compute inner cells",Main::cells.size(),"SpatialCells");
   computedCells=Main::cells.size();
   profile::start("Wait for receives");
   // Calculate derivatives for boundary cells when transfers have completed:
   mpiGrid.wait_neighbour_data_update_receives();
   profile::stop("Wait for receives");
   profile::start("Compute border cells");
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
   profile::stop("Compute border cells",Main::cells.size(),"SpatialCells");
   computedCells+=Main::cells.size();
   profile::start("Wait for sends");
   mpiGrid.wait_neighbour_data_update_sends();
   profile::stop("Wait for sends");
   profile::stop("calcSpatDerivatives",computedCells,"SpatialCells");
}

void calculateSpatialFluxes(dccrg<SpatialCell>& mpiGrid) {
   profile::start("calcSpatFluxes");
   profile::start("Start data exchange");
   unsigned int computedCells;
   /* TODO: update N_blocks first if needed?
   SpatialCell::base_address_identifier = 5;
   mpiGrid.update_remote_neighbour_data();*/

   // Start neighbour data exchange:
   SpatialCell::base_address_identifier = 1;
   mpiGrid.start_remote_neighbour_data_update();
   profile::stop("Start data exchange");
   profile::start("Compute inner cells");
   // Calculate fluxes for inner cells:
   Main::cells = mpiGrid.get_cells_with_local_neighbours();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   // Calculate fluxes for boundary cells when transfers have completed:
   profile::stop("Compute inner cells",Main::cells.size(),"SpatialCells");
   computedCells=Main::cells.size();
   profile::start("Wait for receives");
   mpiGrid.wait_neighbour_data_update_receives();
   profile::stop("Wait for receives");
   profile::start("Compute border cells"); 
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   profile::stop("Compute border cells",Main::cells.size(),"SpatialCells");
   computedCells+=Main::cells.size();
   profile::start("Wait for sends");
   mpiGrid.wait_neighbour_data_update_sends();
   profile::stop("Wait for sends");
   profile::stop("calcSpatFluxes",computedCells,"SpatialCells");
}

void calculateSpatialPropagation(dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) {
   profile::start("calcSpatProp");
   profile::start("Start data exchange");
   unsigned int computedCells;
   /* TODO: update N_blocks first if needed?
   SpatialCell::base_address_identifier = 5;
   mpiGrid.update_remote_neighbour_data();*/

   // Start neighbour data exchange:
   SpatialCell::base_address_identifier = 2;
   mpiGrid.start_remote_neighbour_data_update();
   profile::stop("Start data exchange");
   profile::start("Compute inner cells");
   // Propagate inner cells:
   Main::cells = mpiGrid.get_cells_with_local_neighbours();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   profile::stop("Compute inner cells",Main::cells.size(),"SpatialCells");
   computedCells=Main::cells.size();
   profile::start("Wait for receives");
   // Propagate boundary cells when transfers have completed:
   mpiGrid.wait_neighbour_data_update_receives();
   profile::stop("Wait for receives");
   profile::start("Compute border cells");
   
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   profile::stop("Compute border cells",Main::cells.size(),"SpatialCells");
   computedCells+=Main::cells.size();
   profile::start("Wait for sends");
   mpiGrid.wait_neighbour_data_update_sends();
   profile::stop("Wait for sends");
   profile::stop("calcSpatProp",computedCells,"SpatialCells");
}



#endif // #ifndef PARGRID

