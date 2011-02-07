#ifdef PARGRID
#ifndef MAIN_H
#define MAIN_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "mpilogger.h"
#include "pargrid.h"
#include "definitions.h"
#include "cell_spatial.h"
#include "project.h"
#include "timer.h"
#include "vlsreader.h"

// NOTE: If preprocessor flag PROFILE is undefined, the compiler should optimize out Timer:: calls.

MPILogger mpilogger;

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);

namespace Main {
   std::vector<ID::type> cells;
   std::vector<const SpatialCell*> nbrPtrs(6,NULL);
   SpatialCell* cellPtr;
   
   cuint calcAcc           = Timer::create("Computing: vel. propagation  (total) : ");
   cuint calcSpatDerivs    = Timer::create("Computing: spat. derivatives (total) : ");
   cuint spatDerivsMPIRecv = Timer::create("MPI Recv : spat. derivs              : ");
   cuint spatDerivsMPISend = Timer::create("MPI Send : spat. derivs              : ");
   cuint calcSpatFluxes    = Timer::create("Computing: spat. fluxes      (total) : ");
   cuint spatFluxesMPIRecv = Timer::create("MPI Recv : spat. fluxes              : ");
   cuint spatFluxesMPISend = Timer::create("MPI Send : spat. fluxes              : ");
   cuint calcSpatProp      = Timer::create("Computing: spat. propag      (total) : ");
   cuint spatPropMPIRecv   = Timer::create("MPI Recv : spat. propag              : ");
   cuint spatPropMPISend   = Timer::create("MPI Send : spat. propag              : ");
}

void initialLoadBalance(ParGrid<SpatialCell>& mpiGrid) { 
   /*
   int myrank;
   int N_processes;
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   for (int i=0; i<N_processes; ++i) {
      if (i == myrank) {
	 mpiGrid.getCells(Main::cells);
	 std::cerr << "Proc #" << myrank << " has cells:" << std::endl;
	 for (uint j=0; j<Main::cells.size(); ++j) std::cerr << "\t" << Main::cells[j] << std::endl;
      }
      mpiGrid.barrier();
   }
    */
}
/*
bool initFromRestartFile(const std::string& fname,const int& myrank,ParGrid<SpatialCell>& mpiGrid) {
   const int MASTER_RANK = 0;
   bool success = true;
   bool masterSuccess = false;
   VlsReader vlsReader;
   long long unsigned int N_cells = 0;
   int N_processes;
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   
   // Master process counts the number of cells in restart file:
   if (myrank == MASTER_RANK) {
      if (vlsReader.open(fname) == false) {success = false;}
      if (success == true && vlsReader.readHeader() == false) {success = false;}
      
      // Do some sanity checks based on file header:
      unsigned char bytesPerCRD;
      unsigned char bytesPerGID;
      unsigned char hasNeighbourList;
      if (vlsReader.getHeaderElement(VlsHeader::BYTES_PER_CELL_CRD,bytesPerCRD) == false) success = false;
      if (vlsReader.getHeaderElement(VlsHeader::BYTES_PER_CELL_GID,bytesPerGID) == false) success = false;
      if (vlsReader.getHeaderElement(VlsHeader::CONTAINS_SPAT_NBR_LIST,hasNeighbourList) == false) success = false;
      
      if (sizeof(ID::type) < bytesPerGID) {
	 std::cerr << "(Main) Restart ERROR: File has wider global IDs than local computer!" << std::endl;
	 std::cerr << "\tFile has " << bytesPerGID << " byte IDs, this computer has " << sizeof(ID::type) << " byte wide IDs" << std::endl;
	 success = false;
      }
      if (sizeof(Real) < bytesPerCRD) {
	 std::cerr << "(Main) Restart WARNING: File has wider floating points than local computer!" << std::endl;
      }
      if (hasNeighbourList == 0) {
	 std::cerr << "(Main) Restart ERROR: File does not contain neighbour lists!" << std::endl;
	 success = false;
      }
      
      // Count cells:
      if (success == true) while (vlsReader.readSpatCellCoordEntry() == true) ++N_cells;
      vlsReader.close();
   }
   mpiGrid.barrier();
   
   // Master process has to let everyone know if cell count succeeded:
   if (MPI_Scatter(&success,1,MPI_CHAR,&masterSuccess,1,MPI_CHAR,MASTER_RANK,MPI_COMM_WORLD) != MPI_SUCCESS) {
      std::cerr << "Failed to scatter" << std::endl;
   }
   if (masterSuccess == false) {
      if (myrank == MASTER_RANK) std::cerr << "(MAIN) Restart: Master failed to count cells!" << std::endl;
      return false;
   }
   
   // Create an MPI Datatype for sendinf cell data from file:
   struct CellData {
      int receiver;
      ID::type cellID;
      Real coords[6];
      unsigned char refLevel;
      ID::type nbrs[24];
   } cellData;
   
   MPI_Datatype types[5];
   types[0] = MPI_INT;                   // Rank of receiving process
   types[1] = MPI_Type<ID::type>();      // Cell ID
   types[2] = MPI_Type<Real>();          // xcrd,ycrd,zcrd,dx,dy,dz
   types[3] = MPI_Type<unsigned char>(); // Cell refinement level
   types[4] = MPI_Type<ID::type>();      // Neighbour IDs (can have max. 24 nbrs)
   
   int blocklen[] = {1,1,6,1,24};
   
   // There is something weird going on with MPI_Get_address, so let's calculate
   // the offsets manually (prohibited to use void pointers in arithmetics):
   MPI_Aint disp[5];
   disp[0] = 0;
   disp[1] = reinterpret_cast<char*>(&cellData.cellID) - reinterpret_cast<char*>(&cellData.receiver);
   disp[2] = reinterpret_cast<char*>(cellData.coords) - reinterpret_cast<char*>(&cellData.receiver);
   disp[3] = reinterpret_cast<char*>(&cellData.refLevel) - reinterpret_cast<char*>(&cellData.receiver);
   disp[4] = reinterpret_cast<char*>(cellData.nbrs) - reinterpret_cast<char*>(&cellData.receiver);
   
   MPI_Datatype mpiDataType;
   MPI_Type_create_struct(5,blocklen,disp,types,&mpiDataType);
   MPI_Type_commit(&mpiDataType);
   MPI_Status mpiStatus;
   
   // Master process starts handing out cells:
   if (myrank == MASTER_RANK) {
      if (vlsReader.open(fname) == false) success = false;
      if (success == true && vlsReader.readHeader() == false) success = false;
      unsigned long long int cellsPassed = 0;
      for (int proc=0; proc<N_processes; ++proc) {
	 // Count the number of cells process proc receives:
	 unsigned long long int cellsToProc = N_cells/N_processes;
	 if (proc < N_cells%N_processes) ++cellsToProc;
	 
	 // Send cellsToProc cells to process proc:
	 unsigned long long int cellCounter = 0;
	 while (cellCounter < cellsToProc) {
	    // Read cell data from file:
	    if (success == true && vlsReader.readSpatCellCoordEntry() == false) success = false;
	    
	    // Check that everything is OK.
	    if (success == false) break;
	    
	    cellData.receiver  = proc;
	    cellData.cellID    = vlsReader.getSpatCellGID<ID::type>();
	    cellData.coords[0] = vlsReader.getCrdX<Real>();
	    cellData.coords[1] = vlsReader.getCrdY<Real>();
	    cellData.coords[2] = vlsReader.getCrdZ<Real>();
	    cellData.coords[3] = vlsReader.getDx<Real>();
	    cellData.coords[4] = vlsReader.getDy<Real>();
	    cellData.coords[5] = vlsReader.getDz<Real>();
	    cellData.refLevel  = vlsReader.getRefinementLevel();
	    
	    for (int i=0; i<24; ++i) cellData.nbrs[i] = std::numeric_limits<ID::type>::max(); // Important!
	    for (int i=0; i<vlsReader.getNumberOfSpatialNbrs(); ++i)
	      cellData.nbrs[i] = vlsReader.getNeighbourID(i);
	    
	    if (myrank == proc) {
	       // Cell is for master process:
	       mpiGrid.addCell(cellData.cellID,cellData.coords[0],cellData.coords[1],cellData.coords[2],
			       cellData.coords[3],cellData.coords[4],cellData.coords[5],
			       cellData.refLevel,cellData.nbrs);
	    } else {
	       // Send cell to process proc:
   	       if (MPI_Send(&cellData.receiver,1,mpiDataType,proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
		  mpilogger << "(MAIN) Restart: error occurred while sending cells!" << std::endl << write;
	       }
	    }
	    ++cellCounter;
	    ++cellsPassed;
	 }
	 if (success == false) break;
      }
      // All cells have been passed. Send invalid cellID to every process:
      cellData.cellID = std::numeric_limits<ID::type>::max();
      for (int proc=0; proc<N_processes; ++proc) {
	 if (proc == myrank) continue;
	 if (MPI_Send(&cellData.receiver,1,mpiDataType,proc,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
	    mpilogger << "(MAIN) Restart: error occurred while sending exit signal!" << std::endl << write;
	 }
      }
      if (cellsPassed != N_cells) {
	 std::cerr << "(MAIN) Restart ERROR: master did not pass all cells!" << std::endl;
	 success = false;
      }
      vlsReader.close();
   } else {
      // Keep receiving cells until master process sends an invalid cell ID:
      if (MPI_Recv(&cellData.receiver,1,mpiDataType,MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	 mpilogger << "(MAIN) Restart: error occurred while receiving cells!" << std::endl << write;
      }
      while (cellData.cellID != std::numeric_limits<ID::type>::max()) {
	 // Add cell to ParGrid:
	 mpiGrid.addCell(cellData.cellID,cellData.coords[0],cellData.coords[1],cellData.coords[2],
			 cellData.coords[3],cellData.coords[4],cellData.coords[5],
			 cellData.refLevel,cellData.nbrs);
	 
	 if (MPI_Recv(&cellData.receiver,1,mpiDataType,MASTER_RANK,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus) != MPI_SUCCESS) {
	    mpilogger << "(MAIN) Restart: error occurred while receiving cells!" << std::endl << write;
	 }
      }
   }
   MPI_Type_free(&mpiDataType);
   mpiGrid.barrier();

   // Master process lets everyone know if cell data was passed successfully:
   if (MPI_Scatter(&success,1,MPI_CHAR,&masterSuccess,1,MPI_CHAR,MASTER_RANK,MPI_COMM_WORLD) != MPI_SUCCESS)
     std::cerr << "(Main) Restart ERROR: Master failed to scatter!" << std::endl;
   if (masterSuccess == false) {
      mpilogger << "(MAIN) Restart ERROR: Master send me success=false!" << std::endl << write;
      return false;
   }
   
   // Now we have distributed cell IDs, cell coordinates, and existing cell 
   // neighbour IDs to ParGrid. Init ParGrid:
   if (mpiGrid.initialize() == false) success = false;
   
   // As a final thing remaining, each process needs now to open the 
   // restart file and read its data:
   if (vlsReader.open(fname) == false) {success = false;}
   if (success == true && vlsReader.readHeader() == false) {success = false;}
   if (success == true) {
      while (vlsReader.readSpatCellCoordEntry() == true) {
	 // The following code also reads remote cell data, but 
	 // it is probably not worth correcting (mpiGrid[cellID] returns a 
	 // valid pointer to localCells or remoteCells).
	 const ID::type cellID = vlsReader.getSpatCellGID<ID::type>();
	 SpatialCell* cellptr = mpiGrid[cellID];
	 if (cellptr == NULL) continue;
	 
	 // Read distribution function data:
	 
      }
   }
   
   vlsReader.close();
   return success;
}
*/
bool findNeighbours(std::vector<const SpatialCell*>& nbrPtr,const ParGrid<SpatialCell>& mpiGrid,const ID::type& CELLID) {
   /*
   std::vector<ID::type> nbrs;
   mpiGrid.getNeighbours(nbrs,CELLID);
   if (nbrs.size() != 24) {
      mpilogger << "findNeighbours ERROR: Failed to get neighbour indices!" << std::endl << write;
      return false;
   }
   // Get a pointer to each neighbour. If the neighbour does not exists, 
   // a NULL pointer is inserted:
   if (nbrs[ 0] != std::numeric_limits<ID::type>::max()) nbrPtr[0] = mpiGrid[nbrs[ 0]];
   else nbrPtr[0] = NULL;
   if (nbrs[ 4] != std::numeric_limits<ID::type>::max()) nbrPtr[1] = mpiGrid[nbrs[ 4]];
   else nbrPtr[1] = NULL;
   if (nbrs[ 8] != std::numeric_limits<ID::type>::max()) nbrPtr[2] = mpiGrid[nbrs[ 8]];
   else nbrPtr[2] = NULL;
   if (nbrs[12] != std::numeric_limits<ID::type>::max()) nbrPtr[3] = mpiGrid[nbrs[12]];
   else nbrPtr[3] = NULL;
   if (nbrs[16] != std::numeric_limits<ID::type>::max()) nbrPtr[4] = mpiGrid[nbrs[16]];
   else nbrPtr[4] = NULL;
   if (nbrs[20] != std::numeric_limits<ID::type>::max()) nbrPtr[5] = mpiGrid[nbrs[20]];
   else nbrPtr[5] = NULL;
   */
   ID::type nbrID;
   for (int i=0; i<6; ++i) nbrPtr[i] = NULL;
   nbrID = mpiGrid.getNeighbour(CELLID, 0);
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[0] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID, 8);
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[1] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,16);
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[2] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,24);
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[3] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,32);
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[4] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,40);
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[5] = mpiGrid[nbrID];
   
   return true;
}

void calculateSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, Real& dt) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcSimParameters(mpiGrid, t, dt);
}

void calculateCellParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, ID::type cell) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcCellParameters(mpiGrid[cell]->cpu_cellParams,t);
}

void calculateAcceleration(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcAcc);

   // Calculate acceleration for all cells (inner + boundary):
   mpiGrid.getCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (Main::cellPtr != NULL) cpu_acceleration(*Main::cellPtr);
   }

   Timer::stop(Main::calcAcc);
}

void calculateSpatialDerivatives(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatDerivs);

   // Start neighbour data exchange:
   mpiGrid.startNeighbourExchange(0);
   // Calculate derivatives for inner cells:
   mpiGrid.getInnerCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }

   #ifdef PARGRID_WAITANY
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitAnyReceive() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #elif PARGRID_WAITSOME
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitSomeReceives() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #else
      // Calculate derivatives for boundary cells when transfers have completed:
      Timer::start(Main::spatDerivsMPIRecv);
      mpiGrid.waitAllReceives();
      Timer::stop(Main::spatDerivsMPIRecv);
   
      mpiGrid.getBoundaryCells(Main::cells);
      for (size_t c=0; c<Main::cells.size(); ++c) {
	 Main::cellPtr = mpiGrid[Main::cells[c]];
	 if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	    mpilogger << "Failed to find neighbours." << std::endl << write;
	    continue;
	 }
	 if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
      }
   #endif
   // Wait for all sends to complete:
   Timer::start(Main::spatDerivsMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatDerivsMPISend);
   Timer::stop(Main::calcSpatDerivs);
}

void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatFluxes);

   // Start neighbour data exchange:
   mpiGrid.startNeighbourExchange(1);
   // Calculate fluxes for inner cells:
   mpiGrid.getInnerCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   
   #ifdef PARGRID_WAITSOME
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitSomeReceives() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #else
      // Calculate fluxes for boundary cells when transfers have completed:
      Timer::start(Main::spatFluxesMPIRecv);
      mpiGrid.waitAllReceives();
      Timer::stop(Main::spatFluxesMPIRecv);
   
      mpiGrid.getBoundaryCells(Main::cells);
      for (size_t c=0; c<Main::cells.size(); ++c) {
	 Main::cellPtr = mpiGrid[Main::cells[c]];
	 if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	    mpilogger << "Failed to find neighbours." << std::endl << write; 
	    continue;
	 }
	 if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
      }
   #endif
   Timer::start(Main::spatFluxesMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatFluxesMPISend);
   Timer::stop(Main::calcSpatFluxes);
}

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatProp);

   // Start neighbour data exchange:
   mpiGrid.startNeighbourExchange(2);
   // Start neighbour data exchange:
   mpiGrid.getInnerCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   #ifdef PARGRID_WAITSOME
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitSomeReceives() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #else
      // Propagate boundary cells when transfers have completed:
      Timer::start(Main::spatPropMPIRecv);
      mpiGrid.waitAllReceives();
      Timer::stop(Main::spatPropMPIRecv);
   
      mpiGrid.getBoundaryCells(Main::cells);
      for (size_t c=0; c<Main::cells.size(); ++c) {
	 Main::cellPtr = mpiGrid[Main::cells[c]];
	 if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	    mpilogger << "Failed to find neighbours." << std::endl << write;
	    continue;
	 }
	 if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
      }
   #endif
   Timer::start(Main::spatPropMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatPropMPISend);
   Timer::stop(Main::calcSpatProp);
}

void writeTimers() {Timer::print();}

#endif // #ifndef MAIN_H
#endif // #ifdef PARGRID
