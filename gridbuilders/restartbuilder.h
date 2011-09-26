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

#ifndef RESTARTBUILDER_H
#define RESTARTBUILDER_H

#include <map>
#include "../gridbuilder.h"
#include "../vlsvreader2.h"

class RestartBuilder: public GridBuilder {
 public:
   RestartBuilder();
   virtual ~RestartBuilder();
   
   // **************************************************
   // ***** VIRTUAL FUNCTIONS OF CLASS GRIDBUILDER *****
   // *****       DEFINED BY RESTARTBUILDER        *****
   // **************************************************
   
   bool addCellBlockDataRequests(VirtualCell::ID& totalCells,VirtualCell::ID& blockOffset,VirtualCell::ID* cellIDs,uint* blocksPerCell,
				 Real** avgsBuffer,Real** blockParamsBuffer,uint** nbrsVelBuffer);
   bool addCellBlockNumberRequests(VirtualCell::ID& totalCells,VirtualCell::ID& cellOffset,VirtualCell::ID* cellIDs,uint* N_blocks);
   bool addCellNbrRequests(VirtualCell::ID& totalCells,VirtualCell::ID& totalNbrs,VirtualCell::ID& cellOffset,
			   VirtualCell::ID& nbrOffset,VirtualCell::ID* cellIDs,
			   uchar* nbrsPerCell,Real* coords,VirtualCell::ID* nbrIDs,uchar* nbrTypes);
   bool addCellParamsRequests(VirtualCell::ID& totalCells,VirtualCell::ID& cellOffset,VirtualCell::ID* cellIDs,Real* cellParams);
   bool doInitialLoadBalance();
   bool finalize();
   bool getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(VirtualCell::ID& N_cells);
   bool initialize(MPI_Comm comm,const int& masterRank);
   bool processCellBlockDataRequests();
   bool processCellBlockNumberRequests();
   bool processCellNbrRequests();
   bool processCellParamsRequests();
   bool waitCellBlockDataRequests();
   bool waitCellBlockNumberRequests();
   bool waitCellNbrRequests();
   bool waitCellParamsRequests();
   
 private:
   VirtualCell::ID* cellIDs;
   uchar* cellNbrs;
   MPI_Comm comm;                             /**< MPI communicator used to read the restart file in parallel.*/
   std::string fileName;                      /**< Name of the restart file.*/
   bool initialized;                          /**< If true, RestartBuilder has initialized correctly.*/
   int masterRank;                            /**< MPI rank of the master process.*/
   int myRank;                                /**< MPI rank of this process.*/
   VirtualCell::ID N_cells;                   /**< Number of spatial cells in restart file.*/
   std::map<std::string,std::string> options; /**< Option,value pairs which can be queried from RestartBuilder.*/
   int processes;                             /**< Number of MPI processes in communicator comm.*/
   VLSVParReader vlsvReader;                  /**< VLSV parallel file reader.*/
};

#endif
