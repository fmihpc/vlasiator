#ifndef ASCIIREADER_H
#define ASCIIREADER_H

#include <vector>
#include "mpibuilder.h"

class AsciiBuilder: public MPIBuilder {
 public:
   AsciiBuilder();
   ~AsciiBuilder();
   
   bool finalize();
   bool getCellBlockData(const VirtualCell::ID& cellID,cuint& N_blocks,Real* blocks,Real* blockParams,uint* nbrsVel);
   bool getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs);
   bool getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* nbrTypes);
   bool getCellNumberOfBlocks(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,uint* N_blocks);
   bool getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(VirtualCell::ID& N_cells);
   bool initialize(MPI_Comm comm,const int& MASTER_RANK);
   
 protected:
   bool initialized;
   uint N_cells;
   uint N_dimensions;
   uint vx_blocks;
   uint vy_blocks;
   uint vz_blocks;
};

#endif
