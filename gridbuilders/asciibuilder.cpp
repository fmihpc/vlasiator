#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "asciibuilder.h"
#include "asciireaders.h"
#include "../parameters.h"
#include "../mpiconversion.h"

using namespace std;

static fstream filein;

AsciiBuilder::AsciiBuilder(): MPIBuilder(),initialized(false) { }

AsciiBuilder::~AsciiBuilder() {filein.close();}

bool AsciiBuilder::finalize() {filein.close(); return true;}

bool AsciiBuilder::getParameter(const std::string& parameterName,std::string& value) {
   return false;
}

bool AsciiBuilder::getTotalNumberOfCells(VirtualCell::ID& N_cells) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   N_cells = this->N_cells;
   return initialized;
}

bool AsciiBuilder::initialize(MPI_Comm comm,const int& MASTER_RANK) {
   initialized = MPIBuilder::initialize(comm,MASTER_RANK);

   if (mpiRank == mpiMasterRank) {
      vx_blocks = 10;
      vy_blocks = 10;
      vz_blocks = 10;
   
      // Open input file. Get the filename from somewhere.
      filein.open("vlasov.grid");
      if (filein.good() == false) return false;
      initialized = true;
      
      // Parse grid dimensionality:
      vector<string> line;
      if (readline(filein,line) == false) {initialized = false; return false;}
      if (line.size() != 2) {initialized = false; return false;}
      if (line[0] != "DIMENSIONS") {initialized = false; return false;}
      N_dimensions = atoi(line[1].c_str());
      line.clear();

      // Parse the number of cells:
      if (readline(filein,line) == false) {initialized = false; return false;}
      if (line.size() != 2) {initialized = false; return false;}
      if (line[0] != "CELLS") {initialized = false; return false;}
      N_cells = atoi(line[1].c_str());
      
      filein.close();
   }
   
   // Master broadcasts some globally needed parameter values:
   if (MPI_Bcast(&vx_blocks,1,MPI_Type<uint>(),MASTER_RANK,comm) != MPI_SUCCESS) initialized = false;
   if (MPI_Bcast(&vy_blocks,1,MPI_Type<uint>(),MASTER_RANK,comm) != MPI_SUCCESS) initialized = false;
   if (MPI_Bcast(&vz_blocks,1,MPI_Type<uint>(),MASTER_RANK,comm) != MPI_SUCCESS) initialized = false;
   
   // Set some parameter values (temp solution):
   Parameters::vxblocks_ini = vx_blocks;
   Parameters::vyblocks_ini = vy_blocks;
   Parameters::vzblocks_ini = vz_blocks;
   
   return initialized;
}

bool AsciiBuilder::getCellBlockData(const VirtualCell::ID& cellID,cuint& N_blocks,Real* blocks,Real* blockParams,uint* nbrsVel) {
   bool success = true;
   for (uint b=0; b<N_blocks; ++b) {
      for (uint j=0; j<SIZE_VELBLOCK; ++j)    blocks[b*SIZE_VELBLOCK+j]         = 0.0;
      for (uint j=0; j<SIZE_BLOCKPARAMS; ++j) blockParams[b*SIZE_BLOCKPARAMS+j] = 0.0;
      for (uint j=0; j<SIZE_NBRS_VEL; ++j)    nbrsVel[b*SIZE_NBRS_VEL+j]        = 0.0;
   }
   return success;
}

bool AsciiBuilder::getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs) {
   bool success = true;
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   
   cellIDs.resize(N_cells);
   N_nbrs.resize(N_cells);
   
   vector<string> line;
   filein.open("vlasov.grid");
   if (filein.good() == false) return false;
   if (readline(filein,line) == false) return false; // DIMENSIONS
   if (readline(filein,line) == false) return false; // CELLS
   
   uint counter = 0;
   line.clear();
   while (readline(filein,line) == true) {
      if (line.size() < 1 + 2*N_dimensions) {filein.close(); return false;}
      
      cellIDs[counter] = atoi(line[0].c_str());
      N_nbrs[counter] = line.size() - 1 - 2*N_dimensions;
      
      ++counter;
      line.clear();
   }   
   filein.close();
   return success;
}

bool AsciiBuilder::getCellNumberOfBlocks(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,uint* N_blocks) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false; // Only master rank has correct class parameters
   for (VirtualCell::ID i=0; i<N_cells; ++i) 
     N_blocks[i] = vx_blocks*vy_blocks*vz_blocks;
   return true;
}

bool AsciiBuilder::getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* nbrTypes) {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   bool success = true;

   VirtualCell::ID coordsCounter = 0;
   VirtualCell::ID nbrCounter = 0;
   
   for (VirtualCell::ID i=0; i<N_cells; ++i) {      
      vector<string> line;
      filein.open("vlasov.grid");
      if (filein.good() == false) {filein.close(); return false;}
      if (readline(filein,line) == false) {filein.close(); return false;} // DIMENSIONS
      if (readline(filein,line) == false) {filein.close(); return false;} // CELLS
      
      line.clear();
      bool found = false;
      while (readline(filein,line) == true) {
	 VirtualCell::ID cellIDread = atoi(line[0].c_str());
	 if (cellIDread != cellIDs[i]) {
	    line.clear();
	    continue; 
	 }
	 cuint N_nbrs = line.size() - 1 - 2*N_dimensions;
	 
	 unsigned int offset = 1;
	 for (int k=0; k<N_dimensions; ++k) {coords[coordsCounter+k] = atof(line[k+offset].c_str());}
	 offset += N_dimensions;
	 coordsCounter += N_dimensions;
	 for (int k=0; k<N_dimensions; ++k) {coords[coordsCounter+k] = atof(line[k+offset].c_str());}
	 offset += N_dimensions;
	 coordsCounter += N_dimensions;
	 
	 while (offset < line.size()) {
	    spatNbrIDs[nbrCounter  ] = atoi(line[offset].substr(0,line[offset].find(':')).c_str());
	    nbrTypes[nbrCounter]     = atoi(line[offset].substr(line[offset].find(':')+1,line[offset].size()-line[offset].find(':')-1).c_str());
	    ++offset;
	    ++nbrCounter;
	 }
	 found = true;
	 break;
      }      
      filein.close();
      if (found == false) success = false;
   }
   return success;
}

bool AsciiBuilder::getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams) {
   bool success = true;
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   
   for (VirtualCell::ID i=0; i<N_cells; ++i) {
      vector<string> line;
      filein.open("vlasov.grid");
      if (filein.good() == false) {filein.close(); return false;}
      if (readline(filein,line) == false) {filein.close(); return false;} // DIMENSIONS
      if (readline(filein,line) == false) {filein.close(); return false;} // CELLS
      
      line.clear();
      bool found = false;
      while (readline(filein,line) == true) {
	 VirtualCell::ID cellIDread = atoi(line[0].c_str());
	 if (cellIDread != cellIDs[i]) {
	    line.clear();
	    continue;
	 }
	 Real* tmp = new Real[2*N_dimensions];
	 for (uint k=0; k<2*N_dimensions; ++k) tmp[k] = atof(line[k+1].c_str());
	 
	 cellParams[i*SIZE_CELLPARAMS+CellParams::XCRD] = tmp[0];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::YCRD] = tmp[1];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::ZCRD] = tmp[2];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::DX  ] = tmp[3];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::DY  ] = tmp[4];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::DZ  ] = tmp[5];
	 
	 cellParams[i*SIZE_CELLPARAMS+CellParams::EX  ]  = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::EY  ]  = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::EZ  ]  = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::BX  ]  = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::BY  ]  = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::BZ  ]  = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::RHO  ] = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::RHOVX] = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::RHOVY] = 0.0;
	 cellParams[i*SIZE_CELLPARAMS+CellParams::RHOVZ] = 0.0;
	 found = true;
	 delete tmp;
	 tmp = NULL;
      }
      filein.close();
      if (found == false) success = false;
   }
   return success;
}

// Register AsciiBuilder:

GridBuilder* newAsciiBuilder() {return new AsciiBuilder;}

class Dummy {
 public:
   Dummy() {
      GridBuilderFactory::registerBuilder(newAsciiBuilder);
   }
   ~Dummy() { }
};

static Dummy dummy;
