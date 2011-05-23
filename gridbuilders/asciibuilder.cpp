#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "asciibuilder.h"
#include "asciireaders.h"
#include "../parameters.h"
#include "../mpiconversion.h"
#include "../project.h"

using namespace std;

static fstream filein;

static string inputFile = "msphere.grid";

AsciiBuilder::AsciiBuilder(): MPIBuilder(),initialized(false) { }

AsciiBuilder::~AsciiBuilder() {filein.close();}

bool AsciiBuilder::calculatesAnalyticInitialState() {return false;}

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
   typedef Parameters P;
   if (mpiRank == mpiMasterRank) {
      //vx_blocks = 10;
      //vy_blocks = 10;
      //vz_blocks = 10;

      P::add("gridbuilder.x_min","Minimum value of the x-coordinate.",xmin,0.0);
      P::add("gridbuilder.x_max","Minimum value of the x-coordinate.",xmax,1.0);
      P::add("gridbuilder.y_min","Minimum value of the y-coordinate.",ymin,0.0);
      P::add("gridbuilder.y_max","Minimum value of the y-coordinate.",ymax,1.0);
      P::add("gridbuilder.z_min","Minimum value of the z-coordinate.",zmin,0.0);
      P::add("gridbuilder.z_max","Minimum value of the z-coordinate.",zmax,1.0);
      P::add("gridbuilder.x_length","Number of cells in x-direction in initial grid.",xsize,1);
      P::add("gridbuilder.y_length","Number of cells in y-direction in initial grid.",ysize,1);
      P::add("gridbuilder.z_length","Number of cells in z-direction in initial grid.",zsize,1);
      P::add("gridbuilder.vx_min","Minimum value for velocity block vx-coordinates.",vx_min,-numeric_limits<Real>::max());
      P::add("gridbuilder.vx_max","Maximum value for velocity block vx-coordinates.",vx_max,+numeric_limits<Real>::max());
      P::add("gridbuilder.vy_min","Minimum value for velocity block vy-coordinates.",vy_min,-numeric_limits<Real>::max());
      P::add("gridbuilder.vy_max","Maximum value for velocity block vy-coordinates.",vy_max,+numeric_limits<Real>::max());
      P::add("gridbuilder.vz_min","Minimum value for velocity block vz-coordinates.",vz_min,-numeric_limits<Real>::max());
      P::add("gridbuilder.vz_max","Maximum value for velocity block vz-coordinates.",vz_max,+numeric_limits<Real>::max());
      P::add("gridbuilder.vx_length","Number of cells in vx-direction in initial grid.",vx_blocks,1);
      P::add("gridbuilder.vy_length","Number of cells in vy-direction in initial grid.",vy_blocks,1);
      P::add("gridbuilder.vz_length","Number of cells in vz-direction in initial grid.",vz_blocks,1);
      P::parse();
      dx = (xmax-xmin)/xsize;
      dy = (ymax-ymin)/ysize;
      dz = (zmax-zmin)/zsize;
      
      // Open input file. Get the filename from somewhere.
      //filein.open("vlasov.grid");
      filein.open(inputFile.c_str());
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
   return true;
   /*
   for (uint b=0; b<N_blocks; ++b) {
      for (uint j=0; j<SIZE_VELBLOCK; ++j)    blocks[b*SIZE_VELBLOCK+j]         = 0.0;
      for (uint j=0; j<SIZE_BLOCKPARAMS; ++j) blockParams[b*SIZE_BLOCKPARAMS+j] = 0.0;
      for (uint j=0; j<SIZE_NBRS_VEL; ++j)    nbrsVel[b*SIZE_NBRS_VEL+j]        = 0.0;
   }*/
   namespace VC = VirtualCell;
   const VC::ID K = cellID / (ysize*xsize);
   const VC::ID J = (cellID - K*ysize*xsize)/xsize;
   const VC::ID I = cellID - K*ysize*xsize - J*xsize;
   
   creal dvx_block = (vx_max-vx_min)/vx_blocks; // Size of velocity block in  vx
   creal dvy_block = (vy_max-vy_min)/vy_blocks; //                            vy
   creal dvz_block = (vz_max-vz_min)/vz_blocks; //                            vz
   creal dvx_blockCell = dvx_block / WID;       // Size of a cell in block in vx
   creal dvy_blockCell = dvy_block / WID;       //                            vy
   creal dvz_blockCell = dvz_block / WID;       //                            vz
   for (uint kv=0; kv<vz_blocks; ++kv) for (uint jv=0; jv<vy_blocks; ++jv) for (uint iv=0; iv<vx_blocks; ++iv) {
      cuint BLOCKID = kv*vy_blocks*vx_blocks + jv*vx_blocks + iv;
      // Calculate velocity block parameters:
      creal vx_block = vx_min + iv*dvx_block;
      creal vy_block = vy_min + jv*dvy_block;
      creal vz_block = vz_min + kv*dvz_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::VXCRD] = vx_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::VYCRD] = vy_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::VZCRD] = vz_block;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::DVX  ] = dvx_blockCell;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::DVY  ] = dvy_blockCell;
      blockParams[BLOCKID*SIZE_BLOCKPARAMS + BlockParams::DVZ  ] = dvz_blockCell;
      
      for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	 creal vx_cell = vx_block + ic*dvx_blockCell;
	 creal vy_cell = vy_block + jc*dvy_blockCell;
	 creal vz_cell = vz_block + kc*dvz_blockCell;
	 blocks[BLOCKID*SIZE_VELBLOCK + kc*WID2+jc*WID+ic] =
	   calcPhaseSpaceDensity(xmin+I*dx,ymin+J*dy,zmin+K*dz,dx,dy,dz,vx_cell,
				 vy_cell,vz_cell,dvx_blockCell,dvy_blockCell,dvz_blockCell);
      }
      cuint vxneg_nbr = iv-1;
      cuint vxpos_nbr = iv+1;
      cuint vyneg_nbr = jv-1;
      cuint vypos_nbr = jv+1;
      cuint vzneg_nbr = kv-1;
      cuint vzpos_nbr = kv+1;
      uint state = 0;
      if (iv == 0)           state = state | NbrsVel::VX_NEG_BND; // Check for boundaries
      if (iv == vx_blocks-1) state = state | NbrsVel::VX_POS_BND;
      if (jv == 0)           state = state | NbrsVel::VY_NEG_BND;
      if (jv == vy_blocks-1) state = state | NbrsVel::VY_POS_BND;
      if (kv == 0)           state = state | NbrsVel::VZ_NEG_BND;
      if (kv == vz_blocks-1) state = state | NbrsVel::VZ_POS_BND;
      
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::STATE] = state;
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::MYIND] = BLOCKID;
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::VXNEG] = velBlockIndex(vxneg_nbr,jv,kv);
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::VXPOS] = velBlockIndex(vxpos_nbr,jv,kv);
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::VYNEG] = velBlockIndex(iv,vyneg_nbr,kv);
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::VYPOS] = velBlockIndex(iv,vypos_nbr,kv);
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::VZNEG] = velBlockIndex(iv,jv,vzneg_nbr);
      nbrsVel[BLOCKID*SIZE_NBRS_VEL + NbrsVel::VZPOS] = velBlockIndex(iv,jv,vzpos_nbr);
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
   filein.open(inputFile.c_str());
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
   cerr << "getCellNbrData starting" << endl;
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   bool success = true;
   
   VirtualCell::ID coordsCounter = 0;
   VirtualCell::ID nbrCounter = 0;

   // Read the whole file into a vector:
   vector<string> line;
   map<VirtualCell::ID,vector<string> > inFile;
   filein.open(inputFile.c_str());
   if (filein.good() == false) {filein.close(); return false;}
   if (readline(filein,line) == false) {filein.close(); return false;} // DIMENSIONS
   if (readline(filein,line) == false) {filein.close(); return false;} // CELLS
   line.clear();
   while(readline(filein,line) == true) {
      VirtualCell::ID cellIDread = atoi(line[0].c_str());
      inFile[cellIDread] = line;
      line.clear();
   }
   filein.close();
   
   // Process requests:
   map<VirtualCell::ID,vector<string> >::iterator it;
   for (VirtualCell::ID i=0; i<N_cells; ++i) {
      it = inFile.find(cellIDs[i]);
      if (it == inFile.end()) return false;
      
      unsigned int offset = 1;
      for (int k=0; k<N_dimensions; ++k) {coords[coordsCounter+k] = atof(it->second[k+offset].c_str());}
      offset += N_dimensions;
      coordsCounter += N_dimensions;
      for (int k=0; k<N_dimensions; ++k) {coords[coordsCounter+k] = atof(it->second[k+offset].c_str());}
      offset += N_dimensions;
      coordsCounter += N_dimensions;
      offset += 6; // E,B fields

      while (offset < it->second.size()) {
	 spatNbrIDs[nbrCounter  ] = atoi(it->second[offset].substr(0,it->second[offset].find(':')).c_str());
	 nbrTypes[nbrCounter]     = atoi(it->second[offset].substr(it->second[offset].find(':')+1,it->second[offset].size()-it->second[offset].find(':')-1).c_str());
	 ++offset;
	 ++nbrCounter;
      }
   }
   return success;
   /*
   for (VirtualCell::ID i=0; i<N_cells; ++i) {      
      vector<string> line;
      filein.open(inputFile.c_str());
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

	 offset += 6; // E,B fields
	 
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
    */
}

bool AsciiBuilder::getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams) {
   cerr << "getCellParams starting" << endl;
   bool success = true;
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return false;
   
   // Read the whole file into a map for fast searching:
   vector<string> line;
   map<VirtualCell::ID,vector<string> > inFile;
   filein.open(inputFile.c_str());
   if (filein.good() == false) {filein.close(); return false;}
   if (readline(filein,line) == false) {filein.close(); return false;} // DIMENSIONS
   if (readline(filein,line) == false) {filein.close(); return false;} // CELLS
   line.clear();
   while(readline(filein,line) == true) {
      VirtualCell::ID cellIDread = atoi(line[0].c_str());
      inFile[cellIDread] = line;
      line.clear();
   }
   filein.close();
   
   map<VirtualCell::ID,vector<string> >::iterator it;
   for (VirtualCell::ID i=0; i<N_cells; ++i) {
      it = inFile.find(cellIDs[i]);
      if (it == inFile.end()) return false;
      
      cellParams[i*SIZE_CELLPARAMS+CellParams::XCRD] = atof(it->second[1 + 0].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::YCRD] = atof(it->second[1 + 1].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::ZCRD] = atof(it->second[1 + 2].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::DX  ] = atof(it->second[1 + 3].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::DY  ] = atof(it->second[1 + 4].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::DZ  ] = atof(it->second[1 + 5].c_str());
      
      cellParams[i*SIZE_CELLPARAMS+CellParams::EX  ]  = atof(it->second[2*N_dimensions+1 + 0].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::EY  ]  = atof(it->second[2*N_dimensions+1 + 1].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::EZ  ]  = atof(it->second[2*N_dimensions+1 + 2].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::BX  ]  = atof(it->second[2*N_dimensions+1 + 3].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::BY  ]  = atof(it->second[2*N_dimensions+1 + 4].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::BZ  ]  = atof(it->second[2*N_dimensions+1 + 5].c_str());
      cellParams[i*SIZE_CELLPARAMS+CellParams::RHO  ] = 0.0;
      cellParams[i*SIZE_CELLPARAMS+CellParams::RHOVX] = 0.0;
      cellParams[i*SIZE_CELLPARAMS+CellParams::RHOVY] = 0.0;
      cellParams[i*SIZE_CELLPARAMS+CellParams::RHOVZ] = 0.0;
   }
   cerr << "\t Finished" << endl;
   return success;
   /*
   for (VirtualCell::ID i=0; i<N_cells; ++i) {
      vector<string> line;
      filein.open(inputFile.c_str());
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
	 
	 for (uint k=0; k<6; ++k) tmp[k] = atof(line[1+2*N_dimensions+k].c_str());
	 
	 cellParams[i*SIZE_CELLPARAMS+CellParams::EX  ]  = tmp[0];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::EY  ]  = tmp[1];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::EZ  ]  = tmp[2];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::BX  ]  = tmp[3];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::BY  ]  = tmp[4];
	 cellParams[i*SIZE_CELLPARAMS+CellParams::BZ  ]  = tmp[5];
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
   */
}

uint AsciiBuilder::velBlockIndex(cuint& iv,cuint& jv,cuint& kv) {
   return kv*vy_blocks*vx_blocks + jv*vx_blocks + iv;
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
