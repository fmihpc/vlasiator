#include <cstdlib>
#include <iostream>

#include <silo.h>

#include "silowriter.h"
#include "common.h"
#include "cell_spatial.h"

using namespace std;

static DBfile* fileptr = NULL;

static int* nodeList = NULL;
static uint blockCntr = 0;
static real* x = NULL;
static real* y = NULL;
static real* z = NULL;
static real* avgs = NULL;

void allocateArrays(cuint& BLOCKS) {
   blockCntr = 0;
   nodeList = new int[BLOCKS*SIZE_VELBLOCK*8];
   x = new real[BLOCKS*SIZE_VELBLOCK*8];
   y = new real[BLOCKS*SIZE_VELBLOCK*8];
   z = new real[BLOCKS*SIZE_VELBLOCK*8];
   avgs = new real[BLOCKS*SIZE_VELBLOCK];
}

void allocateCellArrays(cuint& CELLS) {
   blockCntr = 0;
   nodeList = new int[CELLS*8];
   x = new real[CELLS*8];
   y = new real[CELLS*8];
   z = new real[CELLS*8];
   avgs = new real[CELLS];
}

void deallocateArrays() {
   delete nodeList;
   delete x;
   delete y;
   delete z;
   delete avgs;
   
   nodeList = NULL;
   x = NULL;
   y = NULL;
   z = NULL;
   avgs = NULL;
}

bool openOutputFile(const std::string& fileName,const std::string& fileInfo) {
   fileptr = DBCreate(fileName.c_str(),DB_CLOBBER,DB_LOCAL,fileInfo.c_str(),DB_PDB);
   if (fileptr == NULL) return false;
   else return true;
}

bool closeOutputFile() {
   if (DBClose(fileptr) < 0) return false;
   return true;
}

bool freeBlocks() {
   deallocateArrays();
   return true;
}

bool freeCells() {
   deallocateArrays();
   return true;
}

bool reserveSpatialCells(cuint& CELLS) {
   allocateCellArrays(CELLS);
   return true;
}

bool reserveVelocityBlocks(cuint& BLOCKS) {
   allocateArrays(BLOCKS);
   return true;
}
/*
bool writeBlocks(const std::string& gridName,const std::string& varName) {
   if (fileptr == NULL) return false;
   writeReservedBlocks(gridName);
   if (DBPutUcdvar1(fileptr,varName.c_str(),gridName.c_str(),avgs,blockCntr*SIZE_VELBLOCK,NULL,0,DB_FLOAT,DB_ZONECENT,NULL) < 0) return false;
   return true;
}
*/
bool writeReservedBlocks(const std::string& gridName) {
   if (fileptr == NULL) return false;
   cuint Ncrds = WID+1; // Each dimension has 5 coordinate values per block
   
   int Ndims = 3;
   int NshapeTypes = 1;
   
   int Nzones = blockCntr*SIZE_VELBLOCK;
   int shapeTypes[] = {DB_ZONETYPE_HEX}; // Hexahedrons only
   int shapeSizes[] = {8}; // Each block contains WID3 cells
   int shapeCnt[] = {Nzones};
   int Nshapes = 1;
   
   DBPutZonelist2(fileptr,"zonelist",Nzones,Ndims,nodeList,Nzones*8,0,0,0,shapeTypes,shapeSizes,shapeCnt,Nshapes,NULL);
   
   char* crdNames[3];
   crdNames[0] = const_cast<char*>("X-axis");
   crdNames[1] = const_cast<char*>("Y-axis");
   crdNames[2] = const_cast<char*>("Z-axis");
   real* crdArrays[] = {x,y,z};
   
   if (DBPutUcdmesh(fileptr,gridName.c_str(),Ndims,crdNames,crdArrays,Nzones*8,Nzones,"zonelist",NULL,DB_FLOAT,NULL) < 0) return false;
   return true;
}
/*
bool addVelocityBlock(real* blockParams,real* block) {
   if (fileptr == NULL) return false;

   for (uint i=0; i<WID3; ++i) avgs[blockCntr*SIZE_VELBLOCK + i] = block[i];
   addVelocityGridBlock3D(blockParams);
}*/

bool writeSpatialCells(const std::string& gridName,const std::string& varName) {
   if (fileptr == NULL) return false;
   
   int Ndims = 3;
   int NshapeTypes = 1;
   
   int Nzones = blockCntr;
   int shapeTypes[] = {DB_ZONETYPE_HEX}; // Hexahedrons only
   int shapeSizes[] = {8}; 
   int shapeCnt[] = {Nzones};
   int Nshapes = 1;
   
   DBPutZonelist2(fileptr,"zonelist",Nzones,Ndims,nodeList,Nzones*8,0,0,0,shapeTypes,shapeSizes,shapeCnt,Nshapes,NULL);
   
   char* crdNames[3];
   crdNames[0] = const_cast<char*>("X-axis");
   crdNames[1] = const_cast<char*>("Y-axis");
   crdNames[2] = const_cast<char*>("Z-axis");
   real* crdArrays[] = {x,y,z};
      
   DBPutUcdmesh(fileptr,gridName.c_str(),Ndims,crdNames,crdArrays,Nzones*8,Nzones,"zonelist",NULL,DB_FLOAT,NULL);
   
   if (DBPutUcdvar1(fileptr,varName.c_str(),gridName.c_str(),avgs,blockCntr,NULL,0,DB_FLOAT,DB_ZONECENT,NULL) < 0) return false;
   return true;
}

bool addSpatialCell(real* cellParams,creal& avg) {
   if (fileptr == NULL) return false;
   
   x[blockCntr*8 + 0] = cellParams[CellParams::XCRD];
   x[blockCntr*8 + 1] = cellParams[CellParams::XCRD] + cellParams[CellParams::DX];
   x[blockCntr*8 + 2] = cellParams[CellParams::XCRD] + cellParams[CellParams::DX];
   x[blockCntr*8 + 3] = cellParams[CellParams::XCRD];
   x[blockCntr*8 + 4] = cellParams[CellParams::XCRD];
   x[blockCntr*8 + 5] = cellParams[CellParams::XCRD] + cellParams[CellParams::DX];
   x[blockCntr*8 + 6] = cellParams[CellParams::XCRD] + cellParams[CellParams::DX];
   x[blockCntr*8 + 7] = cellParams[CellParams::XCRD];
   
   y[blockCntr*8 + 0] = cellParams[CellParams::YCRD];
   y[blockCntr*8 + 1] = cellParams[CellParams::YCRD];
   y[blockCntr*8 + 2] = cellParams[CellParams::YCRD] + cellParams[CellParams::DY];
   y[blockCntr*8 + 3] = cellParams[CellParams::YCRD] + cellParams[CellParams::DY];
   y[blockCntr*8 + 4] = cellParams[CellParams::YCRD];
   y[blockCntr*8 + 5] = cellParams[CellParams::YCRD];
   y[blockCntr*8 + 6] = cellParams[CellParams::YCRD] + cellParams[CellParams::DY];
   y[blockCntr*8 + 7] = cellParams[CellParams::YCRD] + cellParams[CellParams::DY];
   
   z[blockCntr*8 + 0] = cellParams[CellParams::ZCRD];
   z[blockCntr*8 + 1] = cellParams[CellParams::ZCRD];
   z[blockCntr*8 + 2] = cellParams[CellParams::ZCRD];
   z[blockCntr*8 + 3] = cellParams[CellParams::ZCRD];
   z[blockCntr*8 + 4] = cellParams[CellParams::ZCRD] + cellParams[CellParams::DZ];
   z[blockCntr*8 + 5] = cellParams[CellParams::ZCRD] + cellParams[CellParams::DZ];
   z[blockCntr*8 + 6] = cellParams[CellParams::ZCRD] + cellParams[CellParams::DZ];
   z[blockCntr*8 + 7] = cellParams[CellParams::ZCRD] + cellParams[CellParams::DZ];

   for (uint i=0; i<8; ++i) nodeList[blockCntr*8+i] = blockCntr*8+i;
   
   avgs[blockCntr] = avg;
   ++blockCntr;
   return true;
}

bool addVelocityGridBlock3D(real* blockParams) {
   if (fileptr == NULL) return false;
   
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      cuint cellIndex = k*WID2+j*WID+i;
      // Write vx-coordinates of the cell:
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 0] = blockParams[BlockParams::VXCRD] +   i  *blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 1] = blockParams[BlockParams::VXCRD] + (i+1)*blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 2] = blockParams[BlockParams::VXCRD] + (i+1)*blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 3] = blockParams[BlockParams::VXCRD] +   i  *blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 4] = blockParams[BlockParams::VXCRD] +   i  *blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 5] = blockParams[BlockParams::VXCRD] + (i+1)*blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 6] = blockParams[BlockParams::VXCRD] + (i+1)*blockParams[BlockParams::DVX];
      x[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 7] = blockParams[BlockParams::VXCRD] +   i  *blockParams[BlockParams::DVX];
      // vy-coordinates:
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 0] = blockParams[BlockParams::VYCRD] +   j  *blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 1] = blockParams[BlockParams::VYCRD] +   j  *blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 2] = blockParams[BlockParams::VYCRD] + (j+1)*blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 3] = blockParams[BlockParams::VYCRD] + (j+1)*blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 4] = blockParams[BlockParams::VYCRD] +   j  *blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 5] = blockParams[BlockParams::VYCRD] +   j  *blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 6] = blockParams[BlockParams::VYCRD] + (j+1)*blockParams[BlockParams::DVY];
      y[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 7] = blockParams[BlockParams::VYCRD] + (j+1)*blockParams[BlockParams::DVY];
      // vz-coordinates:
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 0] = blockParams[BlockParams::VZCRD] +   k  *blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 1] = blockParams[BlockParams::VZCRD] +   k  *blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 2] = blockParams[BlockParams::VZCRD] +   k  *blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 3] = blockParams[BlockParams::VZCRD] +   k  *blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 4] = blockParams[BlockParams::VZCRD] + (k+1)*blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 5] = blockParams[BlockParams::VZCRD] + (k+1)*blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 6] = blockParams[BlockParams::VZCRD] + (k+1)*blockParams[BlockParams::DVZ];
      z[blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + 7] = blockParams[BlockParams::VZCRD] + (k+1)*blockParams[BlockParams::DVZ];
      // Nodelist entries for the cell:
      for (uint counter=0; counter<8; ++counter) 
	nodeList[blockCntr*SIZE_VELBLOCK*8+cellIndex*8+counter] = blockCntr*SIZE_VELBLOCK*8 + cellIndex*8 + counter;
   }
   ++blockCntr;
   return true;
}

bool writeVelocityBlockGrid3D(const std::string& gridName,cuint& N_BLOCKS,real* blockParams) {
   if (fileptr == NULL) return false;

   allocateArrays(N_BLOCKS);
   for (uint b=0; b<N_BLOCKS; ++b) {
      addVelocityGridBlock3D(&(blockParams[b*SIZE_BLOCKPARAMS]));
   }
   
   writeReservedBlocks(gridName);
   deallocateArrays();
   return true;
}

bool writeVelocityBlockGridScalar3D(const std::string& varName,const std::string& gridName,cuint& N_BLOCKS,real* data) {
   if (fileptr == NULL) return false;
      
   if (DBPutUcdvar1(fileptr,varName.c_str(),gridName.c_str(),data,N_BLOCKS*SIZE_VELBLOCK,NULL,0,DB_FLOAT,DB_ZONECENT,NULL) < 0) return false;   
   return true;
}


