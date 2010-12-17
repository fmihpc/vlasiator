#ifndef SILOWRITER_H
#define SILOWRITER_H

#include "definitions.h"

bool openOutputFile(const std::string& fileName,const std::string& fileInfo);
bool closeOutputFile();

//bool addVelocityBlock(real* blockParams,real* block);
bool addSpatialCell(Real* cellParams);
bool addVelocityGridBlock3D(Real* blockParams);
bool freeBlocks();
bool freeCells();
bool reserveSpatialCells(cuint& CELLS);
bool reserveVelocityBlocks(cuint& BLOCKS);
//bool writeBlocks(const std::string& gridName,const std::string& varName);
bool writeReservedBlocks(const std::string& gridName);
bool writeSpatialCells(const std::string& gridName);
bool writeVelocityBlockGrid3D(const std::string& gridName,cuint& N_BLOCKS,Real* blockParams);
bool writeVelocityBlockGridScalar3D(const std::string& varName,const std::string& gridName,cuint& N_BLOCKS,Real* data);



#endif
