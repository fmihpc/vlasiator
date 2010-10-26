#ifndef WRITEVARS_H
#define WRITEVARS_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include "definitions.h"
#include "grid.h"
#include "devicegrid.h"

// Variable saving needs to be completely rewritten

//void writeSpatialCells(Grid& grid,DeviceGrid& deviceGrid,cuint& STEP);
//void writeBlocks(Grid& grid,DeviceGrid& deviceGrid,cuint& STEP,const std::vector<uint>& spaIndices);
void writeSpatialCells(Grid& grid,cuint& STEP);
void writeBlocks(Grid& grid,cuint& STEP,const std::vector<uint>& spaIndices);

void openOutputFile(const std::string& fileName);
void closeOutputFile(const std::string& fileName);

void writeVelocityBlockGrid(const std::string& gridName,cuint& BLOCKS,real* blockParams);
void writeVelocityBlockScalar(const std::string& varName,const std::string& gridName,cuint& BLOCKS,real* array);

#endif
