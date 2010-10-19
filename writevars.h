#ifndef WRITEVARS_H
#define WRITEVARS_H

#include <vector>
#include "definitions.h"
#include "grid.h"
#include "devicegrid.h"

// Variable saving needs to be completely rewritten

void writeSpatialCells(Grid& grid,DeviceGrid& deviceGrid,cuint& STEP);
void writeBlocks(Grid& grid,DeviceGrid& deviceGrid,cuint& STEP,const std::vector<uint>& spaIndices);

#endif
