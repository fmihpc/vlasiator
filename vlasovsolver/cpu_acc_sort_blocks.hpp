
/*
This file is part of Vlasiator.
Copyright 2013-2015 Finnish Meteorological Institute

*/

#ifndef CPU_SORT_BLOCKS_FOR_ACC_H
#define CPU_SORT_BLOCKS_FOR_ACC_H

#include <vector>

#include "../common.h"
#include "../spatial_cell.hpp"

void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell, 
                               const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                               const uint dimension,
                               uint* blocks,
                               std::vector<uint> & columnBlockOffsets,
                               std::vector<uint> & columnNumBlocks,
                               std::vector<uint> & setColumnOffsets,
                               std::vector<uint> & setNumColumns);

#endif
