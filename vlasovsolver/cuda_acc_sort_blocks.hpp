/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */


#ifndef CUDA_SORT_BLOCKS_FOR_ACC_H
#define CUDA_SORT_BLOCKS_FOR_ACC_H

#include <vector>

#include "../common.h"
#include "../spatial_cell.hpp"

#include "include/splitvector/splitvec.h"

struct ColumnOffsets : public Managed {
   split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)

   ColumnOffsets(uint nColumns) {
      columnBlockOffsets.resize(nColumns);
      columnNumBlocks.resize(nColumns);
      setColumnOffsets.resize(nColumns);
      setNumColumns.resize(nColumns);
      std::cerr<<"Initializing columnOffsets with capacity "<<nColumns<<" = "<<setNumColumns.size()<<std::endl;
      columnBlockOffsets.clear();
      columnNumBlocks.clear();
      setColumnOffsets.clear();
      setNumColumns.clear();
   }
};

void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               const vmesh::VelocityMesh* vmesh,
                               const uint dimension,
                               uint* blocksGID,
                               uint* blocksLID,
                               ColumnOffsets* columnData,
                               // std::vector<uint> & columnBlockOffsets,
                               // std::vector<uint> & columnNumBlocks,
                               // std::vector<uint> & setColumnOffsets,
                               // std::vector<uint> & setNumColumns
                               cudaStream_t stream
   );

#endif
