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
      columnBlockOffsets.clear();
      columnNumBlocks.clear();
      setColumnOffsets.clear();
      setNumColumns.clear();

   }
   void dev_attachToStream(cudaStream_t stream = 0) {
      if (stream==0) {
         stream = cuda_getStream();
      }
      HANDLE_ERROR( cudaStreamAttachMemAsync(stream,this, 0,cudaMemAttachSingle) );
      columnBlockOffsets.streamAttach(stream);
      columnNumBlocks.streamAttach(stream);
      setColumnOffsets.streamAttach(stream);
      setNumColumns.streamAttach(stream);
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,columnBlockOffsets.data(), 0,cudaMemAttachSingle) );
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,columnNumBlocks.data(), 0,cudaMemAttachSingle) );
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,setColumnOffsets.data(), 0,cudaMemAttachSingle) );
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,setNumColumns.data(), 0,cudaMemAttachSingle) );
   }
   void dev_detachFromStream() {
      cudaStream_t stream = 0;
      HANDLE_ERROR( cudaStreamAttachMemAsync(stream,this, 0,cudaMemAttachGlobal) );
      columnBlockOffsets.streamAttach(0,cudaMemAttachGlobal);
      columnNumBlocks.streamAttach(0,cudaMemAttachGlobal);
      setColumnOffsets.streamAttach(0,cudaMemAttachGlobal);
      setNumColumns.streamAttach(0,cudaMemAttachGlobal);
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,columnBlockOffsets.data(), 0,cudaMemAttachGlobal) );
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,columnNumBlocks.data(), 0,cudaMemAttachGlobal) );
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,setColumnOffsets.data(), 0,cudaMemAttachGlobal) );
      // HANDLE_ERROR( cudaStreamAttachMemAsync(stream,setNumColumns.data(), 0,cudaMemAttachGlobal) );
   }
};

void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               vmesh::VelocityMesh* vmesh,
                               const vmesh::LocalID nBlocks,
                               const uint dimension,
                               vmesh::GlobalID *blocksID_mapped,
                               vmesh::GlobalID *blocksID_mapped_sorted,
                               vmesh::GlobalID *blocksGID,
                               vmesh::LocalID *blocksLID_unsorted,
                               vmesh::LocalID *blocksLID,
                               vmesh::LocalID *dev_columnNBlocks,
                               ColumnOffsets* columnData,
                               // std::vector<uint> & columnBlockOffsets,
                               // std::vector<uint> & columnNumBlocks,
                               // std::vector<uint> & setColumnOffsets,
                               // std::vector<uint> & setNumColumns
                               const uint cuda_async_queue_id,
                               cudaStream_t stream
   );

#endif
