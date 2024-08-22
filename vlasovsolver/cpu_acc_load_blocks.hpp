#ifndef CPU_ACC_LOAD_BLOCKS_H
#define CPU_ACC_LOAD_BLOCKS_H

#include "../common.h"
#include "../spatial_cell_wrapper.hpp"
#include "vec.h"


//index in the temporary and padded column data values array. Each
//column has an empty block in ether end.
#define i_pcolumnv(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )
#define i_pcolumnv_b(planeVectorIndex, k, k_block, num_k_blocks) ( planeVectorIndex * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

void loadColumnBlockData(
   const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,
   vmesh::GlobalID* blocks,
   vmesh::LocalID n_blocks,
   const int dimension,
   Vec* __restrict__ values);



#endif
