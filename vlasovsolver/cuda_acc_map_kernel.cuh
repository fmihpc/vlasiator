#include "vec.h"
#include "../common.h"

struct Column
{
   int valuesOffset;                              // Source data values
   size_t targetBlockOffsets[MAX_BLOCKS_PER_DIM]; // Target data array offsets
   int nblocks;                                   // Number of blocks in this column
   int minBlockK,maxBlockK;                       // Column parallel coordinate limits
   int kBegin;                                    // Actual un-sheared starting block index
   int i,j;                                       // Blocks' perpendicular coordinates
};

extern Realf* acceleration_1_wrapperCaller
(
  Realf *blockData,
  Column *columns,
  Vec *values,
  uint cell_indices_to_id[],
  int totalColumns,
  int valuesSizeRequired,
  int bdsw3,
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Realv v_min,
  Realv i_dv,
  Realv dv,
  Realv minValue
);
extern Realf* acceleration_1_wrapper
(
  Realf *blockData,
  Column *columns,
  Vec *values,
  uint cell_indices_to_id[],
  int totalColumns,
  int valuesSizeRequired,
  int bdsw3,
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Realv v_min,
  Realv i_dv,
  Realv dv,
  Realv minValue
);

#define DIMS 1
#define BLOCKS 256
#define THREADS 4
#define CUDASIZE 1
