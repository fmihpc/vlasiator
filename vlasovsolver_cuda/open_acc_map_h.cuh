#include "../vlasovsolver/vec.h"
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

struct Acceleration_1_struct
{
  Realf *blockData;
  int totalColumns;
  Column *columns;
  Vec *values;
  uint cell_indices_to_id[];
  Realv intersection;
  Realv intersection_di;
  Realv intersection_dj;
  Realv intersection_dk;
  Real minValue;
  Realv dv;
  Realv v_min;
};

extern Acceleration_1_struct acceleration_1_wrapperCaller
(
  int bdsw3,
  Realf *blockData,
  int totalColumns,
  Column *columns,
  int valuesSizeRequired,
  Vec *values,
  uint cell_indices_to_id[],
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Real minValue,
  Realv dv,
  Realv v_min
);
extern Acceleration_1_struct acceleration_1_wrapper
(
  int bdsw3,
  Realf *blockData,
  int totalColumns,
  Column *columns,
  int valuesSizeRequired,
  Vec *values,
  uint cell_indices_to_id[],
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Real minValue,
  Realv dv,
  Realv v_min
);
