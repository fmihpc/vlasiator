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

extern double* acceleration_1_wrapperCaller
(
  double *blockData,
  Column *columns,
  Vec *values,
  uint cell_indices_to_id[],
  int totalColumns,
  int valuesSizeRequired,
  int bdsw3,
  double intersection,
  double intersection_di,
  double intersection_dj,
  double intersection_dk,
  double v_min,
  double i_dv,
  double dv,
  double minValue
);
extern double* acceleration_1_wrapper
(
  double *blockData,
  Column *columns,
  Vec *values,
  uint cell_indices_to_id[],
  int totalColumns,
  int valuesSizeRequired,
  int bdsw3,
  double intersection,
  double intersection_di,
  double intersection_dj,
  double intersection_dk,
  double v_min,
  double i_dv,
  double dv,
  double minValue
);

#define DIMS 4
#define BLOCKS 2
#define THREADS 4
#define CUDASIZE 4
