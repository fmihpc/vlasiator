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

extern void acceleration_1_glue(
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
  Realv minValue,
  const uint cuda_async_queue_id
);

#define DIMS 1
#ifndef CUDABLOCKS
#  define CUDABLOCKS (64)
#endif
#ifndef CUDATHREADS
#  define CUDATHREADS (32) // NVIDIA: 32 AMD: 64
#endif

extern void cuda_acc_allocate_memory (
   uint cpuThreadID,
   uint maxBlockCount
   );
extern void cuda_acc_deallocate_memory (
   uint cpuThreadID
   );

// Device data variables, to be allocated in good time. Made into a long array so that each thread has their own pointer.
extern Realf *dev_blockData[];
extern Column *dev_columns[];
extern int *dev_cell_indices_to_id[];
extern Vec *dev_values[];

// Host data variables, to be allocated and pinned in good time. Made into a long array so that each thread has their own pointer.
extern Column *host_columns[];
extern Vec *host_values[];
extern bool isCudaAllocated;
extern uint cudaMaxBlockCount;
extern float cudaAllocationMultiplier;

