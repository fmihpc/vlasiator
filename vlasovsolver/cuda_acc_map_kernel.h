#include "vec.h"
#include "../common.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

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
   Realf* dev_blockData,
   Vec* dev_blockDataOrdered,
   uint dev_cell_indices_to_id[],
   Column* dev_columns,
   const int totalColumns,
   const Realv intersection,
   const Realv intersection_di,
   const Realv intersection_dj,
   const Realv intersection_dk,
   const Realv v_min,
   const Realv i_dv,
   const Realv dv,
   const Realv minValue,
   const int bdsw3,
   const int cudablocks,
   const int cudathreads,
   cudaStream_t stream
);

extern void reorder_blocks_by_dimension_glue(
   Realf* dev_blockData,
   Vec* dev_blockDataOrdered,
   uint dev_cell_indices_to_id[],
   const uint totalColumns,
   uint* dev_LIDlist,
   uint* dev_columnNumBlocks,
   uint* dev_columnBlockOffsets,
   const int cudablocks, 
   const int cudathreads,
   cudaStream_t stream);

extern void cuda_acc_allocate (
   uint maxBlockCount
   );
extern void cuda_acc_allocate_memory (
   uint cpuThreadID,
   uint maxBlockCount
   );
extern void cuda_acc_deallocate_memory (
   uint cpuThreadID
   );

// Device data variables, to be allocated in good time. Made into an array so that each thread has their own pointer.
//extern Realf *dev_blockData[]; //now inside velocity_block_container.h
extern Vec *dev_blockDataOrdered[];
extern Column *dev_columns[];
extern uint *dev_cell_indices_to_id[];
//extern uint *dev_GIDlist[];
extern uint *dev_LIDlist[];
extern uint *dev_columnNumBlocks[];
extern uint *dev_columnBlockOffsets[];

// Host data variables, to be allocated and pinned in good time. Made into a long array so that each thread has their own pointer.
extern Column *host_columns[];
extern uint *host_GIDlist[];
extern uint *host_LIDlist[];

static const double CUDA_ACC_SAFECTY_FACTOR = 0.75;

extern uint cuda_acc_allocatedSize;
extern uint cuda_acc_allocatedColumns;
extern Real cudaAllocationMultiplier;


