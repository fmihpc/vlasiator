#include "cuda_header.cuh"
#include "open_acc_map_h.cuh"
#include "../vlasovsolver/vec.h"
#include "../definitions.h"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NPP_MAXABS_32F ( 3.402823466e+38f )
#define NPP_MINABS_32F ( 1.175494351e-38f )
#define NPP_MAXABS_64F ( 1.7976931348623158e+308 )
#define NPP_MINABS_64F ( 2.2250738585072014e-308 )

//#define i_pcolumnv_cuda(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}
/*
__device__ Vec minmod(const Vec slope1, const Vec slope2)
{
  const Vec zero(0.0);
  Vec slope = select(abs(slope1) < abs(slope2), slope1, slope2);
  return select(slope1 * slope2 <= 0, zero, slope);
}
__device__ Vec maxmod(const Vec slope1, const Vec slope2)
{
  const Vec zero(0.0);
  Vec slope = select(abs(slope1) > abs(slope2), slope1, slope2);
  return select(slope1 * slope2 <= 0, zero, slope);
}
__device__ Vec slope_limiter_sb(const Vec &l, const Vec &m, const Vec &r)
{
  Vec a = r-m;
  Vec b = m-l;
  const Vec slope1 = minmod(a, 2*b);
  const Vec slope2 = minmod(2*a, b);
  return maxmod(slope1, slope2);
}
__device__ Vec slope_limiter(const Vec &l, const Vec &m,const Vec &r)
{
   return slope_limiter_sb(l,m,r);
}
__device__ void compute_plm_coeff(const Vec * const values, uint k, Vec a[2], const Realv threshold)
{
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realv scale = 1./threshold;
  //Vec v_1 = values[k - 1] * scale;
  //Vec v_2 = values[k] * scale;
  //Vec v_3 = values[k + 1] * scale;
  //Vec d_cv = slope_limiter(v_1, v_2, v_3) * threshold;
  const Vec d_cv = slope_limiter( values[k-1]*scale, values[k]*scale, values[k+1]*scale)*threshold;
  a[0] = values[k] - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}
*/
__global__ void acceleration_1
(
  double *dev_blockData,
  Column *dev_columns,
  Vec dev_target_density,
  Veci dev_target_cell,
  const int blockK,
  const uint column
)
{
  int index = threadIdx.x + blockIdx.x*blockDim.x;
  if(index == 0)
  {
    for (int target_i=0; target_i < VECL; ++target_i)
    {
      const Realf tval = dev_target_density[target_i];
      const uint tcell = dev_target_cell[target_i];
      (&dev_blockData[dev_columns[column].targetBlockOffsets[blockK]])[tcell] += tval;
    }
  }
}

Realf* acceleration_1_wrapper
(
  Realf *blockData,
  Column *columns,
  Vec target_density,
  Veci target_cell,
  const int totalColumns,
  const int blockK,
  const uint column,
  const int bdsw3
)
{
  double *dev_blockData;
  HANDLE_ERROR( cudaMalloc((void**)&dev_blockData, bdsw3*sizeof(double)) );
  HANDLE_ERROR( cudaMemcpy(dev_blockData, blockData, bdsw3*sizeof(double), cudaMemcpyHostToDevice) );

  Column *dev_columns;
  HANDLE_ERROR( cudaMalloc((void**)&dev_columns, totalColumns*sizeof(Column)) );
  HANDLE_ERROR( cudaMemcpy(dev_columns, columns, totalColumns*sizeof(Column), cudaMemcpyHostToDevice) );

  /*
  Vec *dev_target_density;
  HANDLE_ERROR( cudaMalloc((void**)&dev_target_density, VECL*sizeof(Vec)) );
  HANDLE_ERROR( cudaMemcpy(dev_target_density, target_density, VECL*sizeof(Vec), cudaMemcpyHostToDevice) );

  Vec *dev_target_cell;
  HANDLE_ERROR( cudaMalloc((void**)&dev_target_cell, VECL*sizeof(Vec)) );
  HANDLE_ERROR( cudaMemcpy(dev_target_cell, target_cell, VECL*sizeof(Vec), cudaMemcpyHostToDevice) );
  */
  acceleration_1<<<BLOCKS, THREADS>>>
  (
    dev_blockData,
    dev_columns,
    target_density,
    target_cell,
    blockK,
    column
  );

  cudaDeviceSynchronize();
  HANDLE_ERROR( cudaMemcpy(blockData, dev_blockData, bdsw3*sizeof(double), cudaMemcpyDeviceToHost) );

  HANDLE_ERROR( cudaFree(dev_blockData) );
  HANDLE_ERROR( cudaFree(dev_columns) );
  //HANDLE_ERROR( cudaFree(dev_target_density) );
  //HANDLE_ERROR( cudaFree(dev_target_cell) );

  return blockData;
}
