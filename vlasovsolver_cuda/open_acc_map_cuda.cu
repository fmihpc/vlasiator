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

#define i_pcolumnv_cuda(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}
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
__device__ Vec slope_limiter_sb(Vec &l, Vec &m, Vec &r)
{
  Vec a = r-m;
  Vec b = m-l;
  const Vec slope1 = minmod(a, 2*b);
  const Vec slope2 = minmod(2*a, b);
  return maxmod(slope1, slope2);
}
__device__ Vec slope_limiter(Vec &l,Vec &m,Vec &r)
{
   return slope_limiter_sb(l,m,r);
}
//changed to Real it was Realv
__device__ void compute_plm_coeff(const Vec * const values, uint k, Vec a[2], const Real threshold)
{
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realv scale = 1./threshold;

  Vec v_1 = values[k - 1] * scale;
  Vec v_2 = values[k] * scale;
  Vec v_3 = values[k + 1] * scale;
  Vec d_cv = slope_limiter(v_1, v_2, v_3) * threshold;
  //const Vec d_cv = slope_limiter( values[k-1]*scale, values[k]*scale, values[k+1]*scale)*threshold;
  a[0] = values[k] - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}

__global__ void acceleration_1
(
  double *dev_blockData,
  Column *dev_columns,
  Vec *dev_values,
  int *dev_cell_indices_to_id,
  int totalColumns,
  double intersection,
  double intersection_di,
  double intersection_dj,
  double intersection_dk,
  double v_min,
  double i_dv,
  double dv,
  double minValue,
  int acc_semilag_flag,
  int bdsw3
)
{
  int index = threadIdx.x + blockIdx.x*blockDim.x;
  if(index == 0)
  {
    //printf("CUDA 1\n");
    for( uint column=0; column < totalColumns; column++)
    {
      //printf("CUDA 2\n");
      // i,j,k are relative to the order in which we copied data to the values array.
      // After this point in the k,j,i loops there should be no branches based on dimensions
      // Note that the i dimension is vectorized, and thus there are no loops over i
      // Iterate through the perpendicular directions of the column
       for (uint j = 0; j < WID; j += VECL/WID)
       {
         //printf("CUDA 3; VECL = %d\n", VECL);
          const vmesh::LocalID nblocks = dev_columns[column].nblocks;
          // create vectors with the i and j indices in the vector position on the plane.
          #if VECL == 4
            const Veci i_indices = Veci(0, 1, 2, 3);
            const Veci j_indices = Veci(j, j, j, j);
          #elif VECL == 8
            const Veci i_indices = Veci(0, 1, 2, 3, 0, 1, 2, 3);
            const Veci j_indices = Veci(j, j, j, j, j + 1, j + 1, j + 1, j + 1);
          #elif VECL == 16
            const Veci i_indices = Veci(0, 1, 2, 3,
                                        0, 1, 2, 3,
                                        0, 1, 2, 3,
                                        0, 1, 2, 3);
            const Veci j_indices = Veci(j, j, j, j,
                                        j + 1, j + 1, j + 1, j + 1,
                                        j + 2, j + 2, j + 2, j + 2,
                                        j + 3, j + 3, j + 3, j + 3);
          #endif
          const Veci  target_cell_index_common =
             i_indices * dev_cell_indices_to_id[0] +
             j_indices * dev_cell_indices_to_id[1];
          // intersection_min is the intersection z coordinate (z after
          // swaps that is) of the lowest possible z plane for each i,j
          // index (i in vector)
          const Vec intersection_min =
             intersection +
             (dev_columns[column].i * WID + to_realv(i_indices)) * intersection_di +
             (dev_columns[column].j * WID + to_realv(j_indices)) * intersection_dj;

          /*compute some initial values, that are used to set up the
           * shifting of values as we go through all blocks in
           * order. See comments where they are shifted for
           * explanations of their meaning*/
          Vec v_r0( (WID * dev_columns[column].kBegin) * dv + v_min);
          Vec lagrangian_v_r0((v_r0-intersection_min)/intersection_dk);
          /* compute location of min and max, this does not change for one
          column (or even for this set of intersections, and can be used
          to quickly compute max and min later on*/
          //TODO, these can be computed much earlier, since they are
          //identiacal for each set of intersections
          int minGkIndex=0, maxGkIndex=0; // 0 for compiler
          {
            Realv maxV = (sizeof(Realv) == 4) ? NPP_MINABS_32F : NPP_MINABS_64F;
            Realv minV = (sizeof(Realv) == 4) ? NPP_MAXABS_32F : NPP_MAXABS_64F;
             for(int i = 0; i < VECL; i++)
             {
                if (lagrangian_v_r0[i] > maxV)
                {
                   maxV = lagrangian_v_r0[i];
                   maxGkIndex = i;
                }
                if (lagrangian_v_r0[i] < minV)
                {
                   minV = lagrangian_v_r0[i];
                   minGkIndex = i;
                }
             }
          }
          // loop through all blocks in column and compute the mapping as integrals.
          for (uint k=0; k < WID * nblocks; ++k)
          {
             //printf("CUDA 4\n");
             // Compute reconstructions
             // values + i_pcolumnv(n_cblocks, -1, j, 0) is the starting point of the column data for fixed j
             // k + WID is the index where we have stored k index, WID amount of padding.
             //if(acc_semilag_flag==0)
             //{
              Vec a[2];
              compute_plm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), k + WID, a, minValue);
              //}
              /*
              if(acc_semilag_flag==1)
              {
                Vec *a = new Vec[3];
                //compute_ppm_coeff(values + dev_columns[column].valuesOffset  + i_pcolumnv(j, 0, -1, nblocks), h4, k + WID, a, minValue);
              }
              if(acc_semilag_flag==2)
              {
                Vec *a = new Vec[5];
                //compute_pqm_coeff(values + dev_columns[column].valuesOffset  + i_pcolumnv(j, 0, -1, nblocks), h8, k + WID, a, minValue);
              }
              */
             // set the initial value for the integrand at the boundary at v = 0
             // (in reduced cell units), this will be shifted to target_density_1, see below.
             Vec target_density_r(0.0);
             // v_l, v_r are the left and right velocity coordinates of source cell.
             Vec v_r = v_r0  + (k+1)* dv;
             Vec v_l = v_r0  + k* dv;
             // left(l) and right(r) k values (global index) in the target
             // Lagrangian grid, the intersecting cells. Again old right is new left.

             // I keep only this version with Fallback, because the version with Agner requires another call to CPU
             Veci lagrangian_gk_l = truncate_to_int((v_l-intersection_min)/intersection_dk);
             Veci lagrangian_gk_r = truncate_to_int((v_r-intersection_min)/intersection_dk);
              /*
              if(VECTORCLASS_H_DEVICE >= 20000)
              {
                lagrangian_gk_r = truncatei((v_l-intersection_min)/intersection_dk);
                lagrangian_gk_r = truncatei((v_r-intersection_min)/intersection_dk);
              }
              else
              {
                lagrangian_gk_l = truncate_to_int((v_l-intersection_min)/intersection_dk);
                lagrangian_gk_r = truncate_to_int((v_r-intersection_min)/intersection_dk);
              }
              */
             //limits in lagrangian k for target column. Also take into
             //account limits of target column
             int minGk = max(int(lagrangian_gk_l[minGkIndex]), int(dev_columns[column].minBlockK * WID));
             int maxGk = min(int(lagrangian_gk_r[maxGkIndex]), int((dev_columns[column].maxBlockK + 1) * WID - 1));
             // Run along the column and perform the polynomial reconstruction
             for(int gk = dev_columns[column].minBlockK * WID; gk <= dev_columns[column].maxBlockK * WID; gk++)
             {
                //printf("\nminGk = %d; gk = %d; maxGk = %d;\n", minGk, gk, maxGk);
                if(gk < minGk || gk > maxGk)
                { continue; }
                //printf("CUDA 6\n");
                const int blockK = gk/WID;
                const int gk_mod_WID = (gk - blockK * WID);
                //the block of the Lagrangian cell to which we map
                //const int target_block(target_block_index_common + blockK * block_indices_to_id[2]);
                //cell indices in the target block  (TODO: to be replaced by
                //compile time generated scatter write operation)
                const Veci target_cell(target_cell_index_common + gk_mod_WID * dev_cell_indices_to_id[2]);
                //the velocity between which we will integrate to put mass
                //in the targe cell. If both v_r and v_l are in same cell
                //then v_1,v_2 should be between v_l and v_r.
                //v_1 and v_2 normalized to be between 0 and 1 in the cell.
                //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
                const Vec v_norm_r = (  min(  max( (gk + 1) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;
                /*shift, old right is new left*/
                const Vec target_density_l = target_density_r;
                // compute right integrand
                if(acc_semilag_flag==0)
                  target_density_r = v_norm_r * ( a[0] + v_norm_r * a[1] );
                //if(acc_semilag_flag==1)
                //  target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * a[2] ) );
                //if(acc_semilag_flag==2)
                //  target_density_r =
                //    v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * ( a[2] + v_norm_r * ( a[3] + v_norm_r * a[4] ) ) ) );
                //store values, one element at a time. All blocks have been created by now.
                //TODO replace by vector version & scatter & gather operation
                const Vec target_density = target_density_r - target_density_l;
                //printf("CUDA 7\n");
                //int a = dev_columns[column].targetBlockOffsets[blockK];
                for (int target_i=0; target_i < VECL; ++target_i)
                {
                  //printf("CUDA 8\n");
                  // do the conversion from Realv to Realf here, faster than doing it in accumulation
                  const Realf tval = target_density[target_i];
                  const uint tcell = target_cell[target_i];
                  //printf("CUDA: tval = %.2f; tcell = %d;\n", tval, tcell);
                  //printf("&dev_blockData[a] = %.2f; tcell = %d\n", dev_blockData[dev_columns[column].targetBlockOffsets[blockK]], tcell );
                  (&dev_blockData[dev_columns[column].targetBlockOffsets[blockK]])[tcell] += tval;
                  //for(uint cell=0; cell<bdsw3; cell++)
                  //{
                  //  printf("blockData[cell] = %.2f\n", blockData[cell]);
                  //}
                  //for(uint aa=0; aa<bdsw3; aa++)
                  //{
                  //    printf("dev_blockData[%d] = %.2f\n", aa, dev_blockData[aa]);
                  //}
                  //printf("CUDA 11\n");
                }  // for-loop over vector elements
             } // for loop over target k-indices of current source block
          } // for-loop over source blocks
       } //for loop over j index
    } //for loop over columns
  }
}

double* acceleration_1_wrapper
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
)
{
  int acc_semilag_flag = 0;
  #ifdef ACC_SEMILAG_PLM
    acc_semilag_flag = 0;
  #endif
  #ifdef ACC_SEMILAG_PPM
    acc_semilag_flag = 1;
  #endif
  #ifdef ACC_SEMILAG_PQM
    acc_semilag_flag = 2;
  #endif

  double *dev_blockData;
  HANDLE_ERROR( cudaMalloc((void**)&dev_blockData, bdsw3*sizeof(double)) );
  HANDLE_ERROR( cudaMemcpy(dev_blockData, blockData, bdsw3*sizeof(double), cudaMemcpyHostToDevice) );

  Column *dev_columns;
  HANDLE_ERROR( cudaMalloc((void**)&dev_columns, totalColumns*sizeof(Column)) );
  HANDLE_ERROR( cudaMemcpy(dev_columns, columns, totalColumns*sizeof(Column), cudaMemcpyHostToDevice) );

  int *dev_cell_indices_to_id;
  HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id, 3*sizeof(int)) );
  HANDLE_ERROR( cudaMemcpy(dev_cell_indices_to_id, cell_indices_to_id, 3*sizeof(int), cudaMemcpyHostToDevice) );

  Vec *dev_values;
  HANDLE_ERROR( cudaMalloc((void**)&dev_values, valuesSizeRequired*sizeof(Vec)) );
  HANDLE_ERROR( cudaMemcpy(dev_values, values, valuesSizeRequired*sizeof(Vec), cudaMemcpyHostToDevice) );

  acceleration_1<<<BLOCKS, THREADS>>>
  (
    dev_blockData,
    dev_columns,
    dev_values,
    dev_cell_indices_to_id,
        totalColumns,
        intersection,
        intersection_di,
        intersection_dj,
        intersection_dk,
        v_min,
        i_dv,
        dv,
        minValue,
        acc_semilag_flag,
        bdsw3
  );

  cudaDeviceSynchronize();
  HANDLE_ERROR( cudaMemcpy(blockData, dev_blockData, bdsw3*sizeof(double), cudaMemcpyDeviceToHost) );

  HANDLE_ERROR( cudaFree(dev_blockData) );
  HANDLE_ERROR( cudaFree(dev_cell_indices_to_id) );
  HANDLE_ERROR( cudaFree(dev_columns) );
  HANDLE_ERROR( cudaFree(dev_values) );

  return blockData;
}
