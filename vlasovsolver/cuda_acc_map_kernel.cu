#include "cuda_header.cuh"
#include "cuda_acc_map_kernel.cuh"
#include "vec.h"
#include "../definitions.h"
#include "cuda_face_estimates.cuh"
#include "cuda_1d_pqm.cuh"
#include "cuda_1d_ppm.cuh"
#include "cuda_1d_plm.cuh"

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
#define i_pcolumnv_cuda_realf(j, k, k_block, num_k_blocks, index) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) * VECL + (k) * VECL + ( k_block + 1 ) * WID * VECL + (index) )

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

__global__ void acceleration_1
(
  Realf *dev_blockData,
  Column *dev_columns,
  Vec *dev_values,
  int *dev_cell_indices_to_id,
  int totalColumns,
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Realv v_min,
  Realv i_dv,
  Realv dv,
  Realv minValue,
  int bdsw3
)
{
   //int index = threadIdx.x + blockIdx.x*blockDim.x;
   //if(totalColumns > 256)
   //  printf("totalColumns = %d;\n", totalColumns);
   //int column = threadIdx.x + blockIdx.x * blockDim.x;
   //if(column < totalColumns)
   //for( uint column=0; column < totalColumns; column++)

#ifdef CUDA_REALF
   printf("CUDA REALF\n");

   // Optimization idea for future: reverse-sort columns based on
   // block counts
   int index = threadIdx.x;
   int block = blockIdx.x;
   int nBlocks = gridDim.x;
   // How many columns max per block?
   int maxcolumns = (int)ceil((Realv)totalColumns / (Realv)nBlocks);
   
   //Realf * dev_values_realf = reinterpret_cast<Realf*>(dev_values);
//   Realf * dev_values_realf = &dev_values[0][0];
//   printf("Maxcolumns %d totalcolumns %d nblocks %d dev_values %d %d %f %f\n", maxcolumns, totalColumns, nBlocks, &dev_values[0][0], &dev_values_realf[0], dev_values[0][0], dev_values_realf[0]);

   for (uint blockC = 0; blockC < maxcolumns; ++blockC) {
      int column = blockC*nBlocks + block;
      if (column >= totalColumns) continue;
      /* New threading with each warp/wavefront working on one vector */
      Realf v_r0 = ( (WID * dev_columns[column].kBegin) * dv + v_min);

      // i,j,k are relative to the order in which we copied data to the values array.
      // After this point in the k,j,i loops there should be no branches based on dimensions
      // Note that the i dimension is vectorized, and thus there are no loops over i
      // Iterate through the perpendicular directions of the column
      for (uint j = 0; j < WID; j += VECL/WID) {
         // If VECL=WID2 (WID=4, VECL=16, or WID=8, VECL=64, j==0)
         // This loop is still needed for e.g. Warp=VECL=32, WID2=64
         const vmesh::LocalID nblocks = dev_columns[column].nblocks;

         uint i_indices = index % WID;
         uint j_indices = j + index/WID;

         int target_cell_index_common =
            i_indices * dev_cell_indices_to_id[0] +
            j_indices * dev_cell_indices_to_id[1];
         const Realf intersection_min =
            intersection +
            (dev_columns[column].i * WID + (Realv)i_indices) * intersection_di +
            (dev_columns[column].j * WID + (Realv)j_indices) * intersection_dj;

         Realf lagrangian_v_r0 = ((v_r0-intersection_min)/intersection_dk);
         // loop through all perpendicular slices in column and compute the mapping as integrals.
         for (uint k=0; k < WID * nblocks; ++k)
         {
            // Compute reconstructions
#ifdef ACC_SEMILAG_PLM
            Realv a[2];
            //compute_plm_coeff(dev_values_realf + dev_columns[column].valuesOffset * VECL + i_pcolumnv_cuda_realf(j, 0, -1, nblocks, index), (k + WID), a, minValue);
            compute_plm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), (k + WID), a, minValue, index);
#endif
#ifdef ACC_SEMILAG_PPM
            Realv a[3];
          //compute_ppm_coeff(dev_values_realf + dev_columns[column].valuesOffset * VECL + i_pcolumnv_cuda_realf(j, 0, -1, nblocks, index), h4, (k + WID), a, minValue);
            compute_ppm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h4, (k + WID), a, minValue, index);
#endif
#ifdef ACC_SEMILAG_PQM
            Realv a[5];
          //compute_pqm_coeff(dev_values_realf + dev_columns[column].valuesOffset * VECL + i_pcolumnv_cuda_realf(j, 0, -1, nblocks, index), h8, (k + WID), a, minValue);
            compute_pqm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h8, (k + WID), a, minValue, index);
#endif

            // set the initial value for the integrand at the boundary at v = 0
            // (in reduced cell units), this will be shifted to target_density_1, see below.
            Realf target_density_r = 0.0;

            const Realv v_r = v_r0  + (k+1)* dv;
            const Realv v_l = v_r0  + k* dv;
            const int lagrangian_gk_l = trunc((v_l-intersection_min)/intersection_dk);
            const int lagrangian_gk_r = trunc((v_r-intersection_min)/intersection_dk);
            
            //limits in lagrangian k for target column. Also take into
            //account limits of target column
            // TODO: Would it be more efficient to force the same gk-loop for all indexes in the warp?
            const int minGk = max(lagrangian_gk_l, int(dev_columns[column].minBlockK * WID));
            const int maxGk = min(lagrangian_gk_r, int((dev_columns[column].maxBlockK + 1) * WID - 1));
            // Run along the column and perform the polynomial reconstruction
            for(int gk = minGk; gk <= maxGk; gk++) {
               const int blockK = gk/WID;
               const int gk_mod_WID = (gk - blockK * WID);

               //the block of the Lagrangian cell to which we map
               //const int target_block(target_block_index_common + blockK * block_indices_to_id[2]);
               // This already contains the value index via target_cell_index_commom
               const int tcell(target_cell_index_common + gk_mod_WID * dev_cell_indices_to_id[2]);
               //the velocity between which we will integrate to put mass
               //in the targe cell. If both v_r and v_l are in same cell
               //then v_1,v_2 should be between v_l and v_r.
               //v_1 and v_2 normalized to be between 0 and 1 in the cell.
               //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
               const Realf v_norm_r = (  min(  max( (gk + 1) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;

               /*shift, old right is new left*/
               const Realf target_density_l = target_density_r;

                // compute right integrand
                #ifdef ACC_SEMILAG_PLM
                   target_density_r = v_norm_r * ( a[0] + v_norm_r * a[1] );
                #endif
                #ifdef ACC_SEMILAG_PPM
                   target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * a[2] ) );
                #endif
                #ifdef ACC_SEMILAG_PQM
                   target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * ( a[2] + v_norm_r * ( a[3] + v_norm_r * a[4] ) ) ) );
                #endif

                //store values, one element at a time. All blocks have been created by now.
                const Realf tval = target_density_r - target_density_l;
                (&dev_blockData[dev_columns[column].targetBlockOffsets[blockK]])[tcell] += tval;
            } // for loop over target k-indices of current source block
         } // for-loop over source blocks
      } //for loop over j index
   } //for loop over columns

#else // NOT CUDA_REALF
   printf("old CUDA\n");
   int column = threadIdx.x + blockIdx.x * blockDim.x;
   //for( uint column=0; column < totalColumns; column++)
   if (column < totalColumns) {
      // i,j,k are relative to the order in which we copied data to the values array.
      // After this point in the k,j,i loops there should be no branches based on dimensions
      // Note that the i dimension is vectorized, and thus there are no loops over i
      // Iterate through the perpendicular directions of the column
      for (uint j = 0; j < WID; j += VECL/WID) {
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
         #elif VECL == 64 // assumes WID=8
            const Veci i_indices = Veci(0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7,
                                        0, 1, 2, 3, 4, 5, 6, 7);
            const Veci j_indices = Veci(j, j, j, j, j, j, j, j,
                                        j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1,
                                        j + 2, j + 2, j + 2, j + 2, j + 2, j + 2, j + 2, j + 2,
                                        j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, 
                                        j + 4, j + 4, j + 4, j + 4, j + 4, j + 4, j + 4, j + 4, 
                                        j + 5, j + 5, j + 5, j + 5, j + 5, j + 5, j + 5, j + 5, 
                                        j + 6, j + 6, j + 6, j + 6, j + 6, j + 6, j + 6, j + 6, 
                                        j + 7, j + 7, j + 7, j + 7, j + 7, j + 7, j + 7, j + 7, 
);
         #endif

         /* array for converting block indices to id using a dot product, 
            depends on Cartesian direction*/
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
            // Compute reconstructions
#ifdef ACC_SEMILAG_PLM
            Vec a[2];
            compute_plm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), k + WID, a, minValue);
#endif
#ifdef ACC_SEMILAG_PPM
            Vec a[3];
            compute_ppm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h4, k + WID, a, minValue);
#endif
#ifdef ACC_SEMILAG_PQM
            Vec a[5];
            compute_pqm_coeff(dev_values + dev_columns[column].valuesOffset  + i_pcolumnv_cuda(j, 0, -1, nblocks), h8, k + WID, a, minValue);
#endif

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
            //limits in lagrangian k for target column. Also take into
            //account limits of target column
            int minGk = max(int(lagrangian_gk_l[minGkIndex]), int(dev_columns[column].minBlockK * WID));
            int maxGk = min(int(lagrangian_gk_r[maxGkIndex]), int((dev_columns[column].maxBlockK + 1) * WID - 1));
            // Run along the column and perform the polynomial reconstruction
            for(int gk = dev_columns[column].minBlockK * WID; gk <= dev_columns[column].maxBlockK * WID; gk++)
            {
               if(gk < minGk || gk > maxGk)
               { continue; }

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
                #ifdef ACC_SEMILAG_PLM
                  target_density_r = v_norm_r * ( a[0] + v_norm_r * a[1] );
                  //printf("ACC_SEMILAG_PLM\n");
                #endif
                #ifdef ACC_SEMILAG_PPM
                  target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * a[2] ) );
                  //printf("ACC_SEMILAG_PPM\n");
                #endif
                #ifdef ACC_SEMILAG_PQM
                  target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * ( a[2] + v_norm_r * ( a[3] + v_norm_r * a[4] ) ) ) );
                #endif

                //store values, one element at a time. All blocks have been created by now.
                const Vec target_density = target_density_r - target_density_l;
                for (int target_i=0; target_i < VECL; ++target_i)
                {

                  // do the conversion from Realv to Realf here, faster than doing it in accumulation
                  const Realf tval = target_density[target_i];
                  const uint tcell = target_cell[target_i];

                  (&dev_blockData[dev_columns[column].targetBlockOffsets[blockK]])[tcell] += tval;

                }  // for-loop over vector elements
            } // for loop over target k-indices of current source block
         } // for-loop over source blocks
      } //for loop over j index
   } //for loop over columns
#endif // NOT CUDA_REALF
}

Realf* acceleration_1_wrapper
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
)
{
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

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

#ifdef CUDA_REALF
  int blocks = BLOCKS;
  if (THREADS != VECL) printf("CUDA ERROR! VECL does not match thread count.\n");
#else
  int blocks = (totalColumns / THREADS) + 1;
#endif

  acceleration_1<<<blocks, THREADS>>>
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
        bdsw3
  );
  cudaDeviceSynchronize();
  HANDLE_ERROR( cudaMemcpy(blockData, dev_blockData, bdsw3*sizeof(double), cudaMemcpyDeviceToHost) );

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTime;
  cudaEventElapsedTime(&elapsedTime, start, stop);
  //printf("%.3f ms\n", elapsedTime);

  HANDLE_ERROR( cudaFree(dev_blockData) );
  HANDLE_ERROR( cudaFree(dev_cell_indices_to_id) );
  HANDLE_ERROR( cudaFree(dev_columns) );
  HANDLE_ERROR( cudaFree(dev_values) );

  return blockData;
}
