#include "cuda_header.cuh"
#include "open_acc_map_h.cuh"
#include "../vlasovsolver/vec.h"
#include "../definitions.h"
//#include "../vlasovsolver/cpu_face_estimates.hpp"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum face_estimate_order {h4, h5, h6, h8};

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
__device__ Vec slope_limiter_sb(const Vec &l, const Vec &m, const Vec &r)
{
  Vec a = r-m;
  Vec b = m-l;
  const Vec slope1 = minmod(a, 2*b);
  const Vec slope2 = minmod(2*a, b);
  return maxmod(slope1, slope2);
}
__device__ Vec slope_limiter(const Vec &l, const Vec &m, const Vec &r)
{
   return slope_limiter_sb(l,m,r);
}
__device__ void slope_limiter(const Vec& l,const Vec& m, const Vec& r, Vec& slope_abs, Vec& slope_sign)
{
   const Vec slope = slope_limiter(l,m,r);
   slope_abs = abs(slope);
   slope_sign = select(slope > 0, Vec(1.0), Vec(-1.0));
}
__device__ void compute_h3_left_face_derivative(const Vec * const values, uint k, Vec &fv_l)
{
  /*compute left value*/
  fv_l = 1.0/12.0 * (15 * (values[k] - values[k - 1]) - (values[k + 1] - values[k - 2]));
}
__device__ void compute_h4_left_face_value(const Vec * const values, uint k, Vec &fv_l)
{
  //compute left value
  fv_l = 1.0/12.0 * ( - 1.0 * values[k - 2]
                      + 7.0 * values[k - 1]
                      + 7.0 * values[k]
                      - 1.0 * values[k + 1]);
}
__device__ void compute_h4_left_face_derivative(const Vec * const values, uint k, Vec &fd_l)
{
  fd_l = 1.0/12.0 * (15.0 * (values[k] - values[k - 1]) - (values[k + 1] - values[k - 2]));
}
__device__ void compute_h5_face_values(const Vec * const values, uint k, Vec &fv_l, Vec &fv_r)
{
  //compute left values
  fv_l = 1.0/60.0 * (- 3.0 * values[k - 2]
                     + 27.0 * values[k - 1]
                     + 47.0 * values[k ]
                     - 13.0 * values[k + 1]
                     + 2.0 * values[k + 2]);
  fv_r = 1.0/60.0 * ( 2.0 * values[k - 2]
                     - 13.0 * values[k - 1]
                     + 47.0 * values[k]
                     + 27.0 * values[k + 1]
                     - 3.0 * values[k + 2]);
}
__device__ void compute_h5_left_face_derivative(const Vec * const values, uint k, Vec &fd_l)
{
  fd_l = 1.0/180.0 * (245 * (values[k] - values[k - 1])
                     - 25 * (values[k + 1] - values[k - 2])
                     + 2 * (values[k + 2] - values[k - 3]));
}
__device__ void compute_h6_left_face_value(const Vec * const values, uint k, Vec &fv_l)
{
  //compute left value
   fv_l = 1.0/60.0 * (values[k - 3]
            - 8.0 * values[k - 2]
            + 37.0 * values[k - 1]
            + 37.0 * values[k ]
            - 8.0 * values[k + 1]
            + values[k + 2]);
}
__device__ void compute_h7_left_face_derivative(const Vec * const values, uint k, Vec &fd_l){
    fd_l = 1.0/5040.0 * (
       + 9.0 * values[k - 4]
       - 119.0 * values[k - 3]
       + 889.0 * values[k - 2]
       - 7175.0 * values[k - 1]
       + 7175.0 * values[k]
       - 889.0 * values[k + 1]
       + 119.0 * values[k + 2]
       - 9.0 * values[k + 3]);
}
__device__ void compute_h8_left_face_value(const Vec * const values, uint k, Vec &fv_l)
{
   fv_l = 1.0/840.0 * (
              - 3.0 * values[k - 4]
              + 29.0 * values[k - 3]
              - 139.0 * values[k - 2]
              + 533.0 * values[k - 1]
              + 533.0 * values[k]
              - 139.0 * values[k + 1]
              + 29.0 * values[k + 2]
              - 3.0 * values[k + 3]);
}

__device__ void compute_filtered_face_values(const Vec * const values, uint k, face_estimate_order order, Vec &fv_l, Vec &fv_r, const Realv threshold)
{
  switch(order)
  {
      case h4:
         compute_h4_left_face_value(values, k, fv_l);
         compute_h4_left_face_value(values, k + 1, fv_r);
         break;
      case h5:
         compute_h5_face_values(values, k, fv_l, fv_r);
         break;
      default:
      case h6:
          compute_h6_left_face_value(values, k, fv_l);
          compute_h6_left_face_value(values, k + 1, fv_r);
         break;
      case h8:
          compute_h8_left_face_value(values, k, fv_l);
          compute_h8_left_face_value(values, k + 1, fv_r);
         break;
  }
  Vec slope_abs, slope_sign;
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realv scale = 1./threshold;
  slope_limiter(values[k -1]*scale, values[k]*scale, values[k + 1]*scale, slope_abs, slope_sign);
  slope_abs = slope_abs*threshold;

  //check for extrema, flatten if it is
  Vecb is_extrema = (slope_abs == Vec(0.0));
  if(horizontal_or(is_extrema))
  {
     fv_r = select(is_extrema, values[k], fv_r);
     fv_l = select(is_extrema, values[k], fv_l);
  }
  //Fix left face if needed; boundary value is not bounded
  Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
  if(horizontal_or (filter))
  {
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     fv_l = select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
  }
  //Fix  face if needed; boundary value is not bounded
  filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
  if(horizontal_or (filter))
  {
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     fv_r = select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
  }
}

__device__ void compute_filtered_face_values_derivatives(const Vec * const values,uint k, face_estimate_order order, Vec &fv_l, Vec &fv_r, Vec &fd_l, Vec &fd_r, const Realv threshold)
{
   switch(order)
   {
       case h4:
          compute_h4_left_face_value(values, k, fv_l);
          compute_h4_left_face_value(values, k + 1, fv_r);
          compute_h3_left_face_derivative(values, k, fd_l);
          compute_h3_left_face_derivative(values, k + 1, fd_r);
          break;
       case h5:
          compute_h5_face_values(values, k, fv_l, fv_r);
          compute_h4_left_face_derivative(values, k, fd_l);
          compute_h4_left_face_derivative(values, k + 1, fd_r);
          break;
       default:
       case h6:
          compute_h6_left_face_value(values, k, fv_l);
          compute_h6_left_face_value(values, k + 1, fv_r);
          compute_h5_left_face_derivative(values, k, fd_l);
          compute_h5_left_face_derivative(values, k + 1, fd_r);
          break;
       case h8:
          compute_h8_left_face_value(values, k, fv_l);
          compute_h8_left_face_value(values, k + 1, fv_r);
          compute_h7_left_face_derivative(values, k, fd_l);
          compute_h7_left_face_derivative(values, k + 1, fd_r);
          break;
   }
   Vec slope_abs,slope_sign;
   // scale values closer to 1 for more accurate slope limiter calculation
   const Realv scale = 1./threshold;
   slope_limiter(values[k -1]*scale, values[k]*scale, values[k + 1]*scale, slope_abs, slope_sign);
   slope_abs = slope_abs*threshold;
   //check for extrema, flatten if it is
   Vecb is_extrema = (slope_abs == Vec(0.0));
   if(horizontal_or(is_extrema))
   {
      fv_r = select(is_extrema, values[k], fv_r);
      fv_l = select(is_extrema, values[k], fv_l);
      fd_l = select(is_extrema, 0.0 , fd_l);
      fd_r = select(is_extrema, 0.0 , fd_r);
   }
   //Fix left face if needed; boundary value is not bounded or slope is not consistent
   Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
   if(horizontal_or (filter))
   {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
      fd_l=select(filter, slope_sign * slope_abs, fd_l);
   }
   //Fix right face if needed; boundary value is not bounded or slope is not consistent
   filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
   if(horizontal_or (filter))
   {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
      fd_r=select(filter, slope_sign * slope_abs, fd_r);
   }
}
__device__ void filter_pqm_monotonicity(Vec *values, uint k, Vec &fv_l, Vec &fv_r, Vec &fd_l, Vec &fd_r){
   const Vec root_outside = Vec(100.0); //fixed values give to roots clearly outside [0,1], or nonexisting ones*/
   /*second derivative coefficients, eq 23 in white et al.*/
   Vec b0 =   60.0 * values[k] - 24.0 * fv_r - 36.0 * fv_l + 3.0 * (fd_r - 3.0 * fd_l);
   Vec b1 = -360.0 * values[k] + 36.0 * fd_l - 24.0 * fd_r + 168.0 * fv_r + 192.0 * fv_l;
   Vec b2 =  360.0 * values[k] + 30.0 * (fd_r - fd_l) - 180.0 * (fv_l + fv_r);
   /*let's compute sqrt value to be used for computing roots. If we
    take sqrt of negaitve numbers, then we instead set a value that
    will make the root to be +-100 which is well outside range
    of[0,1]. We do not catch FP exceptions, so sqrt(negative) are okish (add
    a max(val_to_sqrt,0) if not*/
   const Vec val_to_sqrt = b1 * b1 - 4 * b0 * b2;
/*
#ifdef VEC16F_AGNER
   //this sqrt gives 10% more perf on acceleration on KNL. Also fairly
   //accurate with AVX512ER. On Xeon it is not any faster, and less accurate.
   const Vec sqrt_val = select(val_to_sqrt < 0.0,
                               b1 + 200.0 * b2,
                               val_to_sqrt * approx_rsqrt(val_to_sqrt));
#else
*/
   const Vec sqrt_val = select(val_to_sqrt < 0.0,
                               b1 + 200.0 * b2,
                               sqrt(val_to_sqrt));
//#endif
   //compute roots. Division is safe with vectorclass (=inf)
   const Vec root1 = (-b1 + sqrt_val) / (2 * b2);
   const Vec root2 = (-b1 - sqrt_val) / (2 * b2);

   /*PLM slope, MC limiter*/
   Vec plm_slope_l = 2.0 * (values[k] - values[k - 1]);
   Vec plm_slope_r = 2.0 * (values[k + 1] - values[k]);
   Vec slope_sign = plm_slope_l + plm_slope_r; //it also has some magnitude, but we will only use its sign.
   /*first derivative coefficients*/
   const Vec c0 = fd_l;
   const Vec c1 = b0;
   const Vec c2 = b1 / 2.0;
   const Vec c3 = b2 / 3.0;
   //compute both slopes at inflexion points, at least one of these
   //is with [0..1]. If the root is not in this range, we
   //simplify later if statements by setting it to the plm slope
   //sign
   Vec root1_slope = select(root1 >= 0.0 && root1 <= 1.0,
                             c0  + root1 * ( c1 + root1 * (c2 + root1 * c3 ) ),
                             slope_sign);
   Vec root2_slope = select(root2 >= 0.0 && root2 <= 1.0,
                            c0  + root2 * ( c1 + root2 * (c2 + root2 * c3 ) ),
                            slope_sign);
   Vecb fixInflexion = root1_slope * slope_sign < 0.0 || root2_slope * slope_sign < 0.0;
   if (horizontal_or (fixInflexion) )
   {
      Realv valuesa[VECL];
      Realv fva_l[VECL];
      Realv fva_r[VECL];
      Realv fda_l[VECL];
      Realv fda_r[VECL];
      Realv slope_signa[VECL];
      values[k].store(valuesa);
      fv_l.store(fva_l);
      fd_l.store(fda_l);
      fv_r.store(fva_r);
      fd_r.store(fda_r);
      slope_sign.store(slope_signa);
      //todo store and then load data to avoid inserts (is it beneficial...?)
      //serialized the handling of inflexion points, these do not happen for smooth regions
      for(uint i = 0;i < VECL; i++)
      {
         if(fixInflexion[i])
         {
            //need to collapse, at least one inflexion point has wrong
            //sign.
            if(fabs(plm_slope_l[i]) <= fabs(plm_slope_r[i]))
            {
               //collapse to left edge (eq 21)
               fda_l[i] =  1.0 / 3.0 * ( 10 * valuesa[i] - 2.0 * fva_r[i] - 8.0 * fva_l[i]);
               fda_r[i] =  -10.0 * valuesa[i] + 6.0 * fva_r[i] + 4.0 * fva_l[i];
               //check if PLM slope is consistent (eq 28 & 29)
               if (slope_signa[i] * fda_l[i] < 0)
               {
                  fda_l[i] =  0;
                  fva_r[i] =  5 * valuesa[i] - 4 * fva_l[i];
                  fda_r[i] =  20 * (valuesa[i] - fva_l[i]);
               }
               else if (slope_signa[i] * fda_r[i] < 0)
               {
                  fda_r[i] =  0;
                  fva_l[i] =  0.5 * (5 * valuesa[i] - 3 * fva_r[i]);
                  fda_l[i] =  10.0 / 3.0 * (-valuesa[i] + fva_r[i]);
               }
            }
            else
            {
               //collapse to right edge (eq 21)
               fda_l[i] =  10.0 * valuesa[i] - 6.0 * fva_l[i] - 4.0 * fva_r[i];
               fda_r[i] =  1.0 / 3.0 * ( - 10.0 * valuesa[i] + 2 * fva_l[i] + 8 * fva_r[i]);
               //check if PLM slope is consistent (eq 28 & 29)
               if (slope_signa[i] * fda_l[i] < 0)
               {
                  fda_l[i] =  0;
                  fva_r[i] =  0.5 * ( 5 * valuesa[i] - 3 * fva_l[i]);
                  fda_r[i] =  10.0 / 3.0 * (valuesa[i] - fva_l[i]);
               }
               else if (slope_signa[i] * fda_r[i] < 0)
               {
                  fda_r[i] =  0;
                  fva_l[i] =  5 * valuesa[i] - 4 * fva_r[i];
                  fda_l[i] =  20.0 * ( - valuesa[i] + fva_r[i]);
               }
            }
         }
      }
      fv_l.load(fva_l);
      fd_l.load(fda_l);
      fv_r.load(fva_r);
      fd_r.load(fda_r);
   }
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

__device__ void compute_ppm_coeff(const Vec * const values, face_estimate_order order, uint k, Vec a[3], const Realv threshold)
{
  Vec fv_l; //left face value
  Vec fv_r; //right face value
  compute_filtered_face_values(values, k, order, fv_l, fv_r, threshold);
  //Coella et al, check for monotonicity
  const Vec one_sixth(1.0/6.0);
  Vec m_face = fv_l;
  Vec p_face = fv_r;
  m_face = select((p_face - m_face) * (values[k] - 0.5 * (m_face + p_face)) >
                  (p_face - m_face) * (p_face - m_face) * one_sixth,
                  3 * values[k] - 2 * p_face, m_face);
  p_face = select(-(p_face - m_face) * (p_face - m_face) * one_sixth >
                  (p_face - m_face) * (values[k] - 0.5 * (m_face + p_face)),
                  3 * values[k] - 2 * m_face, p_face);
  //Fit a second order polynomial for reconstruction see, e.g., White
  //2008 (PQM article) (note additional integration factors built in,
  //contrary to White (2008) eq. 4
  a[0] = m_face;
  a[1] = 3.0 * values[k] - 2.0 * m_face - p_face;
  a[2] = (m_face + p_face - 2.0 * values[k]);
}

__device__ void compute_pqm_coeff(Vec *values, face_estimate_order order, uint k, Vec a[5], const Realv threshold)
{
   Vec fv_l; /*left face value*/
   Vec fv_r; /*right face value*/
   Vec fd_l; /*left face derivative*/
   Vec fd_r; /*right face derivative*/
   compute_filtered_face_values_derivatives(values, k, order, fv_l, fv_r, fd_l, fd_r, threshold);
   filter_pqm_monotonicity(values, k, fv_l, fv_r, fd_l, fd_r);
   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0] = fv_l;
   a[1] = fd_l/2.0;
   a[2] =  10.0 * values[k] - 4.0 * fv_r - 6.0 * fv_l + 0.5 * (fd_r - 3 * fd_l);
   a[3] = -15.0 * values[k]  + 1.5 * fd_l - fd_r + 7.0 * fv_r + 8 * fv_l;
   a[4] =   6.0 * values[k] +  0.5 * (fd_r - fd_l) - 3.0 * (fv_l + fv_r);
}

__global__ void acceleration_1
(
  double *dev_blockData,
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
            /*
              if(column==0 && j==0 && k==0)
              for(int a_1 = 0; a_1 < VECL; a_1++)
              {
                  printf("CUDA a[0]: column = %d; j = %d; k = %d; a[%d] = %.12e\n", column, j, k, a_1, a[0][a_1]);
                  printf("CUDA a[1]: column = %d; j = %d; k = %d; a[%d] = %.12e\n", column, j, k, a_1, a[1][a_1]);
              }
            */
              #ifdef ACC_SEMILAG_PLM
                Vec a[2];
                compute_plm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), k + WID, a, minValue);
                //printf("ACC_SEMILAG_PLM\n");
              #endif
              #ifdef ACC_SEMILAG_PPM
                Vec a[3];
                compute_ppm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h4, k + WID, a, minValue);
                //printf("ACC_SEMILAG_PPM\n");
              #endif
              #ifdef ACC_SEMILAG_PQM
                Vec a[5];
                compute_pqm_coeff(dev_values + dev_columns[column].valuesOffset  + i_pcolumnv_cuda(j, 0, -1, nblocks), h8, k + WID, a, minValue);
                //printf("ACC_SEMILAG_PQM\n");
              #endif
              /*
              for(int a_1 = 0; a_1 < VECL; a_1++)
              {
                printf("CUDA: a[0][%d] = %d; a[1][%d] = %d; a[2][%d] = %d;\n", a[0][a_1], a[1][a_1], a[2][a_1]);
              }
              */
             // set the initial value for the integrand at the boundary at v = 0
             // (in reduced cell units), this will be shifted to target_density_1, see below.
             Vec target_density_r(0.0);
             /*
             if(column==0 && j==0 && k==0)
             {
               for(int a_1 = 0; a_1 < VECL; a_1++)
               {
                 printf("CUDA 0: target_density_r [%d] = %.2f\n", a_1, target_density_r[a_1]);
               }
             }
             */
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
                /*
                if(column==0 && j==0 && k==0)
                {
                  printf("CUDA 1: minGk = %d; gk = %d; maxGk = %d;\n", minGk, gk, maxGk);
                }
                */
                if(gk < minGk || gk > maxGk)
                { continue; }
                //printf("CUDA 6\n");
                const int blockK = gk/WID;
                const int gk_mod_WID = (gk - blockK * WID);
                /*
                if(column==0 && j==0 && k==0)
                {
                 printf("CUDA 2: blockK = %d; gk_mod_WID = %d;\n", blockK, gk_mod_WID);
                }
                */
                //the block of the Lagrangian cell to which we map
                //const int target_block(target_block_index_common + blockK * block_indices_to_id[2]);
                //cell indices in the target block  (TODO: to be replaced by
                //compile time generated scatter write operation)
                const Veci target_cell(target_cell_index_common + gk_mod_WID * dev_cell_indices_to_id[2]);
                /*
                if(column==0 && j==0 && k==0)
                {
                  for(int a_1 = 0; a_1 < VECL; a_1++)
                  {
                    printf("CUDA 3: target_cell [%d] = %d\n", a_1, target_cell[a_1]);
                  }
                }
                */
                //the velocity between which we will integrate to put mass
                //in the targe cell. If both v_r and v_l are in same cell
                //then v_1,v_2 should be between v_l and v_r.
                //v_1 and v_2 normalized to be between 0 and 1 in the cell.
                //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
                const Vec v_norm_r = (  min(  max( (gk + 1) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;
                /*
                if(column==0 && j==0 && k==0)
                {
                  printf("CUDA CHECK 1: gk = %d; intersection_dk = %d; intersection_min = %d; i_dv = %d;\n", gk, intersection_dk, intersection_min, i_dv);
                  for(int a_1 = 0; a_1 < VECL; a_1++)
                  {
                    printf("CUDA CHECK 2: v_l [%d] = %.2f\n", a_1, v_l[a_1]);
                    printf("CUDA CHECK 3: v_r [%d] = %.2f\n", a_1, v_r[a_1]);
                  }
                }
                */
                /*
                if(column==0 && j==0 && k==0)
                {
                  for(int a_1 = 0; a_1 < VECL; a_1++)
                  {
                    printf("CUDA 4: v_norm_r [%d] = %.2f\n", a_1, v_norm_r[a_1]);
                  }
                }
                */
                /*shift, old right is new left*/
                const Vec target_density_l = target_density_r;
                /*
                if(column==0 && j==0 && k==0)
                {
                  for(int a_1 = 0; a_1 < VECL; a_1++)
                  {
                    printf("CUDA 5: target_density_l [%d] = %.2f\n", a_1, target_density_l[a_1]);
                  }
                }
                */
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
                  //printf("CUDA 8\n");
                  // do the conversion from Realv to Realf here, faster than doing it in accumulation
                  const Realf tval = target_density[target_i];
                  const uint tcell = target_cell[target_i];
                  //if(target_i == ( (VECL) - 1) )
                  //  printf("CUDA: tval = %.2f; tcell = %d;\n", tval, tcell);
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
        bdsw3
  );

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTime;
  cudaEventElapsedTime (&elapsedTime, start, stop);
  //printf("%.3f\n", elapsedTime);
  cudaDeviceSynchronize();
  HANDLE_ERROR( cudaMemcpy(blockData, dev_blockData, bdsw3*sizeof(double), cudaMemcpyDeviceToHost) );

  HANDLE_ERROR( cudaFree(dev_blockData) );
  HANDLE_ERROR( cudaFree(dev_cell_indices_to_id) );
  HANDLE_ERROR( cudaFree(dev_columns) );
  HANDLE_ERROR( cudaFree(dev_values) );

  return blockData;
}
