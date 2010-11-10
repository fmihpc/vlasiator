#ifndef CUDA_COMMON_CU
#define CUDA_COMMON_CU

#include <cstdlib>
#include <iostream>
                                   
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "logger.h"
#include "common.h"
#include "devicegrid.h"
#include "grid.h"
#include "parameters.h"

using namespace std;
                                               
extern Logger logger;
                                               
texture<Real,1,cudaReadModeElementType> texRef_avgs;
texture<Real,1,cudaReadModeElementType> texRef_cellParams;
texture<Real,1,cudaReadModeElementType> texRef_d1x;
texture<Real,1,cudaReadModeElementType> texRef_d1y;
texture<Real,1,cudaReadModeElementType> texRef_d1z;
texture<Real,1,cudaReadModeElementType> texRef_bparms;

texture<uint,1,cudaReadModeElementType> texRef_nbrs;

inline __device__ uint bindex2(uint bix,uint biy) {return biy*gridDim.x+bix;}
inline __device__ uint tindex2(uint tix,uint tiy) {return tiy*WID+tix;}
inline __device__ uint tindex3(uint tix,uint tiy,uint tiz) {return tiz*WID2+tiy*WID+tix;}

inline __device__ void loadVelNbrs(uint MYBLOCK,uint* nbrs,uint* sha_nbr) {
   sha_nbr[threadIdx.y*WID+threadIdx.x] = nbrs[MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x];
}

inline __device__ void transpose_yz_1warp(Real* sha_arr) {
   cuint myz1 = (threadIdx.y+threadIdx.z  )%WID;
   cuint myz2 = (threadIdx.y+threadIdx.z+2)%WID;
   creal tmp1 = sha_arr[myz1*WID2 + threadIdx.y*WID + threadIdx.x];
   creal tmp2 = sha_arr[myz2*WID2 + threadIdx.y*WID + threadIdx.x];
   sha_arr[threadIdx.y*WID2 + myz1*WID + threadIdx.x] = tmp1;
   sha_arr[threadIdx.y*WID2 + myz2*WID + threadIdx.x] = tmp2;
}

/** Transpose given array in yz. This function is called with 
 * 2 warps (=64 threads).
 * This function does not cause shared memory bank conflicts.
 * @param sha_arr Pointer to the array to transpose.
 */
inline __device__ void transpose_yz_2warp(Real* sha_arr) {
   cuint myz = (threadIdx.y+threadIdx.z)%4;
   creal tmp = sha_arr[myz*WID2 + threadIdx.y*WID + threadIdx.x];
   __syncthreads();
   sha_arr[threadIdx.y*WID2 + myz*WID + threadIdx.x] = tmp;
}

inline __device__ void transpose_xz_1warp(Real* sha_arr) {
   cuint myk1 = (threadIdx.z+threadIdx.x  )%WID;
   cuint myk2 = (threadIdx.z+threadIdx.x+2)%WID;
   creal tmp1 = sha_arr[myk1*WID2 + threadIdx.y*WID + threadIdx.x];
   creal tmp2 = sha_arr[myk2*WID2 + threadIdx.y*WID + threadIdx.x];
   sha_arr[threadIdx.x*WID2 + threadIdx.y*WID + myk1] = tmp1;
   sha_arr[threadIdx.x*WID2 + threadIdx.y*WID + myk2] = tmp2;
}

/** Transpose given 3D array in xz. This function is 
 * called with 2 warps (=64 threads).
 * @param sha_arr Pointer to the array to transpose.
 */
inline __device__ void transpose_xz_2warp(Real* sha_arr) {
   cuint myk = (threadIdx.z+threadIdx.x)%4;
   creal tmp = sha_arr[myk*WID2 + threadIdx.y*WID + threadIdx.x];
   __syncthreads();
   sha_arr[threadIdx.x*WID2 + threadIdx.y*WID + myk] = tmp;
}

inline __device__ Real reconstruct_neg(Real avg,Real d1,Real d2) {return avg - d2/24.0f + 0.5f*d1 + d2/8.0f;}
inline __device__ Real reconstruct_pos(Real avg,Real d1,Real d2) {return avg - d2/24.0f - 0.5f*d1 + d2/8.0f;}

inline __device__ Real minmod(Real a,Real b) { // MinMod limiter
   if (a*b < 0.0f) return 0.0f;
   if (fabsf(a) < fabsf(b)) return a;
   else return b;
}

inline __device__ Real MClimiter(Real xl1,Real xcc,Real xr1) { // MC limiter
   creal forw = xr1-xcc;
   creal back = xcc-xl1;
   creal cntr = 0.5f*(xr1-xl1);
   creal slope1 = minmod(2.0f*forw,cntr);
   creal slope2 = minmod(2.0f*back,cntr);
   return minmod(slope1,slope2);
}
      
inline __device__ Real superbee(Real xl1,Real xcc,Real xr1) { // SuperBee limiter
   creal forw = xr1-xcc;
   creal back = xcc-xl1;
   
   Real tmp = fminf(fabsf(back),fabsf(forw));
   tmp = fminf(tmp,0.5f*fabsf(back));
   tmp = fminf(tmp,fabsf(forw));
   
   Real sgnfwd = 1.0f;
   if (forw < 0.0f) sgnfwd = -1.0f;
   Real sgnbck = 1.0f;
   if (back < 0.0f) sgnbck = -1.0f;
   
   return (sgnfwd+sgnbck)*tmp;
}

inline __device__ Real vanLeer(Real xl1,Real xcc,Real xr1) { // Van Leer limiter
   creal EPS = 1.0e-20;
   //creal EPS = 0.333333f*(xl1+xcc+xr1)*1.0e-4f;
   //if (fabsf(xr1-xl1) < EPS) return 0.0f;
   //return fmaxf((xr1-xcc)*(xcc-xl1),0.0f)/(xr1-xl1);
   return fmaxf((xr1-xcc)*(xcc-xl1),0.0f)/(EPS+xr1-xl1);
}
/*
inline __device__ void kernel_stencils(uint MYIND,Real* sha_avg,Real* d1,Real* d2) {
   // Read a five-point stencil into registers:
   Real xl2 = sha_avg[MYIND + 2*WID2 - 2*WID2];
   Real xl1 = sha_avg[MYIND + 2*WID2 -   WID2];
   Real xcc = sha_avg[MYIND + 2*WID2         ];
   Real xr1 = sha_avg[MYIND + 2*WID2 +   WID2];
   Real xr2 = sha_avg[MYIND + 2*WID2 + 2*WID2];
   
   // Calculate 1st and 2nd derivatives:
   //d1[MYIND] = vanLeer(xl1,xcc,xr1);
   d1[MYIND] = 0.0f;
   d2[MYIND] = 0.0f;
}
*/
#endif
