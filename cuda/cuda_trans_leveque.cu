#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <limits>
#include "../common.h"

texture<uint,1,cudaReadModeElementType> texRef_nbrsSpa;

const uint I0_J0 = 0;
const uint I1_J0 = 1;
const uint I2_J0 = 2;
const uint I0_J1 = 3;
const uint I1_J1 = 4;
const uint I2_J1 = 5;
const uint I0_J2 = 6;
const uint I1_J2 = 7;
const uint I2_J2 = 8;

__device__ void dropMutex(int* flag) {
   atomicExch(flag,0);
}

__device__ void raiseMutex(int* flag) {
   while (atomicExch(flag,1) == 1) { }
}

__device__ Real limiter(Real theta_up,Real theta_lo,Real xcc) {
   return superbee(theta_up/theta_lo);
}

__device__ void warpReduce(volatile Real* array,uint MYIND) {
   array[MYIND] += array[MYIND+32];
   array[MYIND] += array[MYIND+16];
   array[MYIND] += array[MYIND+8 ];
   array[MYIND] += array[MYIND+4 ];
   array[MYIND] += array[MYIND+2 ];
   array[MYIND] += array[MYIND+1 ];
}

__device__ void blockVelocityMoments(uint MYBLOCK,uint MYIND,Real* sha_avgs,Real* blockParams,Real* gpuRho,Real* gpuRhovx,Real* gpuRhovy,Real* gpuRhovz) {
   // Reduce contributions to velocity momemts of this block to 
   // given global arrays:
   __shared__ volatile Real velMoment[WID3];
   velMoment[MYIND] = sha_avgs[MYIND];
   __syncthreads();

   if (MYIND < 32) warpReduce(velMoment,MYIND);
   if (MYIND == 0) gpuRho[MYBLOCK] = velMoment[0]*blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   __syncthreads();
   
   velMoment[MYIND] = sha_avgs[MYIND]*(blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX]);
   __syncthreads();
   if (MYIND < 32) warpReduce(velMoment,MYIND);
   if (MYIND == 0) gpuRhovx[MYBLOCK] = velMoment[0]*blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   __syncthreads();
   
   velMoment[MYIND] = sha_avgs[MYIND]*(blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY]);
   __syncthreads();
   if (MYIND < 32) warpReduce(velMoment,MYIND);
   if (MYIND == 0) gpuRhovy[MYBLOCK] = velMoment[0]*blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   __syncthreads();
   
   velMoment[MYIND] = sha_avgs[MYIND]*(blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ]);
   __syncthreads();
   if (MYIND < 32) warpReduce(velMoment,MYIND);
   if (MYIND == 0) gpuRhovz[MYBLOCK] = velMoment[0]*blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];

}

__device__ Real fetchAvgs(uint arrayPosition) {
   // 2D Texture fetches are required for large arrays:
   cuint YCRD = arrayPosition / CUDA_WIDTH;
   cuint XCRD = arrayPosition - YCRD*CUDA_WIDTH;
   return tex2D(texRef_avgs2D,XCRD+HALF,YCRD+HALF);
}

__global__ void cuda_trans_fluxes(Real* avgs,Real* flux,uint offsetNbrs,uint offsetCellParams,uint offsetBlockParams,int* semaphores,Real dt) {
   __shared__ Real dfdt[9*WID3];
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
  
   uint nbrBlock;
   Real xm2,xm1,xcc;
   Real Vx,Vy,dt_per_dx,dt_per_dy,R,corr_wave;
   
   cuint MYIND    = index(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYOFFSET = blockIdx.y*gridDim.x + blockIdx.x;
   cuint MYBLOCK  = tex1Dfetch(texRef_nbrsSpa,offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 13);

   // Fetch block parameters to shared mem:
   const int i = MYIND % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYOFFSET*SIZE_BLOCKPARAMS+i);
   
   // Clear dfdt array:
   for (int i=0; i<9; ++i) dfdt[i*WID3 + MYIND] = ZERO;
   
   // ***** Consider interface between cells (i,j,k) and (i-1,j,k) *****
   dt_per_dx = dt/tex1Dfetch(texRef_cellParams,offsetCellParams + CellParams::DX);
   dt_per_dy = dt/tex1Dfetch(texRef_cellParams,offsetCellParams + CellParams::DY);

   // Case Vx > 0:
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 13);
   xcc = fetchAvgs(nbrBlock*WID3 + MYIND); 
   
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 12);
   xm1 = fetchAvgs(nbrBlock*WID3 + MYIND);
   
   Vx = fmax(ZERO,blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX]);
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];

   R = xcc - xm1;
   
   // Increment waves:
   dfdt[I1_J1*WID3 + MYIND] += Vx*xm1*dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] -= Vx*xm1*dt_per_dx;
   
   // Transverse Increment waves:
   dfdt[I1_J2*WID3 + MYIND] -= HALF*dt_per_dx*Vx*fmax(ZERO,Vy)*R * dt_per_dy; // Vy > 0
   dfdt[I1_J1*WID3 + MYIND] += HALF*dt_per_dx*Vx*fmax(ZERO,Vy)*R * dt_per_dy;
   dfdt[I1_J0*WID3 + MYIND] += HALF*dt_per_dx*Vx*fmin(ZERO,Vy)*R * dt_per_dy; // Vy < 0
   dfdt[I1_J1*WID3 + MYIND] -= HALF*dt_per_dx*Vx*fmin(ZERO,Vy)*R * dt_per_dy;   
   
   // Correction waves:
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 27);
   xm2 = fetchAvgs(nbrBlock*WID3 + MYIND);
   R *= limiter(xm1-xm2,R+EPSILON,xcc);
   corr_wave = HALF*Vx*(ONE - dt_per_dx*Vx)*R;
   dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] -= corr_wave * dt_per_dx;
      
   // Transverse Correction waves:
   dfdt[I1_J2*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy > 0
   dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J2*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J1*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;

   dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   
   // Case Vx < 0:
   R = xcc - xm1;
   Vx = fmin(ZERO,blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX]);
   
   // Increment waves:
   dfdt[I1_J1*WID3 + MYIND] += Vx*xcc*dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] -= Vx*xcc*dt_per_dx;

   // Transverse Increment Waves:
   dfdt[I0_J2*WID3 + MYIND] -= HALF*dt_per_dx*Vx*fmax(ZERO,Vy)*R * dt_per_dy; // Vx > 0
   dfdt[I0_J1*WID3 + MYIND] += HALF*dt_per_dx*Vx*fmax(ZERO,Vy)*R * dt_per_dy;
   dfdt[I0_J0*WID3 + MYIND] += HALF*dt_per_dx*Vx*fmin(ZERO,Vy)*R * dt_per_dy; // Vx < 0
   dfdt[I0_J1*WID3 + MYIND] -= HALF*dt_per_dx*Vx*fmin(ZERO,Vy)*R * dt_per_dy;

   // Correction Waves:
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 14); // x+1 nbr
   xm2 = fetchAvgs(nbrBlock*WID3 + MYIND);
   //R *= limiter(xm1-xcc,R+EPSILON,xcc);
   R *= limiter(xm2-xcc,R+EPSILON,xcc);
   corr_wave = -HALF*Vx*(ONE + dt_per_dx*Vx)*R;
   dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] -= corr_wave * dt_per_dx;

   // Transverse Correction Waves:
   dfdt[I1_J2*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy > 0
   dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J2*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J1*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   
   dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy < 0
   dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
   dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;

   // ***** Consider Interface between cells (i,j,k) and (i,j-1,k) *****
   // Case Vy > 0:
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 10);
   xm1 = fetchAvgs(nbrBlock*WID3 + MYIND);
   
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 28);
   xm2 = fetchAvgs(nbrBlock*WID3 + MYIND);
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = fmax(ZERO,blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY]);
   
   R = xcc - xm1;
   
   // Increment Waves:
   dfdt[I1_J1*WID3 + MYIND] += Vy*xm1*dt_per_dy;
   dfdt[I1_J0*WID3 + MYIND] -= Vy*xm1*dt_per_dy;
   
   // Transverse Increment waves:
   dfdt[I2_J1*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx; // Vx > 0
   dfdt[I1_J1*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx; // Vx < 0
   dfdt[I1_J1*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx;

   // Correction Waves:
   R *= limiter(xm1-xm2,R+EPSILON,xcc);
   corr_wave = HALF*Vy*(ONE - dt_per_dy*Vy)*R;
   dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dy;
   dfdt[I1_J0*WID3 + MYIND] -= corr_wave * dt_per_dy;

   // Transverse Correction Waves:
   dfdt[I2_J1*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx > 0
   dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I2_J0*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I1_J0*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx < 0
   dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;

   // Case Vy < 0:
   Vy = fmin(ZERO,blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY]);
   
   R = xcc - xm1;
   
   // Increment Waves:
   dfdt[I1_J1*WID3 + MYIND] += Vy*xcc*dt_per_dy;
   dfdt[I1_J0*WID3 + MYIND] -= Vy*xcc*dt_per_dy;
   
   // Transverse Increment Waves:
   dfdt[I2_J0*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx; // Vx > 0
   dfdt[I1_J0*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx;
   dfdt[I0_J0*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx; // Vx < 0
   dfdt[I1_J0*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx;

   // Correction Waves:
   nbrBlock = tex1Dfetch(texRef_nbrsSpa, offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 16);
   xm2 = fetchAvgs(nbrBlock*WID3 + MYIND);
   R *= limiter(xm2-xcc,R+EPSILON,xcc);
   corr_wave = -HALF*Vy*(ONE + dt_per_dy*Vy)*R;
   dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dy;
   dfdt[I1_J0*WID3 + MYIND] -= corr_wave * dt_per_dy;

   // Transverse Correction Waves:
   dfdt[I2_J1*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx > 0
   dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I2_J0*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I1_J0*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx < 0
   dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
   dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;

   // Accumulate calculated time derivatives. Accumulates need to be 
   // atomic, as multiple thread blocks may be accumulating changes to the 
   // same block. Array semaphores contains an integer flag for each velocity block.
   uint STATUS = tex1Dfetch(texRef_nbrsSpa,offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 30);

   for (int i=0; i<9; ++i) {
      if (((STATUS >> (9+i)) & 1) == 0) continue; // If the neighbour does not exist, do not copy data
      
      nbrBlock = tex1Dfetch(texRef_nbrsSpa,offsetNbrs + MYOFFSET*SIZE_NBRS_SPA + 9 + i);
      if (MYIND == 0) raiseMutex(&(semaphores[nbrBlock]));
      __syncthreads();
      
      flux[nbrBlock*WID3+MYIND] += dfdt[i*WID3 + MYIND];
      if (MYIND == 0) dropMutex(&(semaphores[nbrBlock]));
   }

}

// Run this with 16 thread blocks of 32 threads (=1 warp)
__global__ void cuda_trans_reduce1(Real* array) {
   volatile __shared__ Real sha_arr[32];
   
   uint tid = threadIdx.x;
   uint index = blockIdx.x*blockDim.x*2 + tid;
   sha_arr[tid] = array[index] + array[index+blockDim.x];
   __syncthreads();
   
   // Unrolled reduction with 16 threads:
   if (tid < 16) {
      sha_arr[tid] += sha_arr[tid+16];
      sha_arr[tid] += sha_arr[tid+8];
      sha_arr[tid] += sha_arr[tid+4];
      sha_arr[tid] += sha_arr[tid+2];
      sha_arr[tid] += sha_arr[tid+1];
   }
   // Thread 0 writes out the result:
   if (tid == 0) array[index] = sha_arr[0];
}

// Second pass of reduction - run with a single thread block of 16 threads
__global__ void cuda_trans_reduce2(Real* array,Real* out) {
   volatile __shared__ Real sha_arr[32];
   
   uint tid = threadIdx.x;
   uint index = tid*64;
   sha_arr[tid] = array[index];
   sha_arr[tid+16] = ZERO;
   __syncthreads();
   
   // Unrolled reduction with 8 threads:
   if (tid < 16) {
      sha_arr[tid] += sha_arr[tid+16];
      sha_arr[tid] += sha_arr[tid+8];
      sha_arr[tid] += sha_arr[tid+4];
      sha_arr[tid] += sha_arr[tid+2];
      sha_arr[tid] += sha_arr[tid+1];
   }
   // Thread 0 writes out the result:
   if (tid == 0) out[0] = sha_arr[0];
}

__global__ void cuda_trans_propagate(Real* avgs,Real* flux,Real* cellParams,uint offsetBlockParams,uint offsetFluxBuffer,Real* fluxBuffer) {
   cuint MYIND   = index(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYBLOCK = blockIdx.y*gridDim.x + blockIdx.x;

   // Fetch old volume averages into a shared mem array:
   __shared__ Real sha_avgs[WID3];
   sha_avgs[MYIND] = avgs[MYBLOCK*WID3 + MYIND];
   
   // Propagate volume averages and store result:
   sha_avgs[MYIND] += flux[MYBLOCK*WID3 + MYIND];
   
   // Add contributions from remote hosts, if applicable:
   if (fluxBuffer != NULL) {
      sha_avgs[MYIND] += fluxBuffer[MYBLOCK*WID3 + MYIND];
   }
   avgs[MYBLOCK*WID3+MYIND] = sha_avgs[MYIND];
}

__global__ void cuda_trans_propagateWithMoments(Real* avgs,Real* flux,Real* cellParams,uint offsetBlockParams,
						Real* gpuRho,Real* gpuRhovx,Real* gpuRhovy,Real* gpuRhovz,uint offsetFluxBuffer,Real* fluxBuffer) {
   cuint MYIND   = index(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYBLOCK = blockIdx.y*gridDim.x + blockIdx.x;
   
   // Fetch block parameters into shared mem array (causes warp serialization):
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   const int i = MYIND % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);   

   // Fetch old volume averages into a shared mem array:
   __shared__ Real sha_avgs[WID3];
   sha_avgs[MYIND] = avgs[MYBLOCK*WID3 + MYIND];
   
   // Propagate volume averages:
   sha_avgs[MYIND]          += flux[MYBLOCK*WID3 + MYIND];
   if (fluxBuffer != NULL) {
      sha_avgs[MYIND] += fluxBuffer[MYBLOCK*WID3 + MYIND];
   }
   avgs[MYBLOCK*WID3+MYIND] = sha_avgs[MYIND];
   __syncthreads();
   
   blockVelocityMoments(MYBLOCK,MYIND,sha_avgs,blockParams,gpuRho,gpuRhovx,gpuRhovy,gpuRhovz);
}

__global__ void cuda_calcVelocityMoments(Real* avgs,Real* cellParams,uint offsetBlockParams,Real* gpuRho,Real* gpuRhovx,Real* gpuRhovy,Real* gpuRhovz) {
   cuint MYIND   = index(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYBLOCK = blockIdx.y*gridDim.x + blockIdx.x;
   
   // Fetch block parameters into shared mem array (causes warp serialization):
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   const int i = MYIND % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   // Fetch old volume averages into a shared mem array:
   __shared__ Real sha_avgs[WID3];
   sha_avgs[MYIND] = avgs[MYBLOCK*WID3 + MYIND];
   __syncthreads();

   blockVelocityMoments(MYBLOCK,MYIND,sha_avgs,blockParams,gpuRho,gpuRhovx,gpuRhovy,gpuRhovz);
}

__global__ void cuda_clearFluxes(Real* flux) {
   cuint MYIND   = index(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYBLOCK = blockIdx.y*gridDim.x + blockIdx.x;
   
   flux[MYBLOCK*WID3 + MYIND] = ZERO;
}

__global__ void cuda_clearVelocityMoments(Real* cellParams) {
   cuint MYCELL = blockIdx.y*gridDim.x + blockIdx.x;
   cuint MYIND = threadIdx.x + CellParams::RHO;
   cellParams[MYCELL*SIZE_CELLPARAMS + MYIND] = ZERO;
}
