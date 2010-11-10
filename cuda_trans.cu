#include "cuda_common.cu"
#include "cudalaunch.h"
#include "project.cu"

__host__ bool bindTexture(texture<uint,1,cudaReadModeElementType>& texRef,uint* arrptr,cuint& BYTES,size_t& offset) {
   bool success = true;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindUnsigned);
   cudaError_t cudaError = cudaBindTexture(&offset,&texRef,arrptr,&channelDesc,BYTES);
   if (cudaError != cudaSuccess) {
      logger << "(CUDA_TRANS): Failed to bind textures uint*!" << endl;
      success = false;
   }
   return success;
}

__host__ bool bindTexture(texture<Real,1,cudaReadModeElementType>& texRef,Real* arrptr,cuint& BYTES,size_t& offset) {
   bool success = true;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
   cudaError_t cudaError = cudaBindTexture(&offset,&texRef,arrptr,&channelDesc,BYTES);
   if (cudaError != cudaSuccess) {
      logger << "(CUDA_TRANS): Failed to bind textures to Real*!" << endl;
      success = false;
      }
   return success;
}

// ------------------ DEVICE FUNCTIONS -------------------
__device__ Real spatDerivs(Real xl1,Real xcc,Real xr1) {
   //return MClimiter(xl1,xcc,xr1);
   return superbee(xl1,xcc,xr1);
   //return vanLeer(xl1,xcc,xr1);
}

// ---------------------- KERNELS ------------------------
__global__ void calcSpatDerivs_1warp(uint OFFSET,Real* d1x,Real* d1y,Real* d1z,Real* d2x,Real* d2y,Real* d2z,uint* nbrs) {
   __shared__ uint sha_nbr[2*SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYIND2 = MYIND1 + 2*WID2;
   
   // Fetch neighbour indices:
   //if (threadIdx.z == 0) loadVelNbrs(MYBLOCK,nbrs,sha_nbr);
   //sha_nbr[threadIdx.y*WID+threadIdx.x] = nbrs[MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x];
   //sha_nbr[threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_nbr[MYIND1] = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages for this block:
   creal xcc_1 = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK + MYIND1);
   creal xcc_2 = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK + MYIND2);   
   
   // Fetch +/- x-neighbours and calculate derivatives:
   Real xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XNEG]*SIZE_VELBLOCK + MYIND1);
   Real xl1_2 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XNEG]*SIZE_VELBLOCK + MYIND2);
   Real xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XPOS]*SIZE_VELBLOCK + MYIND1);
   Real xr1_2 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XPOS]*SIZE_VELBLOCK + MYIND2);
   d1x[MYBLOCK*SIZE_DERIV + MYIND1] = spatDerivs(xl1_1,xcc_1,xr1_1);
   d1x[MYBLOCK*SIZE_DERIV + MYIND2] = spatDerivs(xl1_2,xcc_2,xr1_2);
   
   // Fetch +/- y-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YNEG]*SIZE_VELBLOCK + MYIND1); 
   xl1_2 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YNEG]*SIZE_VELBLOCK + MYIND2);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YPOS]*SIZE_VELBLOCK + MYIND1);
   xr1_2 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YPOS]*SIZE_VELBLOCK + MYIND2);
   d1y[MYBLOCK*SIZE_DERIV + MYIND1] = spatDerivs(xl1_1,xcc_1,xr1_1);
   d1y[MYBLOCK*SIZE_DERIV + MYIND2] = spatDerivs(xl1_2,xcc_2,xr1_2);
   
   // Fetch +/- z-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZNEG]*SIZE_VELBLOCK + MYIND1);
   xl1_2 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZNEG]*SIZE_VELBLOCK + MYIND2);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZPOS]*SIZE_VELBLOCK + MYIND1);
   xr1_2 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZPOS]*SIZE_VELBLOCK + MYIND2);
   d1z[MYBLOCK*SIZE_DERIV + MYIND1] = spatDerivs(xl1_1,xcc_1,xr1_1);
   d1z[MYBLOCK*SIZE_DERIV + MYIND2] = spatDerivs(xl1_2,xcc_2,xr1_2);
}

__global__ void calcSpatDerivs_2warp(uint OFFSET,Real* d1x,Real* d1y,Real* d1z,Real* d2x,Real* d2y,Real* d2z,uint* nbrs) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Fetch neighbour indices:
   //sha_nbr[MYIND] = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_nbr[threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages for this block:
   creal xcc_1 = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK + MYIND);
   
   // Fetch +/- x-neighbours and calculate derivatives:
   Real xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XNEG]*SIZE_VELBLOCK + MYIND);
   Real xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XPOS]*SIZE_VELBLOCK + MYIND);
   d1x[MYBLOCK*SIZE_DERIV + MYIND] = spatDerivs(xl1_1,xcc_1,xr1_1);

   // Fetch +/- y-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YNEG]*SIZE_VELBLOCK + MYIND);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YPOS]*SIZE_VELBLOCK + MYIND);
   d1y[MYBLOCK*SIZE_DERIV + MYIND] = spatDerivs(xl1_1,xcc_1,xr1_1);
   
   // Fetch +/- z-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZNEG]*SIZE_VELBLOCK + MYIND);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZPOS]*SIZE_VELBLOCK + MYIND);
   d1z[MYBLOCK*SIZE_DERIV + MYIND] = spatDerivs(xl1_1,xcc_1,xr1_1);
}

template<uint WARPS> __global__ void calcSpatDerivs_n2warp(uint OFFSET,Real* d1x,Real* d1y,Real* d1z,Real* d2x,Real* d2y,Real* d2z,uint* nbrs) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*WARPS*4];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/4)*gridDim.x + blockIdx.x;
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%4);
   
   // Fetch neighbour indices:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages for this block:
   creal xcc_1 = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK + MYIND);
    
   // Fetch +/- x-neighbours and calculate derivatives:
   Real xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::XNEG]*SIZE_VELBLOCK + MYIND);
   Real xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::XPOS]*SIZE_VELBLOCK + MYIND);
   d1x[MYBLOCK*SIZE_DERIV + MYIND] = spatDerivs(xl1_1,xcc_1,xr1_1);
   
   // Fetch +/- y-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::YNEG]*SIZE_VELBLOCK + MYIND);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::YPOS]*SIZE_VELBLOCK + MYIND);
   d1y[MYBLOCK*SIZE_DERIV + MYIND] = spatDerivs(xl1_1,xcc_1,xr1_1);
   
   // Fetch +/- z-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::ZNEG]*SIZE_VELBLOCK + MYIND);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::ZPOS]*SIZE_VELBLOCK + MYIND);
   d1z[MYBLOCK*SIZE_DERIV + MYIND] = spatDerivs(xl1_1,xcc_1,xr1_1);
}

template<uint WARPS> __global__ void calcSpatDerivs_nwarp(uint OFFSET,Real* d1x,Real* d1y,Real* d1z,Real* d2x,Real* d2y,Real* d2z,uint* nbrs) {
   //__shared__ uint sha_nbr[2*SIZE_NBRS_VEL*WARPS];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*WARPS];
   
   //cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/2)*gridDim.x + blockIdx.x;
   //cuint MYBLOCK = OFFSET + blockIdx.y*gridDim.x + blockIdx.x;
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%2);
   cuint MYIND2 = MYIND1 + 2*WID2;
   
   // Fetch neighbour indices:
   //sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL*2+MYIND1] = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] 
     = tex1Dfetch(texRef_nbrs, MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   __syncthreads();
   
   // Fetch volume averages for this block:
   creal xcc_1 = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK + MYIND1);
   creal xcc_2 = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK + MYIND2);
   
   // Fetch +/- x-neighbours and calculate derivatives:
   Real xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::XNEG]*SIZE_VELBLOCK + MYIND1);
   Real xl1_2 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::XNEG]*SIZE_VELBLOCK + MYIND2);
   Real xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::XPOS]*SIZE_VELBLOCK + MYIND1);
   Real xr1_2 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::XPOS]*SIZE_VELBLOCK + MYIND2);
   d1x[MYBLOCK*SIZE_DERIV + MYIND1] = spatDerivs(xl1_1,xcc_1,xr1_1);
   d1x[MYBLOCK*SIZE_DERIV + MYIND2] = spatDerivs(xl1_2,xcc_2,xr1_2);
   
   // Fetch +/- y-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::YNEG]*SIZE_VELBLOCK + MYIND1);
   xl1_2 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::YNEG]*SIZE_VELBLOCK + MYIND2);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::YPOS]*SIZE_VELBLOCK + MYIND1);
   xr1_2 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::YPOS]*SIZE_VELBLOCK + MYIND2);
   d1y[MYBLOCK*SIZE_DERIV + MYIND1] = spatDerivs(xl1_1,xcc_1,xr1_1);
   d1y[MYBLOCK*SIZE_DERIV + MYIND2] = spatDerivs(xl1_2,xcc_2,xr1_2);
   
   // Fetch +/- z-neighbours and calculate derivatives:
   xl1_1 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::ZNEG]*SIZE_VELBLOCK + MYIND1);
   xl1_2 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::ZNEG]*SIZE_VELBLOCK + MYIND2);
   xr1_1 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::ZPOS]*SIZE_VELBLOCK + MYIND1);
   xr1_2 = tex1Dfetch(texRef_avgs,sha_nbr[(threadIdx.z/2)*SIZE_NBRS_VEL + NbrsVel::ZPOS]*SIZE_VELBLOCK + MYIND2);
   d1z[MYBLOCK*SIZE_DERIV + MYIND1] = spatDerivs(xl1_1,xcc_1,xr1_1);
   d1z[MYBLOCK*SIZE_DERIV + MYIND2] = spatDerivs(xl1_2,xcc_2,xr1_2);
}




__global__ void zFlux_1warp(uint OFFSET,uint SPATCELL,Real* fz,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   if (threadIdx.z == 0) loadVelNbrs(MYBLOCK,nbrs,sha_nbr);

   // Fetch volume averages and z-derivatives for this block and -z neighbour, 
   // reconstruct face values and calculate z-flux:
   Real avg_neg,avg_pos,d1z_neg,d1z_pos,d2z_neg,d2z_pos;
   avg_neg = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZNEG]*SIZE_VELBLOCK+MYIND);
   d1z_neg = tex1Dfetch(texRef_d1z ,sha_nbr[NbrsVel::ZNEG]*SIZE_DERIV+MYIND);
   d2z_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1z_neg,d2z_neg);
   avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   d1z_pos = tex1Dfetch(texRef_d1z ,MYBLOCK*SIZE_DERIV+MYIND);
   d2z_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1z_pos,d2z_pos);
   fz[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxZ(avg_neg,avg_pos,sha_bparms,threadIdx.z  );

   MYIND += 2*WID2;
   avg_neg = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::ZNEG]*SIZE_VELBLOCK+MYIND);
   d1z_neg = tex1Dfetch(texRef_d1z ,sha_nbr[NbrsVel::ZNEG]*SIZE_DERIV+MYIND);
   d2z_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1z_neg,d2z_neg);
   avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   d1z_pos = tex1Dfetch(texRef_d1z ,MYBLOCK*SIZE_DERIV+MYIND);
   d2z_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1z_pos,d2z_pos);
   fz[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxZ(avg_neg,avg_pos,sha_bparms,threadIdx.z+2);
}

__global__ void zFlux_2warp(uint OFFSET,uint SPATCELL,Real* fz,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS*4];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*4];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the velocity block:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1z_pos = tex1Dfetch(texRef_d1z ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2z_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1z_pos,d2z_pos);
   
   cuint MYZNEG = sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::ZNEG];
   Real avg_neg = tex1Dfetch(texRef_avgs,MYZNEG*SIZE_VELBLOCK+MYIND);
   Real d1z_neg = tex1Dfetch(texRef_d1z ,MYZNEG*SIZE_DERIV+MYIND);
   Real d2z_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1z_neg,d2z_neg);
   
   fz[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxZ(avg_neg,avg_pos,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS,threadIdx.z);
}

template<uint WARPS> __global__ void zFlux_n2warp(uint OFFSET,uint SPATCELL,Real* fz,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS*WARPS*4];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*WARPS*4];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/4)*gridDim.x + blockIdx.x;
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%4);
   
   // Load parameters for the velocity block:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1z_pos = tex1Dfetch(texRef_d1z ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2z_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1z_pos,d2z_pos);
   
   cuint MYZNEG = sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::ZNEG];
   Real avg_neg = tex1Dfetch(texRef_avgs,MYZNEG*SIZE_VELBLOCK+MYIND);
   Real d1z_neg = tex1Dfetch(texRef_d1z ,MYZNEG*SIZE_DERIV+MYIND);
   Real d2z_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1z_neg,d2z_neg);
   
   fz[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxZ(avg_neg,avg_pos,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS,threadIdx.z%4);
}

__global__ void yFlux_1warp(uint OFFSET,uint SPATCELL,Real* fy,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);

   // Load parameters for the velocity block:
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   if (threadIdx.z == 0) loadVelNbrs(MYBLOCK,nbrs,sha_nbr);
   
   // Fetch volume averages and y-derivatives for this block and -y neighbour,
   // reconstruct face values and calculate y-flux:
   Real avg_neg = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YNEG]*SIZE_VELBLOCK+MYIND);
   Real d1y_neg = tex1Dfetch(texRef_d1y ,sha_nbr[NbrsVel::YNEG]*SIZE_DERIV+MYIND);
   Real d2y_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1y_neg,d2y_neg);
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1y_pos = tex1Dfetch(texRef_d1y ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2y_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1y_pos,d2y_pos);
   fy[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxY(avg_neg,avg_pos,sha_bparms);

   MYIND += 2*WID2;      
   avg_neg = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::YNEG]*SIZE_VELBLOCK+MYIND);
   d1y_neg = tex1Dfetch(texRef_d1y ,sha_nbr[NbrsVel::YNEG]*SIZE_DERIV+MYIND);
   d2y_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1y_neg,d2y_neg);
   avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   d1y_pos = tex1Dfetch(texRef_d1y ,MYBLOCK*SIZE_DERIV+MYIND);
   d2y_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1y_pos,d2y_pos);
   fy[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxY(avg_neg,avg_pos,sha_bparms);
}

__global__ void yFlux_2warp(uint OFFSET,uint SPATCELL,Real* fy,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS*4];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*4];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the velocity block:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1y_pos = tex1Dfetch(texRef_d1y ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2y_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1y_pos,d2y_pos);
   
   cuint MYYNEG = sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::YNEG];
   Real avg_neg = tex1Dfetch(texRef_avgs,MYYNEG*SIZE_VELBLOCK+MYIND);
   Real d1y_neg = tex1Dfetch(texRef_d1y ,MYYNEG*SIZE_DERIV+MYIND);
   Real d2y_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1y_neg,d2y_neg);
   
   fy[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxY(avg_neg,avg_pos,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}

template<uint WARPS> __global__ void yFlux_n2warp(uint OFFSET,uint SPATCELL,Real* fy,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS*WARPS*4];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*WARPS*4];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/4)*gridDim.x + blockIdx.x;
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%4);
   
   // Load parameters for the velocity block:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1y_pos = tex1Dfetch(texRef_d1y ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2y_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1y_pos,d2y_pos);
   
   cuint MYYNEG = sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::YNEG];
   Real avg_neg = tex1Dfetch(texRef_avgs,MYYNEG*SIZE_VELBLOCK+MYIND);
   Real d1y_neg = tex1Dfetch(texRef_d1y ,MYYNEG*SIZE_DERIV+MYIND);
   Real d2y_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1y_neg,d2y_neg);
    
   fy[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxY(avg_neg,avg_pos,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}

__global__ void xFlux_1warp(uint OFFSET,uint SPATCELL,Real* fx,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the velocity block:
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   if (threadIdx.z == 0) loadVelNbrs(MYBLOCK,nbrs,sha_nbr);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_neg = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XNEG]*SIZE_VELBLOCK+MYIND);
   Real d1x_neg = tex1Dfetch(texRef_d1x ,sha_nbr[NbrsVel::XNEG]*SIZE_DERIV+MYIND);
   Real d2x_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1x_neg,d2x_neg);
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1x_pos = tex1Dfetch(texRef_d1x ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2x_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1x_pos,d2x_pos);
   fx[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxX(avg_neg,avg_pos,sha_bparms);
   
   MYIND += 2*WID2;
   avg_neg = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::XNEG]*SIZE_VELBLOCK+MYIND);
   d1x_neg = tex1Dfetch(texRef_d1x ,sha_nbr[NbrsVel::XNEG]*SIZE_DERIV+MYIND);
   d2x_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1x_neg,d2x_neg);
   avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   d1x_pos = tex1Dfetch(texRef_d1x ,MYBLOCK*SIZE_DERIV+MYIND);
   d2x_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1x_pos,d2x_pos);
   fx[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxX(avg_neg,avg_pos,sha_bparms);
}

__global__ void xFlux_2warp(uint OFFSET,uint SPATCELL,Real* fx,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS*4];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*4];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the velocity block:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1x_pos = tex1Dfetch(texRef_d1x ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2x_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1x_pos,d2x_pos);
   
   cuint MYXNEG = sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::XNEG];
   Real avg_neg = tex1Dfetch(texRef_avgs,MYXNEG*SIZE_VELBLOCK+MYIND);
   Real d1x_neg = tex1Dfetch(texRef_d1x ,MYXNEG*SIZE_DERIV+MYIND);
   Real d2x_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1x_neg,d2x_neg);
   
   fx[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxX(avg_neg,avg_pos,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}

template<uint WARPS> __global__ void xFlux_n2warp(uint OFFSET,uint SPATCELL,Real* fx,Real* blockParams,uint* nbrs) {
   __shared__ Real sha_bparms[SIZE_BLOCKPARAMS*WARPS*4];
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*WARPS*4];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/4)*gridDim.x + blockIdx.x;
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%4);
   
   // Load parameters for the velocity block:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS + threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume averages and x-derivatives for this block and -x neighbour,
   // reconstruct face values and calculate x-flux:
   Real avg_pos = tex1Dfetch(texRef_avgs,MYBLOCK*SIZE_VELBLOCK+MYIND);
   Real d1x_pos = tex1Dfetch(texRef_d1x ,MYBLOCK*SIZE_DERIV+MYIND);
   Real d2x_pos = 0.0f;
   avg_pos = reconstruct_pos(avg_pos,d1x_pos,d2x_pos);
   
   cuint MYXNEG = sha_nbr[threadIdx.z*SIZE_NBRS_VEL + NbrsVel::XNEG];
   Real avg_neg = tex1Dfetch(texRef_avgs,MYXNEG*SIZE_VELBLOCK+MYIND);
   Real d1x_neg = tex1Dfetch(texRef_d1x ,MYXNEG*SIZE_DERIV+MYIND);
   Real d2x_neg = 0.0f;
   avg_neg = reconstruct_neg(avg_neg,d1x_neg,d2x_neg);
   
   fx[MYBLOCK*SIZE_FLUXS+MYIND] = spatialFluxX(avg_neg,avg_pos,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}





__global__ void propagateSpat_1warp(uint OFFSET,uint SPATCELL,Real* avgs,Real* fx,Real* fy,Real* fz,uint* nbrs,Real DT) {
   __shared__ Real sha_avg[SIZE_VELBLOCK];          // Shared mem array for volume averages
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYIND2 = MYIND1 + 2*WID2;
   
   if (threadIdx.z == 0) loadVelNbrs(MYBLOCK,nbrs,sha_nbr);
   
   Real F_neg1,F_neg2,F_pos1,F_pos2;

   // Fetch x-fluxes and calculate the contribution to total time derivative:
   F_neg1 = fx[MYBLOCK*SIZE_FLUXS + MYIND1];
   F_neg2 = fx[MYBLOCK*SIZE_FLUXS + MYIND2];
   F_pos1 = fx[sha_nbr[NbrsVel::XPOS]*SIZE_FLUXS + MYIND1];
   F_pos2 = fx[sha_nbr[NbrsVel::XPOS]*SIZE_FLUXS + MYIND2];
   sha_avg[MYIND1] = (F_neg1-F_pos1)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DX);
   sha_avg[MYIND2] = (F_neg2-F_pos2)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DX);

   // Fetch y-fluxes and calculate the contribution to total time derivative:
   F_neg1 = fy[MYBLOCK*SIZE_FLUXS + MYIND1];
   F_neg2 = fy[MYBLOCK*SIZE_FLUXS + MYIND2];
   F_pos1 = fy[sha_nbr[NbrsVel::YPOS]*SIZE_FLUXS + MYIND1];
   F_pos2 = fy[sha_nbr[NbrsVel::YPOS]*SIZE_FLUXS + MYIND2];
   sha_avg[MYIND1] += (F_neg1-F_pos1)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DY);
   sha_avg[MYIND2] += (F_neg2-F_pos2)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DY);

   // Fetch z-fluxes and calculate the contribution to total time derivative:
   F_neg1 = fz[MYBLOCK*SIZE_FLUXS + MYIND1];
   F_neg2 = fz[MYBLOCK*SIZE_FLUXS + MYIND2];
   F_pos1 = fz[sha_nbr[NbrsVel::ZPOS]*SIZE_FLUXS + MYIND1];
   F_pos2 = fz[sha_nbr[NbrsVel::ZPOS]*SIZE_FLUXS + MYIND2];
   sha_avg[MYIND1] += (F_neg1-F_pos1)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DZ);
   sha_avg[MYIND2] += (F_neg2-F_pos2)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DZ);

   // Store new volume averages:
   avgs[MYBLOCK*SIZE_VELBLOCK+MYIND1] += sha_avg[MYIND1];
   avgs[MYBLOCK*SIZE_VELBLOCK+MYIND2] += sha_avg[MYIND2];
}

__global__ void propagateSpat_2warp(uint OFFSET,uint SPATCELL,Real* avgs,Real* fx,Real* fy,Real* fz,uint* nbrs,Real DT) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*4];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   // Fetch neighbour indices:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   
   Real F_neg,F_pos;
   // Fetch x-fluxes and calculate the contribution to total time derivative:
   F_neg = fx[MYBLOCK*SIZE_FLUXS + MYIND];
   uint MYNBR = sha_nbr[threadIdx.z*SIZE_NBRS_VEL+NbrsVel::XPOS];
   F_pos = fx[MYNBR*SIZE_FLUXS + MYIND];
   Real myavg = (F_neg-F_pos)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DX);
   
   // Fetch y-fluxes and calculate the contribution to total time derivative:
   F_neg = fy[MYBLOCK*SIZE_FLUXS + MYIND];
   MYNBR = sha_nbr[threadIdx.z*SIZE_NBRS_VEL+NbrsVel::YPOS];
   F_pos = fy[MYNBR*SIZE_FLUXS + MYIND];
   myavg += (F_neg-F_pos)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DY);
   
   // Fetch z-fluxes and calculate the contribution to total time derivative:
   F_neg = fz[MYBLOCK*SIZE_FLUXS + MYIND];
   MYNBR = sha_nbr[threadIdx.z*SIZE_NBRS_VEL+NbrsVel::ZPOS];
   F_pos = fz[MYNBR*SIZE_FLUXS + MYIND];
   myavg += (F_neg-F_pos)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DZ);
   
   // Store new volume averages:
   avgs[MYBLOCK*SIZE_VELBLOCK+MYIND] += myavg;
}

template<uint WARPS> __global__ void propagateSpat_n2warp(uint OFFSET,uint SPATCELL,Real* avgs,Real* fx,Real* fy,Real* fz,uint* nbrs,Real DT) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL*4*WARPS];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/4)*gridDim.x + blockIdx.x;
   cuint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%4);
   // Fetch neighbour indices:
   sha_nbr[threadIdx.z*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_nbrs,MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID+threadIdx.x);
   
   Real F_neg,F_pos;
   // Fetch x-fluxes and calculate the contribution to total time derivative:
   F_neg = fx[MYBLOCK*SIZE_FLUXS + MYIND];
   uint MYNBR = sha_nbr[threadIdx.z*SIZE_NBRS_VEL+NbrsVel::XPOS];
   F_pos = fx[MYNBR*SIZE_FLUXS + MYIND];
   Real myavg = (F_neg-F_pos)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DX);
   
   // Fetch y-fluxes and calculate the contribution to total time derivative:
   F_neg = fy[MYBLOCK*SIZE_FLUXS + MYIND];
   MYNBR = sha_nbr[threadIdx.z*SIZE_NBRS_VEL+NbrsVel::YPOS];
   F_pos = fy[MYNBR*SIZE_FLUXS + MYIND];
   myavg += (F_neg-F_pos)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DY);
   
   // Fetch z-fluxes and calculate the contribution to total time derivative:
   F_neg = fz[MYBLOCK*SIZE_FLUXS + MYIND];
   MYNBR = sha_nbr[threadIdx.z*SIZE_NBRS_VEL+NbrsVel::ZPOS];
   F_pos = fz[MYNBR*SIZE_FLUXS + MYIND];
   myavg += (F_neg-F_pos)*DT/tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_BLOCKPARAMS+CellParams::DZ);
       
   // Store new volume averages:
   avgs[MYBLOCK*SIZE_VELBLOCK+MYIND] += myavg;
}




inline void calcSpatialDerivs(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   // Calculate the derivatives in spatial coordinates:
   offset = grid[cell]->velBlockIndex;
   #ifdef SPAT_DERIVS_1WARP // 1 warp calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;                             
   calcSpatDerivs_1warp<<<gridSize,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						 deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef SPAT_DERIVS_NWARP // n warps calculating a single velocity grid block each
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x*WARPS_SPAT_DERIVS);
   gridSize.z = 1;
   
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2*WARPS_SPAT_DERIVS;
   calcSpatDerivs_nwarp<WARPS_SPAT_DERIVS><<<gridSize,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
								    deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   #endif 
   #ifdef SPAT_DERIVS_2WARP // 2 warps calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4;
   calcSpatDerivs_2warp<<<gridSize,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						 deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif 
   #ifdef SPAT_DERIVS_N2WARP // n 2-warp thread blocks calculating a single velocity grid block each
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x*WARPS_SPAT_DERIVS);
   gridSize.z = 1;
   
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4*WARPS_SPAT_DERIVS;
   calcSpatDerivs_n2warp<WARPS_SPAT_DERIVS><<<gridSize,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
								     deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   #endif
   
   // Calculate the leftover blocks:
   #ifdef EXTRA_SPAT_DERIVS_1WARP
   blockSizes.z = 2;
   calcSpatDerivs_1warp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						      deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),deviceGrid.getNbrsVel());
   #endif
   #ifdef EXTRA_SPAT_DERIVS_2WARP
   blockSizes.z = 4;
   calcSpatDerivs_2warp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						      deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),deviceGrid.getNbrsVel());
   #endif
}

inline void calcSpatialFluxesX(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex; 
   #ifdef FLUXES_1WARP // 1 warp calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   xFlux_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFx(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef FLUXES_2WARP // 2 warps calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4;
   xFlux_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFx(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef FLUXES_N2WARP // n 2-warp thread blocks calculating a single velocity grid block each
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x*WARPS_SPAT_DERIVS);
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4*FLUXES_WARPS;
   xFlux_n2warp<FLUXES_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFx(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y*FLUXES_WARPS;
   #endif
   
   // Calculate the leftover blocks:
   #ifdef EXTRA_FLUXES_1WARP
   blockSizes.z = 2;
   xFlux_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFx(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   #endif
   #ifdef EXTRA_FLUXES_2WARP
   blockSizes.z = 4;
   xFlux_2warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFx(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   #endif
}

inline void calcSpatialFluxesY(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef FLUXES_1WARP // 1 warp calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   yFlux_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFy(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef FLUXES_2WARP // 2 warps calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4;
   yFlux_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFy(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef FLUXES_N2WARP // n 2-warp thread blocks calculating a single velocity grid block each
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x*WARPS_SPAT_DERIVS);
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4*FLUXES_WARPS;
   yFlux_n2warp<FLUXES_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFy(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y*FLUXES_WARPS;
   #endif
   
   // Calculate the leftover blocks:
   #ifdef EXTRA_FLUXES_1WARP
   blockSizes.z = 2;
   yFlux_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFy(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   #endif
   #ifdef EXTRA_FLUXES_2WARP
   blockSizes.z = 4;
   yFlux_2warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFy(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   #endif
}

inline void calcSpatialFluxesZ(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef FLUXES_1WARP // 1 warp calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   zFlux_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFz(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef FLUXES_2WARP // 2 warps calculating a single velocity grid block
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4;
   zFlux_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFz(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef FLUXES_N2WARP // n 2-warp thread blocks calculating a single velocity grid block each
   gridSize.x = SPAT_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x*WARPS_SPAT_DERIVS);
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*WARPS_SPAT_DERIVS;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4*FLUXES_WARPS;
   zFlux_n2warp<FLUXES_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFz(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y*FLUXES_WARPS;
   #endif

   // Calculate the leftover blocks:
   #ifdef EXTRA_FLUXES_1WARP
   blockSizes.z = 2;
   zFlux_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFz(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   #endif
   #ifdef EXTRA_FLUXES_2WARP
   blockSizes.z = 4;
   zFlux_2warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getFz(),deviceGrid.getBlockParams(),deviceGrid.getNbrsVel());
   #endif
}

inline void spatialPropagate(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes,creal& DT) {
   offset = grid[cell]->velBlockIndex;
   #ifdef PROPAG_1WARP // 1 warp calculates a single velocity grid block
   gridSize.x = PROPAG_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   propagateSpat_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
						deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
						deviceGrid.getNbrsVel(),DT);
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef PROPAG_2WARP // 2 warps calculate a single velocity grid block
   gridSize.x = PROPAG_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4;
   propagateSpat_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
						deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
						deviceGrid.getNbrsVel(),DT);
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef PROPAG_N2WARP // N 2-warp thread blocks calculate a single velocity grid block each
   gridSize.x = PROPAG_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x*PROPAG_WARPS);
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*PROPAG_WARPS;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 4*PROPAG_WARPS;
   propagateSpat_n2warp<PROPAG_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
							       deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
							       deviceGrid.getNbrsVel(),DT);
   offset += gridSize.x*gridSize.y*PROPAG_WARPS;
   #endif
   
   // Calculate the leftover blocks:
   #ifdef EXTRA_PROPAG_1WARP
   blockSizes.z = 2;
   propagateSpat_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
						     deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
						     deviceGrid.getNbrsVel(),DT);
   #endif
   #ifdef EXTRA_PROPAG_2WARP
   blockSizes.z = 4;
   propagateSpat_2warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
						     deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
						     deviceGrid.getNbrsVel(),DT);
   #endif   
}

bool translation_step1(Grid& grid,DeviceGrid& deviceGrid,creal& DT) {
   bool success = true;

   // Bind textures:
   size_t texOffset;
   if (bindTexture(texRef_avgs,deviceGrid.getBlockArray(),deviceGrid.sizeBlockArray(),texOffset) == false) success = false;
   if (bindTexture(texRef_cellParams,deviceGrid.getCellParams(),deviceGrid.sizeCellParams(),texOffset) == false) success = false;
   if (bindTexture(texRef_nbrs,deviceGrid.getNbrsVel(),deviceGrid.sizeNbrsVel(),texOffset) == false) success = false;
   if (success == false) return success;
   
   // Go though each spatial cell:
   uint offset;
   dim3 gridSize;
   dim3 gridSizeExtra;
   dim3 blockSizes(WID,WID,1);
   for (uint cell=0; cell<grid.size(); ++cell) {
      calcSpatialDerivs(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
   }
   return success;
}

bool translation_step2(Grid& grid,DeviceGrid& deviceGrid,creal& DT) {
   bool success = true;

   // Bind textures:
   size_t texOffset;
   if (bindTexture(texRef_avgs,deviceGrid.getBlockArray(),deviceGrid.sizeBlockArray(),texOffset) == false) success = false;
   if (bindTexture(texRef_nbrs,deviceGrid.getNbrsVel(),deviceGrid.sizeNbrsVel(),texOffset) == false) success = false;
   if (bindTexture(texRef_bparms,deviceGrid.getBlockParams(),deviceGrid.sizeBlockParams(),texOffset) == false) success = false;
   if (bindTexture(texRef_d1x,deviceGrid.getD1x(),deviceGrid.sizeD1x(),texOffset) == false) success = false;
   if (bindTexture(texRef_d1y,deviceGrid.getD1y(),deviceGrid.sizeD1y(),texOffset) == false) success = false;
   if (bindTexture(texRef_d1z,deviceGrid.getD1z(),deviceGrid.sizeD1z(),texOffset) == false) success = false;
   if (success == false) return false;
   
   // Go though each spatial cell:
   uint offset;
   dim3 gridSize;
   dim3 gridSizeExtra;
   dim3 blockSizes(WID,WID,2);
   for (uint cell=0; cell<grid.size(); ++cell) {
      // Calculate x,y,z spatial fluxes:
      calcSpatialFluxesZ(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
      calcSpatialFluxesY(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
      calcSpatialFluxesX(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
   }   
   return success;
}

bool translation_step3(Grid& grid,DeviceGrid& deviceGrid,creal& DT) {
   bool success = true;

   // Bind textures:
   size_t texOffset;
   if (bindTexture(texRef_cellParams,deviceGrid.getCellParams(),deviceGrid.sizeCellParams(),texOffset) == false) success = false;
   if (bindTexture(texRef_nbrs,deviceGrid.getNbrsVel(),deviceGrid.sizeNbrsVel(),texOffset) == false) success = false;
   if (success == false) return success;
   
   // Propagate volume averages in spatial space:
   uint offset;
   dim3 gridSize;
   dim3 gridSizeExtra;
   dim3 blockSizes(WID,WID,2);
   for (uint cell=0; cell<grid.size(); ++cell) {
      spatialPropagate(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes,DT);
   }
   return success;
}

