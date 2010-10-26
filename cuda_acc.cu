#include "cuda_common.cu"
#include "cudalaunch.h"
#include "project.cu"

__host__ bool bindTextureAcc(texture<uint,1,cudaReadModeElementType>& texRef,uint* arrptr,cuint& BYTES,size_t& offset) {
   bool success = true;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindUnsigned);
   cudaError_t cudaError = cudaBindTexture(&offset,&texRef,arrptr,&channelDesc,BYTES);
   if (cudaError != cudaSuccess) {
      logger << "(CUDA_TRANS): Failed to bind textures uint*!" << endl;
      success = false;
   }
   return success;
}

__host__ bool bindTextureAcc(texture<real,1,cudaReadModeElementType>& texRef,real* arrptr,cuint& BYTES,size_t& offset) {
   bool success = true;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
   cudaError_t cudaError = cudaBindTexture(&offset,&texRef,arrptr,&channelDesc,BYTES);
   if (cudaError != cudaSuccess) {
      logger << "(CUDA_TRANS): Failed to bind textures to real*!" << endl;
      success = false;
   }
   return success;
}

// ------------------- DEVICE FUNCTIONS ------------------

__device__ void velDerivs(uint MYIND,real* sha_avg,real* d1,real* d2) {
   // Read a five-point stencil into registers:
   real xl2 = sha_avg[MYIND + 2*WID2 - 2*WID2];
   real xl1 = sha_avg[MYIND + 2*WID2 -   WID2];
   real xcc = sha_avg[MYIND + 2*WID2         ];
   real xr1 = sha_avg[MYIND + 2*WID2 +   WID2];
   real xr2 = sha_avg[MYIND + 2*WID2 + 2*WID2];
   
   // Calculate 1st and 2nd derivatives:
   d1[MYIND] = superbee(xl1,xcc,xr1);
   d2[MYIND] = 0.0f;
}

// ---------------------- KERNELS ------------------------

__global__ void copyVelBlocks_1warp(cuint OFFSET,real* avgs,real* avgnbrx,real* avgnbry,real* avgnbrz,uint* nbrs) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   cuint MYIND2 = MYIND1 + 2*WID2;
   uint myread;
   
   // Creates a divergent branch but it is faster than reading neighbours twice
   if (threadIdx.z==0) loadVelNbrs(MYBLOCK,nbrs,sha_nbr);
   
   // Either calculate ghost cell data using a boundary condition function or copy averages from +/- vz-neighbours:
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VZ_NEG_BND) == NbrsVel::VZ_NEG_BND) {
      avgnbrz[MYBLOCK*SIZE_BOUND+MYIND1] = 0.0f;
   } else {
      avgnbrz[MYBLOCK*SIZE_BOUND+MYIND1] = avgs[sha_nbr[NbrsVel::VZNEG]*SIZE_VELBLOCK+MYIND2];
   }
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VZ_POS_BND) == NbrsVel::VZ_POS_BND) {
      avgnbrz[MYBLOCK*SIZE_BOUND+MYIND2] = 0.0f;
   } else {
      avgnbrz[MYBLOCK*SIZE_BOUND+MYIND2] = avgs[sha_nbr[NbrsVel::VZPOS]*SIZE_VELBLOCK+MYIND1];
   }

   // Either calculate ghost cell data using a boundary condition function or copy averages from +/- vy-neighbours:
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VY_NEG_BND) == NbrsVel::VY_NEG_BND) {
      avgnbry[MYBLOCK*SIZE_BOUND+MYIND1] = 0.0f;
   } else {
      myread = tindex3(threadIdx.x,threadIdx.z+2,threadIdx.y);
      avgnbry[MYBLOCK*SIZE_BOUND+MYIND1] = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::VYNEG]*SIZE_VELBLOCK+myread);
   }   
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VY_POS_BND) == NbrsVel::VY_POS_BND) {
      avgnbry[MYBLOCK*SIZE_BOUND+MYIND2] = 0.0f;
   } else {
      myread = tindex3(threadIdx.x,threadIdx.z  ,threadIdx.y);
      avgnbry[MYBLOCK*SIZE_BOUND+MYIND2] = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::VYPOS]*SIZE_VELBLOCK+myread);
   }
   
   // Either calculate ghost cell data using a boundary condition function or copy averages from +/- vx-neighbours:
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VX_NEG_BND) == NbrsVel::VX_NEG_BND) {
      avgnbrx[MYBLOCK*SIZE_BOUND+MYIND1] = 0.0f;
   } else {
      myread = tindex3(threadIdx.z+2,threadIdx.x,threadIdx.y);
      avgnbrx[MYBLOCK*SIZE_BOUND+MYIND1] = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::VXNEG]*SIZE_VELBLOCK+myread);
   }
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VX_POS_BND) == NbrsVel::VX_POS_BND) {
      avgnbrx[MYBLOCK*SIZE_BOUND+MYIND2] = 0.0f;
   } else {
      myread = tindex3(threadIdx.z  ,threadIdx.x,threadIdx.y);
      avgnbrx[MYBLOCK*SIZE_BOUND+MYIND2] = tex1Dfetch(texRef_avgs,sha_nbr[NbrsVel::VXPOS]*SIZE_VELBLOCK+myread);
   }
}


__global__ void calcVelDerivs_1warp(uint OFFSET,real* avgs,real* avgnbrx,real* avgnbry,real* avgnbrz,
				    real* d1x,real* d1y,real* d1z,real* d2x,real* d2y,real* d2z) {
   __shared__ real sha_avg[SIZE_VELBLOCK+SIZE_BOUND]; // Shared mem array for holding volume averages + ghosts
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Fetch volume averages of the block to the middle of the shared mem array
   sha_avg[MYIND+2*WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND       ];
   sha_avg[MYIND+4*WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND+2*WID2];

   // Fetch z-ghosts and calculate z-derivatives:
   sha_avg[MYIND       ] = avgnbrz[MYBLOCK*SIZE_BOUND+MYIND       ];
   sha_avg[MYIND+6*WID2] = avgnbrz[MYBLOCK*SIZE_BOUND+MYIND+2*WID2];
   velDerivs(MYIND       ,sha_avg,d1z+MYBLOCK*SIZE_DERIV,d2z+MYBLOCK*SIZE_DERIV);
   velDerivs(MYIND+2*WID2,sha_avg,d1z+MYBLOCK*SIZE_DERIV,d2z+MYBLOCK*SIZE_DERIV);

   // Fetch y-ghosts, transpose, and calculate y-derivatives:
   sha_avg[MYIND       ] = avgnbry[MYBLOCK*SIZE_BOUND+MYIND       ];
   sha_avg[MYIND+6*WID2] = avgnbry[MYBLOCK*SIZE_BOUND+MYIND+2*WID2];
   transpose_yz_1warp(sha_avg+2*WID2);
   velDerivs(MYIND       ,sha_avg,d1y+MYBLOCK*SIZE_DERIV,d2y+MYBLOCK*SIZE_DERIV);
   velDerivs(MYIND+2*WID2,sha_avg,d1y+MYBLOCK*SIZE_DERIV,d2y+MYBLOCK*SIZE_DERIV);

   // Fetch x-ghosts, transpose, and calculate x-derivatives:
   sha_avg[MYIND       ] = avgnbrx[MYBLOCK*SIZE_BOUND+MYIND       ];
   sha_avg[MYIND+6*WID2] = avgnbrx[MYBLOCK*SIZE_BOUND+MYIND+2*WID2];
   transpose_xz_1warp(sha_avg+2*WID2);
   velDerivs(MYIND       ,sha_avg,d1x+MYBLOCK*SIZE_DERIV,d2x+MYBLOCK*SIZE_DERIV);
   velDerivs(MYIND+2*WID2,sha_avg,d1x+MYBLOCK*SIZE_DERIV,d2x+MYBLOCK*SIZE_DERIV);
}

__global__ void copyVelDerivs_hwarp(uint OFFSET,real* dx,real* dy,real* dz,real* dxnbr,real* dynbr,real* dznbr,uint* nbrs) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND = tindex2(threadIdx.x,threadIdx.y);
   
   loadVelNbrs(MYBLOCK,nbrs,sha_nbr);
   
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VX_NEG_BND) == NbrsVel::VX_NEG_BND) {
      dxnbr[MYBLOCK*SIZE_BDERI + MYIND] = 0.0f;
   } else {
      dxnbr[MYBLOCK*SIZE_BDERI + MYIND] = dx[sha_nbr[NbrsVel::VXNEG]*SIZE_DERIV + MYIND + 3*WID2];
   }
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VY_NEG_BND) == NbrsVel::VY_NEG_BND) {
      dynbr[MYBLOCK*SIZE_BDERI + MYIND] = 0.0f;
   } else {
      dynbr[MYBLOCK*SIZE_BDERI + MYIND] = dy[sha_nbr[NbrsVel::VYNEG]*SIZE_DERIV + MYIND + 3*WID2];
   }
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VZ_NEG_BND) == NbrsVel::VZ_NEG_BND) {
      dznbr[MYBLOCK*SIZE_BDERI + MYIND] = 0.0f;
   } else {
      dznbr[MYBLOCK*SIZE_BDERI + MYIND] = dz[sha_nbr[NbrsVel::VZNEG]*SIZE_DERIV + MYIND + 3*WID2];
   }
}

// ----------------------------------------------------------------------------------
// ------------ FUNCTIONS FOR CALCULATING Z-COMPONENT OF VELOCITY FLUXES ------------
// ----------------------------------------------------------------------------------

__global__ void vzFlux_1warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbrz,real* d1z,real* d2z,real* d1znbr,real* d2znbr,real* fz,real* blockParams) {
   __shared__ real sha_avg[SIZE_VELBLOCK+WID2];
   __shared__ real sha_d1z[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_d2z[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_cparms[SIZE_CELLPARAMS];
   __shared__ real sha_bparms[SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_cparms[MYIND] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+MYIND);
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   
   // Fetch volume averages and first layer of ghost cells of -z boundary:
   sha_avg[MYIND+  WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND       ];
   sha_avg[MYIND+3*WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_avg[MYIND] = avgnbrz[MYBLOCK*SIZE_BOUND+MYIND+WID2];

   // Fetch 1st derivatives and ghost cells:
   sha_d1z[MYIND+  WID2] = d1z[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d1z[MYIND+3*WID2] = d1z[MYBLOCK*SIZE_DERIV+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_d1z[MYIND] = d1znbr[MYBLOCK*SIZE_BDERI+MYIND];

   // Fetch 2nd derivatives and ghost cells:
   sha_d2z[MYIND+  WID2] = d2z[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d2z[MYIND+3*WID2] = d2z[MYBLOCK*SIZE_DERIV+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_d2z[MYIND] = d2znbr[MYBLOCK*SIZE_BDERI+MYIND];

   // Reconstruct negative and positive side values at the -vz face, and store calculated flux:
   for (int cntr=0; cntr<2; ++cntr) {
      real avg_neg = reconstruct_neg(sha_avg[MYIND     ],sha_d1z[MYIND     ],sha_d2z[MYIND     ]);
      real avg_pos = reconstruct_pos(sha_avg[MYIND+WID2],sha_d1z[MYIND+WID2],sha_d2z[MYIND+WID2]);
      fz[MYBLOCK*SIZE_FLUXS+MYIND] = velocityFluxZ(avg_neg,avg_pos,sha_cparms,sha_bparms);
      MYIND += 2*WID2;
   }
}

template<uint WARPS> __global__ void vzFlux_n1warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbrz,real* d1z,real* d2z,real* d1znbr,real* d2znbr,real* fz,real* blockParams) {
   cuint SIZE_VEL = SIZE_VELBLOCK+WID2;
   cuint SIZE_DER = SIZE_DERIV+SIZE_BDERI;
   
   __shared__ real sha_avg[WARPS*SIZE_VEL];
   __shared__ real sha_d1z[WARPS*SIZE_DER];
   __shared__ real sha_d2z[WARPS*SIZE_DER];
   __shared__ real sha_cparms[WARPS*2*SIZE_CELLPARAMS];
   __shared__ real sha_bparms[WARPS*2*SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/2)*gridDim.x + blockIdx.x;
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%2);
   cuint MYIND2 = MYIND1 + 2*WID2;
   
   // Load parameters for the spatial cell & velocity block:
   sha_cparms[threadIdx.z*SIZE_CELLPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS+ threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS+threadIdx.y*WID+threadIdx.x);

   // Fetch volume averages and first layer of ghost cells of -z boundary:
   sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND1];
   sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND2];
   if (threadIdx.z % 2 == 0) sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1] = avgnbrz[MYBLOCK*SIZE_BOUND+MYIND1+WID2];
   //__syncthreads();
   
   // Fetch 1st derivatives and ghost cells:
   sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2] = d1z[MYBLOCK*SIZE_DERIV+MYIND1];
   sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2] = d1z[MYBLOCK*SIZE_DERIV+MYIND2];
   if (threadIdx.z % 2 == 0) sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND1] = d1znbr[MYBLOCK*SIZE_BDERI+MYIND1];
   //__syncthreads();
   
   // Fetch 2nd derivatives and ghost cells:
   sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2] = d2z[MYBLOCK*SIZE_DERIV+MYIND1];
   sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2] = d2z[MYBLOCK*SIZE_DERIV+MYIND2];
   if (threadIdx.z % 2 == 0) sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND1] = d2znbr[MYBLOCK*SIZE_BDERI+MYIND1];
   //__syncthreads();

   real avg_neg,avg_pos;
   avg_neg = reconstruct_neg(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1     ],sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND1     ],sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND1     ]);
   avg_pos = reconstruct_pos(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1+WID2],sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2],sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2]);
   fz[MYBLOCK*SIZE_FLUXS+MYIND1] = velocityFluxZ(avg_neg,avg_pos,sha_cparms+threadIdx.z*SIZE_CELLPARAMS,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
   
   avg_neg = reconstruct_neg(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2     ],sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND2     ],sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND2     ]);
   avg_pos = reconstruct_pos(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2+WID2],sha_d1z[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2],sha_d2z[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2]);
   fz[MYBLOCK*SIZE_FLUXS+MYIND2] = velocityFluxZ(avg_neg,avg_pos,sha_cparms+threadIdx.z*SIZE_CELLPARAMS,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}

__global__ void vzFlux_2warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbrz,real* d1z,real* d2z,real* d1znbr,real* d2znbr,real* fz,real* blockParams) {
   __shared__ real sha_avg[SIZE_VELBLOCK+WID2];
   __shared__ real sha_d1z[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_d2z[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_cparms[SIZE_CELLPARAMS];
   __shared__ real sha_bparms[SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   
   // Fetch volume averages and first+second derivatives for the block:
   sha_avg[MYIND+  WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND];
   sha_d1z[MYIND+  WID2] = d1z[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d2z[MYIND+  WID2] = d2z[MYBLOCK*SIZE_DERIV+MYIND       ];

   // Fetch ghost cells & block parameters:
   if (threadIdx.z == 1) sha_d1z[threadIdx.y*WID+threadIdx.x] = d1znbr[MYBLOCK*SIZE_BDERI+threadIdx.y*WID+threadIdx.x];
   if (threadIdx.z == 2) sha_d2z[threadIdx.y*WID+threadIdx.x] = d2znbr[MYBLOCK*SIZE_BDERI+threadIdx.y*WID+threadIdx.x];
   if (threadIdx.z == 3) sha_cparms[threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+threadIdx.y*WID+threadIdx.x);
   if (threadIdx.z == 0) sha_avg[threadIdx.y*WID+threadIdx.x] = avgnbrz[MYBLOCK*SIZE_BOUND+WID2+threadIdx.y*WID+threadIdx.x];
   __syncthreads();
   
   real avg_neg = reconstruct_neg(sha_avg[MYIND     ],sha_d1z[MYIND     ],sha_d2z[MYIND     ]);
   real avg_pos = reconstruct_pos(sha_avg[MYIND+WID2],sha_d1z[MYIND+WID2],sha_d2z[MYIND+WID2]);
   fz[MYBLOCK*SIZE_FLUXS+MYIND] = velocityFluxZ(avg_neg,avg_pos,sha_cparms,sha_bparms);
}

// ----------------------------------------------------------------------------------
// ------------ FUNCTIONS FOR CALCULATING Y-COMPONENT OF VELOCITY FLUXES ------------
// ----------------------------------------------------------------------------------

__global__ void vyFlux_1warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbry,real* d1y,real* d2y,real* d1ynbr,real* d2ynbr,real* fy,real* blockParams) {
   __shared__ real sha_avg[SIZE_VELBLOCK+WID2];
   __shared__ real sha_d1y[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_d2y[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_cparms[SIZE_CELLPARAMS];
   __shared__ real sha_bparms[SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_cparms[MYIND] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+MYIND);
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   
   // Fetch volume avgs to the middle of the shared mem array
   sha_avg[MYIND+  WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND       ];
   sha_avg[MYIND+3*WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_avg[MYIND] = avgnbry[MYBLOCK*SIZE_BOUND+MYIND+WID2];
   // Fetch 1st derivatives
   sha_d1y[MYIND+  WID2] = d1y[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d1y[MYIND+3*WID2] = d1y[MYBLOCK*SIZE_DERIV+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_d1y[MYIND] = d1ynbr[MYBLOCK*SIZE_BDERI+MYIND];
   // Fetch 2nd derivatives
   sha_d2y[MYIND+  WID2] = d2y[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d2y[MYIND+3*WID2] = d2y[MYBLOCK*SIZE_DERIV+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_d2y[MYIND] = d2ynbr[MYBLOCK*SIZE_BDERI+MYIND];
   
   // Transpose volume averages in yz-direction
   transpose_yz_1warp(sha_avg+WID2);
   
   // Reconstruct negative and positive side values at the -vy face, and store calculated flux.
   for (int cntr=0; cntr<2; ++cntr) {
      real avg_neg = reconstruct_neg(sha_avg[MYIND     ],sha_d1y[MYIND     ],sha_d2y[MYIND     ]);
      real avg_pos = reconstruct_pos(sha_avg[MYIND+WID2],sha_d1y[MYIND+WID2],sha_d2y[MYIND+WID2]);
      fy[MYBLOCK*SIZE_FLUXS+MYIND] = velocityFluxY(avg_neg,avg_pos,sha_cparms,sha_bparms);
      MYIND += 2*WID2;
   }
}

template<uint WARPS> __global__ void vyFlux_n1warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbry,real* d1y,real* d2y,real* d1ynbr,real* d2ynbr,real* fy,real* blockParams) {
   cuint SIZE_VEL = SIZE_VELBLOCK+WID2;
   cuint SIZE_DER = SIZE_DERIV+SIZE_BDERI;
   
   __shared__ real sha_avg[WARPS*SIZE_VEL];
   __shared__ real sha_d1y[WARPS*SIZE_DER];
   __shared__ real sha_d2y[WARPS*SIZE_DER];
   __shared__ real sha_cparms[WARPS*2*SIZE_CELLPARAMS];
   __shared__ real sha_bparms[WARPS*2*SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/2)*gridDim.x + blockIdx.x;
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%2);
   cuint MYIND2 = MYIND1 + 2*WID2;
   
   // Load parameters for the spatial cell & velocity block:
   sha_cparms[threadIdx.z*SIZE_CELLPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS+ threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS+threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume avgs to the middle of the shared mem array
   sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND1];
   sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND2];
   if (threadIdx.z % 2 == 0) sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1] = avgnbry[MYBLOCK*SIZE_BOUND+MYIND1+WID2];
   
   // Fetch 1st derivatives and ghost cells:
   sha_d1y[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2] = d1y[MYBLOCK*SIZE_DERIV+MYIND1];
   sha_d1y[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2] = d1y[MYBLOCK*SIZE_DERIV+MYIND2];
   if (threadIdx.z % 2 == 0) sha_d1y[(threadIdx.z/2)*SIZE_VEL+MYIND1] = d1ynbr[MYBLOCK*SIZE_BDERI+MYIND1];
   
   // Fetch 2nd derivatives and ghost cells:
   sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2] = d2y[MYBLOCK*SIZE_DERIV+MYIND1];
   sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2] = d2y[MYBLOCK*SIZE_DERIV+MYIND2];
   if (threadIdx.z % 2 == 0) sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND1] = d2ynbr[MYBLOCK*SIZE_BDERI+MYIND1];
   
   transpose_yz_1warp(sha_avg+(threadIdx.z/2)*SIZE_VEL+WID2);
   
   // Reconstruct negative and positive side values at the -vy face, and store calculated flux.
   real avg_neg,avg_pos;
   avg_neg = reconstruct_neg(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1     ],sha_d1y[(threadIdx.z/2)*SIZE_DER+MYIND1     ],sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND1     ]);
   avg_pos = reconstruct_pos(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1+WID2],sha_d1y[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2],sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2]);
   fy[MYBLOCK*SIZE_FLUXS+MYIND1] = velocityFluxY(avg_neg,avg_pos,sha_cparms+threadIdx.z*SIZE_CELLPARAMS,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
   
   avg_neg = reconstruct_neg(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2     ],sha_d1y[(threadIdx.z/2)*SIZE_DER+MYIND2     ],sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND2     ]);
   avg_pos = reconstruct_pos(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2+WID2],sha_d1y[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2],sha_d2y[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2]);
   fy[MYBLOCK*SIZE_FLUXS+MYIND2] = velocityFluxY(avg_neg,avg_pos,sha_cparms+threadIdx.z*SIZE_CELLPARAMS,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}

__global__ void vyFlux_2warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbry,real* d1y,real* d2y,real* d1ynbr,real* d2ynbr,real* fy,real* blockParams) {
   __shared__ real sha_avg[SIZE_VELBLOCK+WID2];
   __shared__ real sha_d1y[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_d2y[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_cparms[SIZE_CELLPARAMS];
   __shared__ real sha_bparms[SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   
   // Fetch volume averages and first+second derivatives for the block:
   sha_avg[MYIND+  WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND];
   sha_d1y[MYIND+  WID2] = d1y[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d2y[MYIND+  WID2] = d2y[MYBLOCK*SIZE_DERIV+MYIND       ];
   
   // Fetch ghost cells & block parameters:
   if (threadIdx.z == 1) sha_d1y[threadIdx.y*WID+threadIdx.x] = d1ynbr[MYBLOCK*SIZE_BDERI+threadIdx.y*WID+threadIdx.x];
   if (threadIdx.z == 2) sha_d2y[threadIdx.y*WID+threadIdx.x] = d2ynbr[MYBLOCK*SIZE_BDERI+threadIdx.y*WID+threadIdx.x];
   if (threadIdx.z == 3) sha_cparms[threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+threadIdx.y*WID+threadIdx.x);
   if (threadIdx.z == 0) sha_avg[threadIdx.y*WID+threadIdx.x] = avgnbry[MYBLOCK*SIZE_BOUND+WID2+threadIdx.y*WID+threadIdx.x];
   __syncthreads();
   
   transpose_yz_2warp(sha_avg+WID2);
   
   // Reconstruct negative and positive side values at the -vy face, and store calculated flux.
   real avg_neg = reconstruct_neg(sha_avg[MYIND     ],sha_d1y[MYIND     ],sha_d2y[MYIND     ]);
   real avg_pos = reconstruct_pos(sha_avg[MYIND+WID2],sha_d1y[MYIND+WID2],sha_d2y[MYIND+WID2]);
   fy[MYBLOCK*SIZE_FLUXS+MYIND] = velocityFluxY(avg_neg,avg_pos,sha_cparms,sha_bparms);
}

// ----------------------------------------------------------------------------------
// ------------ FUNCTIONS FOR CALCULATING X-COMPONENT OF VELOCITY FLUXES ------------
// ----------------------------------------------------------------------------------

__global__ void vxFlux_1warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbrx,real* d1x,real* d2x,real* d1xnbr,real* d2xnbr,real* fx,real* blockParams) {
   __shared__ real sha_avg[SIZE_VELBLOCK+WID2];
   __shared__ real sha_d1x[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_d2x[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_cparms[SIZE_CELLPARAMS];
   __shared__ real sha_bparms[SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_cparms[MYIND] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+MYIND);
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND]; 
   
   // Fetch volume avgs to the middle of the shared mem array
   sha_avg[MYIND+  WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND       ];
   sha_avg[MYIND+3*WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_avg[MYIND] = avgnbrx[MYBLOCK*SIZE_BOUND+MYIND+WID2];
   // Fetch 1st derivatives
   sha_d1x[MYIND+  WID2] = d1x[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d1x[MYIND+3*WID2] = d1x[MYBLOCK*SIZE_DERIV+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_d1x[MYIND] = d1xnbr[MYBLOCK*SIZE_BDERI+MYIND];
   // Fetch 2nd derivatives
   sha_d2x[MYIND+  WID2] = d2x[MYBLOCK*SIZE_DERIV+MYIND       ];
   sha_d2x[MYIND+3*WID2] = d2x[MYBLOCK*SIZE_DERIV+MYIND+2*WID2];
   if (threadIdx.z == 0) sha_d2x[MYIND] = d2xnbr[MYBLOCK*SIZE_BDERI+MYIND];
   
   // Transpose volume averages twice, yz + xz:
   transpose_yz_1warp(sha_avg+WID2);
   transpose_xz_1warp(sha_avg+WID2);
   
   // Reconstruct negative and positive side values at the -vx face, and store calculated flux.
   for (int cntr=0; cntr<2; ++cntr) {
      real avg_neg = reconstruct_neg(sha_avg[MYIND     ],sha_d1x[MYIND     ],sha_d2x[MYIND     ]);
      real avg_pos = reconstruct_pos(sha_avg[MYIND+WID2],sha_d1x[MYIND+WID2],sha_d2x[MYIND+WID2]);
      fx[MYBLOCK*SIZE_FLUXS+MYIND] = velocityFluxX(avg_neg,avg_pos,sha_cparms,sha_bparms);
      MYIND += 2*WID2;
   }
}

template<uint WARPS> __global__ void vxFlux_n1warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbrx,real* d1x,real* d2x,real* d1xnbr,real* d2xnbr,real* fx,real* blockParams) {
   cuint SIZE_VEL = SIZE_VELBLOCK+WID2;
   cuint SIZE_DER = SIZE_DERIV+SIZE_BDERI;
   
   __shared__ real sha_avg[WARPS*SIZE_VEL];
   __shared__ real sha_d1x[WARPS*SIZE_DER];
   __shared__ real sha_d2x[WARPS*SIZE_DER];
   __shared__ real sha_cparms[WARPS*2*SIZE_CELLPARAMS];
   __shared__ real sha_bparms[WARPS*2*SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z/2)*gridDim.x + blockIdx.x;
   cuint MYIND1 = tindex3(threadIdx.x,threadIdx.y,threadIdx.z%2);
   cuint MYIND2 = MYIND1 + 2*WID2;
   
   // Load parameters for the spatial cell & velocity block:
   sha_cparms[threadIdx.z*SIZE_CELLPARAMS + threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+threadIdx.y*WID+threadIdx.x);
   sha_bparms[threadIdx.z*SIZE_BLOCKPARAMS+ threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_bparms,MYBLOCK*SIZE_BLOCKPARAMS+threadIdx.y*WID+threadIdx.x);
   
   // Fetch volume avgs to the middle of the shared mem array
   sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND1];
   sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND2];
   if (threadIdx.z % 2 == 0) sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1] = avgnbrx[MYBLOCK*SIZE_BOUND+MYIND1+WID2];
   // Fetch 1st derivatives
   sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2] = d1x[MYBLOCK*SIZE_DERIV+MYIND1];
   sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2] = d1x[MYBLOCK*SIZE_DERIV+MYIND2];
   if (threadIdx.z % 2 == 0) sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND1] = d1xnbr[MYBLOCK*SIZE_BDERI+MYIND1];
   // Fetch 2nd derivatives
   sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2] = d2x[MYBLOCK*SIZE_DERIV+MYIND1];
   sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2] = d2x[MYBLOCK*SIZE_DERIV+MYIND2];
   if (threadIdx.z % 2 == 0) sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND1] = d2xnbr[MYBLOCK*SIZE_BDERI+MYIND1];
   
   transpose_yz_1warp(sha_avg+(threadIdx.z/2)*SIZE_VEL+WID2);
   transpose_xz_1warp(sha_avg+(threadIdx.z/2)*SIZE_VEL+WID2);
   
   // Reconstruct negative and positive side values at the -vx face, and store calculated fluxes:
   real avg_neg,avg_pos;
   avg_neg = reconstruct_neg(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1     ],sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND1     ],sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND1     ]);
   avg_pos = reconstruct_pos(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND1+WID2],sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2],sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND1+WID2]);
   fx[MYBLOCK*SIZE_FLUXS+MYIND1] = velocityFluxX(avg_neg,avg_pos,sha_cparms+threadIdx.z*SIZE_CELLPARAMS,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
   avg_neg = reconstruct_neg(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2     ],sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND2     ],sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND2     ]);
   avg_pos = reconstruct_pos(sha_avg[(threadIdx.z/2)*SIZE_VEL+MYIND2+WID2],sha_d1x[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2],sha_d2x[(threadIdx.z/2)*SIZE_DER+MYIND2+WID2]);
   fx[MYBLOCK*SIZE_FLUXS+MYIND2] = velocityFluxX(avg_neg,avg_pos,sha_cparms+threadIdx.z*SIZE_CELLPARAMS,sha_bparms+threadIdx.z*SIZE_BLOCKPARAMS);
}

__global__ void vxFlux_2warp(uint OFFSET,uint SPATCELL,real* avgs,real* avgnbrx,real* d1x,real* d2x,real* d1xnbr,real* d2xnbr,real* fx,real* blockParams) {
   __shared__ real sha_avg[SIZE_VELBLOCK+WID2];
   __shared__ real sha_d1x[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_d2x[SIZE_DERIV+SIZE_BDERI];
   __shared__ real sha_cparms[SIZE_CELLPARAMS];
   __shared__ real sha_bparms[SIZE_BLOCKPARAMS];
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Load parameters for the spatial cell & velocity block:
   if (threadIdx.z == 0) sha_bparms[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   
   // Fetch volume averages and first+second derivatives for the block:
   sha_avg[MYIND+WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND];
   sha_d1x[MYIND+WID2] = d1x[MYBLOCK*SIZE_DERIV+MYIND];
   sha_d2x[MYIND+WID2] = d2x[MYBLOCK*SIZE_DERIV+MYIND];
   
   // Fetch ghost cells & block parameters:
   if (threadIdx.z == 1) sha_d1x[threadIdx.y*WID+threadIdx.x] = d1xnbr[MYBLOCK*SIZE_BDERI+threadIdx.y*WID+threadIdx.x];
   if (threadIdx.z == 2) sha_d2x[threadIdx.y*WID+threadIdx.x] = d2xnbr[MYBLOCK*SIZE_BDERI+threadIdx.y*WID+threadIdx.x];
   if (threadIdx.z == 3) sha_cparms[threadIdx.y*WID+threadIdx.x] = tex1Dfetch(texRef_cellParams,SPATCELL*SIZE_CELLPARAMS+threadIdx.y*WID+threadIdx.x);
   if (threadIdx.z == 0) sha_avg[threadIdx.y*WID+threadIdx.x] = avgnbrx[MYBLOCK*SIZE_BOUND+WID2+threadIdx.y*WID+threadIdx.x];
   __syncthreads();
   
   transpose_yz_2warp(sha_avg+WID2);
   transpose_xz_2warp(sha_avg+WID2);
   __syncthreads();
   
   // Reconstruct negative and positive side values at the -vx face, and store calculated flux.
   real avg_neg = reconstruct_neg(sha_avg[MYIND     ],sha_d1x[MYIND     ],sha_d2x[MYIND     ]);
   real avg_pos = reconstruct_pos(sha_avg[MYIND+WID2],sha_d1x[MYIND+WID2],sha_d2x[MYIND+WID2]);
   fx[MYBLOCK*SIZE_FLUXS+MYIND] = velocityFluxX(avg_neg,avg_pos,sha_cparms,sha_bparms);
}

__global__ void copyVelFlux_hwarp(uint OFFSET,real* fx,real* fy,real* fz,real* fxnbr,real* fynbr,real* fznbr,uint* nbrs) {
   __shared__ uint sha_nbr[SIZE_NBRS_VEL]; // Shared memory array for neighbour indices.
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   cuint MYIND = tindex2(threadIdx.x,threadIdx.y);
   
   loadVelNbrs(MYBLOCK,nbrs,sha_nbr);

   // If the block is a boundary block calculate vx-flux using a boundary function, otherwise copy from +vx neighbour:
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VX_POS_BND) == NbrsVel::VX_POS_BND) {
      //fxnbr[MYBLOCK*SIZE_BFLUX + MYIND] = 0.0f;
      fxnbr[MYBLOCK*SIZE_BFLUX + MYIND] = fx[MYBLOCK*SIZE_FLUXS + 3*WID2 + MYIND];
   } else {
      fxnbr[MYBLOCK*SIZE_BFLUX + MYIND] = fx[sha_nbr[NbrsVel::VXPOS]*SIZE_FLUXS + MYIND];
   }
   // If the block is a boundary block calculate vy-flux using a boundary function, otherwise copy from +vy neighbour:
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VY_POS_BND) == NbrsVel::VY_POS_BND) {
      //fynbr[MYBLOCK*SIZE_BFLUX + MYIND] = 0.0f;
      fynbr[MYBLOCK*SIZE_BFLUX + MYIND] = fy[MYBLOCK*SIZE_FLUXS + 3*WID2 + MYIND];
   } else {
      fynbr[MYBLOCK*SIZE_BFLUX + MYIND] = fy[sha_nbr[NbrsVel::VYPOS]*SIZE_FLUXS + MYIND];
   }
   // If the block is a boundary block calculate vz-flux using a boundary function, otherwise copy from +vz neighbour:
   if ((sha_nbr[NbrsVel::STATE] & NbrsVel::VZ_POS_BND) == NbrsVel::VZ_POS_BND) {
      //fznbr[MYBLOCK*SIZE_BFLUX + MYIND] = 0.0f;
      fznbr[MYBLOCK*SIZE_BFLUX + MYIND] = fz[MYBLOCK*SIZE_FLUXS + 3*WID2 + MYIND];
   } else {
      fznbr[MYBLOCK*SIZE_BFLUX + MYIND] = fz[sha_nbr[NbrsVel::VZPOS]*SIZE_FLUXS + MYIND];
   }
}

template<uint WARPS> __global__ void copyVelFlux_nhwarp(uint OFFSET,real* fx,real* fy,real* fz,real* fxnbr,real* fynbr,real* fznbr,uint* nbrs) {
   //__shared__ uint sha_nbr[WARPS*SIZE_NBRS_VEL]; // Shared memory array for neighbour indices.
   cuint SIZE_NBRS = SIZE_NBRS_VEL+1;
   __shared__ uint sha_nbr[WARPS*SIZE_NBRS];
   
   //cuint MYBLOCK = OFFSET + (blockIdx.y*WARPS + threadIdx.z)*gridDim.x + blockIdx.x;
   cuint MYBLOCK = OFFSET + blockIdx.y*gridDim.x*WARPS + blockIdx.x*WARPS + threadIdx.z;
   cuint MYIND = tindex2(threadIdx.x,threadIdx.y);
   
   sha_nbr[threadIdx.z*SIZE_NBRS + threadIdx.y*WID + threadIdx.x] = nbrs[MYBLOCK*SIZE_NBRS_VEL + threadIdx.y*WID + threadIdx.x];

   // If the block is a boundary block calculate vx-flux using a boundary function, otherwise copy from +vx neighbour:
   if ((sha_nbr[threadIdx.z*SIZE_NBRS+NbrsVel::STATE] & NbrsVel::VX_POS_BND) == NbrsVel::VX_POS_BND) {
      fxnbr[MYBLOCK*SIZE_BFLUX + MYIND] = fx[MYBLOCK*SIZE_FLUXS + 3*WID2 + MYIND];
   } else {
      fxnbr[MYBLOCK*SIZE_BFLUX + MYIND] = fx[sha_nbr[threadIdx.z*SIZE_NBRS+NbrsVel::VXPOS]*SIZE_FLUXS + MYIND];
   }
   // If the block is a boundary block calculate vy-flux using a boundary function, otherwise copy from +vy neighbour:
   if ((sha_nbr[threadIdx.z*SIZE_NBRS+NbrsVel::STATE] & NbrsVel::VY_POS_BND) == NbrsVel::VY_POS_BND) {
      fynbr[MYBLOCK*SIZE_BFLUX + MYIND] = fy[MYBLOCK*SIZE_FLUXS + 3*WID2 + MYIND];
   } else {
      fynbr[MYBLOCK*SIZE_BFLUX + MYIND] = fy[sha_nbr[threadIdx.z*SIZE_NBRS+NbrsVel::VYPOS]*SIZE_FLUXS + MYIND];
   }
   // If the block is a boundary block calculate vz-flux using a boundary function, otherwise copy from +vz neighbour:
   if ((sha_nbr[threadIdx.z*SIZE_NBRS+NbrsVel::STATE] & NbrsVel::VZ_POS_BND) == NbrsVel::VZ_POS_BND) {
      fznbr[MYBLOCK*SIZE_BFLUX + MYIND] = fz[MYBLOCK*SIZE_FLUXS + 3*WID2 + MYIND];
   } else {
      fznbr[MYBLOCK*SIZE_BFLUX + MYIND] = fz[sha_nbr[threadIdx.z*SIZE_NBRS+NbrsVel::VZPOS]*SIZE_FLUXS + MYIND];
   }
}

// ---------------------------------------------------------------------------------
// ---------- FUNCTIONS FOR PROPAGATING VOLUME AVERAGES IN VELOCITY SPACE ----------
// ---------------------------------------------------------------------------------

__global__ void propagateVel_1warp(uint OFFSET,uint SPATCELL,real* avgs,real* fx,real* fy,real* fz,
				       real* fxnbr,real* fynbr,real* fznbr,real* blockParams,real DT) {
   __shared__ real sha_avg[SIZE_VELBLOCK];          // Shared mem array for volume averages
   __shared__ real sha_bparams[SIZE_BLOCKPARAMS];   //                      block parameters
   __shared__ real sha_flux[SIZE_FLUXS+SIZE_BFLUX]; //                      flux
   
   cuint MYBLOCK = OFFSET + bindex2(blockIdx.x,blockIdx.y);
   uint MYIND = tindex3(threadIdx.x,threadIdx.y,threadIdx.z);
   
   // Fetch parameters for the block and grid
   if (threadIdx.z == 0) sha_bparams[MYIND] = blockParams[MYBLOCK*SIZE_BLOCKPARAMS+MYIND];
   
   // Fetch volume averages:
   sha_avg[MYIND       ] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND       ];
   sha_avg[MYIND+2*WID2] = avgs[MYBLOCK*SIZE_VELBLOCK+MYIND+2*WID2];
   
   // Fetch vz-flux and calculate its contribution to the volume average:
   sha_flux[MYIND       ] = fz[MYBLOCK*SIZE_FLUXS+MYIND       ];
   sha_flux[MYIND+2*WID2] = fz[MYBLOCK*SIZE_FLUXS+MYIND+2*WID2];
   sha_flux[threadIdx.y*WID+threadIdx.x+WID3] = fznbr[MYBLOCK*SIZE_BFLUX+threadIdx.y*WID+threadIdx.x];
   sha_avg[MYIND       ] += (sha_flux[MYIND       ] - sha_flux[MYIND+  WID2])*DT/sha_bparams[BlockParams::DVZ];
   sha_avg[MYIND+2*WID2] += (sha_flux[MYIND+2*WID2] - sha_flux[MYIND+3*WID2])*DT/sha_bparams[BlockParams::DVZ];
   
   // Transpose in yz, fetch vy-flux, and calculate its contribution to the volume average:
   transpose_yz_1warp(sha_avg);
   sha_flux[MYIND       ] = fy[MYBLOCK*SIZE_FLUXS+MYIND       ];
   sha_flux[MYIND+2*WID2] = fy[MYBLOCK*SIZE_FLUXS+MYIND+2*WID2];
   sha_flux[threadIdx.y*WID+threadIdx.x+WID3] = fynbr[MYBLOCK*SIZE_BFLUX+threadIdx.y*WID+threadIdx.x];
   sha_avg[MYIND       ] += (sha_flux[MYIND       ] - sha_flux[MYIND+  WID2])*DT/sha_bparams[BlockParams::DVY];
   sha_avg[MYIND+2*WID2] += (sha_flux[MYIND+2*WID2] - sha_flux[MYIND+3*WID2])*DT/sha_bparams[BlockParams::DVY];
   
   // Transpose in xz, fetch vx-flux, and calculate its contribution to the volume average:
   transpose_xz_1warp(sha_avg);
   sha_flux[MYIND       ] = fx[MYBLOCK*SIZE_FLUXS+MYIND       ];
   sha_flux[MYIND+2*WID2] = fx[MYBLOCK*SIZE_FLUXS+MYIND+2*WID2];
   sha_flux[threadIdx.y*WID+threadIdx.x+WID3] = fxnbr[MYBLOCK*SIZE_BFLUX+threadIdx.y*WID+threadIdx.x];
   sha_avg[MYIND       ] += (sha_flux[MYIND       ] - sha_flux[MYIND+  WID2])*DT/sha_bparams[BlockParams::DVX];
   sha_avg[MYIND+2*WID2] += (sha_flux[MYIND+2*WID2] - sha_flux[MYIND+3*WID2])*DT/sha_bparams[BlockParams::DVX];
   
   // Invert transposes and store results:
   transpose_xz_1warp(sha_avg);
   transpose_yz_1warp(sha_avg);
   avgs[MYBLOCK*SIZE_VELBLOCK+MYIND       ] = sha_avg[MYIND       ];
   avgs[MYBLOCK*SIZE_VELBLOCK+MYIND+2*WID2] = sha_avg[MYIND+2*WID2];
}

/** Copy volume averages to ghost cells of each velocity grid block using a 
 * method selected at compile time. Different methods can be used to optimize 
 * CUDA performance. The method is selected by setting one (and only one) of 
 * the following preprocessor macros in cudalaunch.h: 
 * VELCOPY_1WARP
 * 
 * Additionally the x-size of CUDA grid is set with the parameter VELCOPY_GSIZEX.
 */
inline void copyVelAverages(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef VELCOPY_1WARP
   gridSize.x = VELCOPY_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 2;
   copyVelBlocks_1warp<<<gridSize,blockSizes>>>(offset,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),
						deviceGrid.getAvgnbry(),deviceGrid.getAvgnbrz(),deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   
   // Copy the leftover blocks:
   #ifdef EXTRA_VELCOPY_1WARP
   blockSizes.z = 2;
   copyVelBlocks_1warp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),
						     deviceGrid.getAvgnbry(),deviceGrid.getAvgnbrz(),deviceGrid.getNbrsVel());
   #endif
}

/** Calculate the derivatives of volume averages in velocity coordinates using a 
 * method selected at compile time. Different methods can be used to optimize CUDA 
 * performance. The method is selected by setting one (and only one) of the 
 * following preprocessor macros in cudalaunch.h:
 * VEL_DERIVS_1WARP
 * 
 * Additionally the x-size of CUDA grid is set with the parameter VEL_DERIVS_GSIZEX.
 */
inline void calcVelDerivatives(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef VEL_DERIVS_1WARP // 1 warp calculates one velocity grid block
   gridSize.x = VEL_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 2;
   calcVelDerivs_1warp<<<gridSize,blockSizes>>>(offset,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),deviceGrid.getAvgnbry(),
						deviceGrid.getAvgnbrz(),deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z());
   offset += gridSize.x*gridSize.y;
   #endif
   
   #ifdef EXTRA_VEL_DERIVS_1WARP
   blockSizes.z = 2;
   calcVelDerivs_1warp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),deviceGrid.getAvgnbry(),
						     deviceGrid.getAvgnbrz(),deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						     deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z());
   #endif
}

/** Copy derivatives of volume averages to ghost cells of each velocity grid block
 * using a method selected at compile time. Different methods can be used to optimize 
 * CUDA performace. The method is selected by setting one (and only one) of the 
 * following preprocessor macros in cudalaunch.h:
 * COPY_VEL_DERIVS_1WARP
 * 
 * Additionally the x-size of CUDA grid is set with parameter COPY_VEL_DERIVS_GSIZEX.
 */
inline void copyVelDerivatives(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef COPY_VEL_DERIVS_HWARP // Half-warp calculates one velocity grid block
   gridSize.x = COPY_VEL_DERIVS_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 1;
   copyVelDerivs_hwarp<<<gridSize,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						deviceGrid.getD1xnbr(),deviceGrid.getD1ynbr(),deviceGrid.getD1znbr(),
						deviceGrid.getNbrsVel());
   copyVelDerivs_hwarp<<<gridSize,blockSizes>>>(offset,deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),
						deviceGrid.getD2xnbr(),deviceGrid.getD2ynbr(),deviceGrid.getD2znbr(),
						deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   
   #ifdef EXTRA_COPY_VEL_DERIVS_HWARP
   blockSizes.z = 1;
   copyVelDerivs_hwarp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getD1x(),deviceGrid.getD1y(),deviceGrid.getD1z(),
						     deviceGrid.getD1xnbr(),deviceGrid.getD1ynbr(),deviceGrid.getD1znbr(),
						     deviceGrid.getNbrsVel());
   
   copyVelDerivs_hwarp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getD2x(),deviceGrid.getD2y(),deviceGrid.getD2z(),
						     deviceGrid.getD2xnbr(),deviceGrid.getD2ynbr(),deviceGrid.getD2znbr(),
						     deviceGrid.getNbrsVel());
   #endif
}

/** Calculate z-components of velocity fluxes using a method selected at 
 * compile time. These methods are used to optimise CUDA performance. 
 * The method is selected by setting one (and only one) of the following 
 * preprocessor macros in cudalaunch.h: VELFLUX_1WARP, VELFLUX_N1WARP, VELFLUX_2WARP.
 *
 * The difference 
 * between these methods is the number of warps used to calculate the fluxes for one 
 * velocity grid block. VELFLUX_1WARP uses one warp per block, VELFLUX_2warp uses two 
 * warps per block, and VELFLUX_N1WARP is a templatized function which uses one warp per 
 * block and simultaneously calculates N blocks.
 * 
 * Additionally the size of CUDA grid is selected by setting the desired value for 
 * parameter VELFLUX_GSIZEX, and if VELFLUX_N1WARP is used the the number of blocks 
 * calculated simultaneously is selected by setting an integer value >= 1 for parameter
 * VELFLUX_WARPS.
 */
inline void calcFluxVZ(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef VELFLUX_1WARP // 1 warp calculates one velocity grid block
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 2;
   vzFlux_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrz(),
					 deviceGrid.getD1z(),deviceGrid.getD2z(),deviceGrid.getD1znbr(),deviceGrid.getD2znbr(),
					 deviceGrid.getFz(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef VELFLUX_2WARP // 2 warps calculate one velocity grid block
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 4;
   vzFlux_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrz(),
					 deviceGrid.getD1z(),deviceGrid.getD2z(),deviceGrid.getD1znbr(),deviceGrid.getD2znbr(),
					 deviceGrid.getFz(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef VELFLUX_N1WARP // N 1-warps calculate one velocity grid block each
   blockSizes.z = 2*VELFLUX_WARPS;
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x * VELFLUX_WARPS);
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*VELFLUX_WARPS;
   gridSizeExtra.y = 1;
   vzFlux_n1warp<VELFLUX_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrz(),
							 deviceGrid.getD1z(),deviceGrid.getD2z(),deviceGrid.getD1znbr(),deviceGrid.getD2znbr(),
							 deviceGrid.getFz(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y*VELFLUX_WARPS;
   #endif
   
   #ifdef EXTRA_VELFLUX_1WARP
   blockSizes.z = 2;
   vzFlux_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrz(),
					      deviceGrid.getD1z(),deviceGrid.getD2z(),deviceGrid.getD1znbr(),deviceGrid.getD2znbr(),
					      deviceGrid.getFz(),deviceGrid.getBlockParams());
   #endif
}

inline void calcFluxVY(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef VELFLUX_1WARP // 1 warp calculates one velocity grid block
   gridSize.x = 8;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   vyFlux_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbry(),
					 deviceGrid.getD1y(),deviceGrid.getD2y(),deviceGrid.getD1ynbr(),deviceGrid.getD2ynbr(),
					 deviceGrid.getFy(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef VELFLUX_2WARP // 2 warps calculate one velocity grid block
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 4;
   vyFlux_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbry(),
					 deviceGrid.getD1y(),deviceGrid.getD2y(),deviceGrid.getD1ynbr(),deviceGrid.getD2ynbr(),
					 deviceGrid.getFy(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef VELFLUX_N1WARP // N 1-warps calculate one velocity grid block each
   blockSizes.z = 2*VELFLUX_WARPS;
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x * VELFLUX_WARPS);
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*VELFLUX_WARPS;
   gridSizeExtra.y = 1;
   vyFlux_n1warp<VELFLUX_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbry(),
							 deviceGrid.getD1y(),deviceGrid.getD2y(),deviceGrid.getD1ynbr(),deviceGrid.getD2ynbr(),
							 deviceGrid.getFy(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y*VELFLUX_WARPS;
   #endif
   
   #ifdef EXTRA_VELFLUX_1WARP
   blockSizes.z = 2;
   vyFlux_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbry(),
					      deviceGrid.getD1y(),deviceGrid.getD2y(),deviceGrid.getD1ynbr(),deviceGrid.getD2ynbr(),
					      deviceGrid.getFy(),deviceGrid.getBlockParams());
   #endif
}

inline void calcFluxVX(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef VELFLUX_1WARP // 1 warp calculates one velocity grid block
   gridSize.x = 8;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   vxFlux_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),
					 deviceGrid.getD1x(),deviceGrid.getD2x(),deviceGrid.getD1xnbr(),deviceGrid.getD2xnbr(),
					 deviceGrid.getFx(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef VELFLUX_2WARP // 2 warps calculate one velocity grid block
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   blockSizes.z = 4;
   vxFlux_2warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),
					 deviceGrid.getD1x(),deviceGrid.getD2x(),deviceGrid.getD1xnbr(),deviceGrid.getD2xnbr(),
					 deviceGrid.getFx(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef VELFLUX_N1WARP // N 1-warps calculate one velocity grid block each
   blockSizes.z = 2*VELFLUX_WARPS;
   gridSize.x = VELFLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x * VELFLUX_WARPS);
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*VELFLUX_WARPS;
   gridSizeExtra.y = 1;
   vxFlux_n1warp<VELFLUX_WARPS><<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),
							 deviceGrid.getD1x(),deviceGrid.getD2x(),deviceGrid.getD1xnbr(),deviceGrid.getD2xnbr(),
							 deviceGrid.getFx(),deviceGrid.getBlockParams());
   offset += gridSize.x*gridSize.y*VELFLUX_WARPS;
   #endif
   
   #ifdef EXTRA_VELFLUX_1WARP
   blockSizes.z = 2;
   vxFlux_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),deviceGrid.getAvgnbrx(),
					      deviceGrid.getD1x(),deviceGrid.getD2x(),deviceGrid.getD1xnbr(),deviceGrid.getD2xnbr(),
					      deviceGrid.getFx(),deviceGrid.getBlockParams());
   #endif
}

/** Copy fluxes to ghost cells of each velocity grid block using a method 
 * selected at compile time. Different methods can be used to optimize CUDA 
 * performance. The method is selected by setting one (and only one) of the 
 * following preprocessor macros in cudalaunch.h:
 * COPY_VEL_FLUX_HWARP, COPY_VEL_FLUX_NHWARP
 * 
 * Additionally the x-size of CUDA grid is set with parameter COPY_VEL_FLUX_GSIZEX. 
 * If COPY_VEL_FLUX_NHWARP is used, then the number or blocks calculated simultaneously 
 * must be speciefied by setting an integer value >=1 for parameter COPY_VEL_FLUX_WARPS.
 */
inline void copyVelFluxes(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes) {
   offset = grid[cell]->velBlockIndex;
   #ifdef COPY_VEL_FLUX_HWARP // Half-warp calculates one velocity grid block
   gridSize.x = COPY_VEL_FLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;

   blockSizes.z = 1;
   copyVelFlux_hwarp<<<gridSize,blockSizes>>>(offset,deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
					      deviceGrid.getFxnbr(),deviceGrid.getFynbr(),deviceGrid.getFznbr(),
					      deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y;
   #endif
   #ifdef COPY_VEL_FLUX_NHWARP // N half-warps calculate one velcoity grid block each
   gridSize.x = COPY_VEL_FLUX_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / (gridSize.x * COPY_VEL_FLUX_WARPS);
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y*COPY_VEL_FLUX_WARPS;
   gridSizeExtra.y = 1;
   blockSizes.z = COPY_VEL_FLUX_WARPS;
   copyVelFlux_nhwarp<COPY_VEL_FLUX_WARPS><<<gridSize,blockSizes>>>(offset,deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
								    deviceGrid.getFxnbr(),deviceGrid.getFynbr(),deviceGrid.getFznbr(),
								    deviceGrid.getNbrsVel());
   offset += gridSize.x*gridSize.y*COPY_VEL_FLUX_WARPS;
   #endif
         
   #ifdef EXTRA_COPY_VEL_FLUX_HWARP
   blockSizes.z = 1;
   copyVelFlux_hwarp<<<gridSizeExtra,blockSizes>>>(offset,deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
						   deviceGrid.getFxnbr(),deviceGrid.getFynbr(),deviceGrid.getFznbr(),
						   deviceGrid.getNbrsVel());
   #endif
}

inline void velPropagate(Grid& grid,DeviceGrid& deviceGrid,cuint& cell,uint& offset,dim3& gridSize,dim3& gridSizeExtra,dim3& blockSizes,creal& DT) {
   offset = grid[cell]->velBlockIndex;
   #ifdef VELPROP_1WARP // 1 warp calculates one velocity grid block
   gridSize.x = VELPROP_GSIZEX;
   gridSize.y = grid[cell]->N_blocks / gridSize.x;
   gridSize.z = 1;
   gridSizeExtra.x = grid[cell]->N_blocks - gridSize.x*gridSize.y;
   gridSizeExtra.y = 1;
   gridSizeExtra.z = 1;
   blockSizes.z = 2;
   propagateVel_1warp<<<gridSize,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
					       deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
					       deviceGrid.getFxnbr(),deviceGrid.getFynbr(),deviceGrid.getFznbr(),
					       deviceGrid.getBlockParams(),DT);
   
   offset += gridSize.x*gridSize.y;
   #endif
   
   #ifdef EXTRA_VELPROP_1WARP
   blockSizes.z = 2;
   propagateVel_1warp<<<gridSizeExtra,blockSizes>>>(offset,grid[cell]->cellIndex,deviceGrid.getBlockArray(),
						    deviceGrid.getFx(),deviceGrid.getFy(),deviceGrid.getFz(),
						    deviceGrid.getFxnbr(),deviceGrid.getFynbr(),deviceGrid.getFznbr(),
						    deviceGrid.getBlockParams(),DT);
   #endif
}

bool acceleration(Grid& grid,DeviceGrid& deviceGrid,creal& DT) {
   bool success = true;
   
   // Bind textures:
   size_t texOffset;
   if (bindTextureAcc(texRef_avgs,deviceGrid.getBlockArray(),deviceGrid.sizeBlockArray(),texOffset) == false) success = false;
   if (bindTextureAcc(texRef_cellParams,deviceGrid.getCellParams(),deviceGrid.sizeCellParams(),texOffset) == false) success = false;
   if (bindTextureAcc(texRef_bparms,deviceGrid.getBlockParams(),deviceGrid.sizeBlockParams(),texOffset) == false) success = false;
   if (success == false) return success;
   
   // Calculate acceleration for each spatial cell
   uint offset;
   dim3 gridSize;
   dim3 gridSizeExtra;
   dim3 blockSizes(WID,WID,2);
   for (uint cell=0; cell<grid.size(); ++cell) {
      // Copy volume averages to ghost cells:
      copyVelAverages(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
      
      // Calculate derivatives in velocity coordinates:
      calcVelDerivatives(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
      
      // Copy derivatives to ghost cells:
      copyVelDerivatives(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);

      // Calculate vx,vy,vz -fluxes:
      calcFluxVZ(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
      calcFluxVY(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);
      calcFluxVX(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);

      // Copy vx,vy,vz -fluxes to ghost cells:
      copyVelFluxes(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes);

      // Propagate volume averages in velocity space:
      velPropagate(grid,deviceGrid,cell,offset,gridSize,gridSizeExtra,blockSizes,DT);

   }
   return success;
}





