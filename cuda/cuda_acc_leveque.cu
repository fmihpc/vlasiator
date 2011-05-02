#include <cuda_runtime.h>

#include "../common.h"

const Real EPSILON = 1.0e-25;
const Real ZERO = 0.0;
const Real FOURTH = 1.0/4.0;
const Real SIXTH = 1.0/6.0;
const Real HALF = 0.5;
const Real ONE  = 1.0;
const Real TWO  = 2.0;

const Real AX = 0.0f;
const Real AY = 0.0f;
const Real AZ = 0.0f;

texture<Real,1,cudaReadModeElementType> texRef_avgs;
texture<Real,1,cudaReadModeElementType> texRef_dfdt;
texture<Real,1,cudaReadModeElementType> texRef_flux;
texture<uint,1,cudaReadModeElementType> texRef_nbrsVel;
texture<Real,1,cudaReadModeElementType> texRef_blockParams;
texture<Real,1,cudaReadModeElementType> texRef_cellParams;

texture<Real,2,cudaReadModeElementType> texRef_avgs2D;
texture<Real,2,cudaReadModeElementType> texRef_dfdt2D;
texture<Real,2,cudaReadModeElementType> texRef_flux2D;

__device__ uint blockIndex(uint bix,uint biy) {return biy*gridDim.x + bix;}
__device__ uint index(uint i,uint j,uint k) {return k*WID2 + j*WID + i;}

__device__ int neg_sign(int i) {return (i >> 31);}
__device__ int pos_sign(int i) {return i/WID;}

__device__ int sgn(int i) {
   int rvalue = (i >> 31);
   rvalue += max(0,i)/WID;
   return rvalue;
}

__device__ Real fetchAvgs(uint MYBLOCK,uint offsetAvgs,uint offsetNbrs,uint STATUS,int i,int j,int k) {
   uint blockOffset = 13 + sgn(k)*9 + sgn(j)*3 + sgn(i);
   blockOffset = offsetAvgs + tex1Dfetch(texRef_nbrsVel,offsetNbrs + MYBLOCK*28 + blockOffset)*WID3;
   blockOffset += index(i-sgn(i)*WID,j-sgn(j)*WID,k-sgn(k)*WID);
   //return tex1Dfetch(texRef_avgs,blockOffset);
   
   // 2D texture fetches are required for large arrays:
   cuint YCRD = blockOffset / CUDA_WIDTH;
   cuint XCRD = blockOffset - YCRD*CUDA_WIDTH;
   return tex2D(texRef_avgs2D,XCRD+HALF,YCRD+HALF);
}

__device__ Real fetchFlux(uint MYBLOCK,uint offsetFlux,uint offsetNbrs,uint STATUS,int i,int j,int k) {
   
   uint blockOffset = 13 + sgn(k)*9 + sgn(j)*3 + sgn(i);
   blockOffset = offsetFlux + tex1Dfetch(texRef_nbrsVel,offsetNbrs + MYBLOCK*28 + blockOffset)*WID3;
   blockOffset += index(i-sgn(i)*WID,j-sgn(j)*WID,k-sgn(k)*WID);
   //return tex1Dfetch(texRef_avgs,blockOffset);
   
   // 2D texture fetches are requiref for large arrays:
   cuint YCRD = blockOffset / CUDA_WIDTH;
   cuint XCRD = blockOffset - YCRD*CUDA_WIDTH;
   return tex2D(texRef_flux2D,XCRD+HALF,YCRD+HALF);
}

__device__ Real accelerationX(Real Vx,Real Vy,Real Vz,Real* cellParams,Real q_per_m) {
   return q_per_m*(cellParams[CellParams::EX] + Vy*cellParams[CellParams::BZ] - Vz*cellParams[CellParams::BY]);
}

__device__ Real accelerationY(Real Vx,Real Vy,Real Vz,Real* cellParams,Real q_per_m) {
   return q_per_m*(cellParams[CellParams::EY] + Vz*cellParams[CellParams::BX] - Vx*cellParams[CellParams::BZ]);
}

__device__ Real accelerationZ(Real Vx,Real Vy,Real Vz,Real* cellParams,Real q_per_m) {
   return q_per_m*(cellParams[CellParams::EZ] + Vx*cellParams[CellParams::BY] - Vy*cellParams[CellParams::BX]);
}

__device__ Real VanLeer(Real theta) {
   return (theta + fabs(theta)) / (ONE + fabs(theta));
}

__device__ Real MClimiter(Real theta) {
   Real phi = fmin(TWO,TWO*theta);
   phi = fmin(HALF*(ONE+theta),phi);
   return fmax(ZERO,phi);
}

__device__ Real superbee(Real theta) {
   Real phi = fmin(TWO,theta);
   phi = fmax(phi,fmin(ONE,TWO*theta));
   return fmax(ZERO,phi);
}

__device__ Real limiter(Real theta) {
   //return VanLeer(theta);
   //return MClimiter(theta);
   return superbee(theta);
}

__global__ void cuda_acc_xfluxes(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   volatile Real Vx,Vy,Vz,R;
   volatile Real Ax,Ay,Az;
   Real myflux,dt_per_dvx;
   cuint STATUS = tex1Dfetch(texRef_nbrsVel,offsetNbrs+MYBLOCK*28+27);
   
   // Fetch block parameters into shared mem array (causes warp serialization):
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   // Fetch cell parameters to shared mem array:
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   // Fetch the 4-point stencil of volume averages for this thread:
   Real xm2,xm1,xcc,xp1;

   // ******************************************************************
   // ***** Calculate flux at interface between (i,j,k), (i-1,j,k) *****
   // ******************************************************************

   // ***** Add contributions from vx-faces *****
   dt_per_dvx = dt/blockParams[BlockParams::DVX];
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   
   // Contribution from Ax > 0:
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z);
   myflux = Ax*xm1;
 
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y,threadIdx.z);
   R = (xcc - xm1)*limiter((xm1-xm2)/(xcc-xm1+EPSILON));
   R = HALF*Ax*(ONE - dt_per_dvx*Ax)*R;
   myflux += R;

   // Contribution from Ax < 0:
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*xcc;
   
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y,threadIdx.z);
   R = (xcc - xm1)*limiter((xp1-xcc)/(xcc-xm1+EPSILON));
   R = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R;
   myflux += R;

   // ***** Add contributions from vy-faces *****
   // correction wave = 0.5*Ax*|Ay|(1-dt/dvy*|Ay|)*dt/dvy*R(limited)
   dt_per_dvx = dt/blockParams[BlockParams::DVY];
   
   // ***** Case Ax > 0, Ay > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z);    
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvx*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-2,threadIdx.z);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE - Ay*dt_per_dvx)*R;

   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE - Ay*dt_per_dvx)*R;
   
   // ***** Case Ax < 0, Ay > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z);   
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvx*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-2,threadIdx.z);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE - Ay*dt_per_dvx)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE - Ay*dt_per_dvx)*R;
   
   // ***** Case Ax > 0, Ay < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));

   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z);
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvx*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+2,threadIdx.z);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE + Ay*dt_per_dvx)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE + Ay*dt_per_dvx)*R;
   
   // ***** Case Ax < 0, Ay < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z);   
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvx*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // 2D Correction Waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+2,threadIdx.z);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE + Ay*dt_per_dvx)*R;

   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Ay*dt_per_dvx*(FOURTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*(ONE + Ay*dt_per_dvx)*R;
   
   // ***** Add contributions from vz-faces *****
   // correction wave = 0.5*Ax*|Az|(1-dt/dvy*|Az|)*dt/dvz*R(limited)
   dt_per_dvx = dt/blockParams[BlockParams::DVZ];

   // ***** Case Ax > 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvx*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // 2D Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE - Az*dt_per_dvx)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE - Az*dt_per_dvx)*R;
   
   // ***** Case Ax < 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));

   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvx*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // 2D Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE - Az*dt_per_dvx)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE - Az*dt_per_dvx)*R;
   
   // ***** Case Ax > 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvx*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE + Az*dt_per_dvx)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));   
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE + Az*dt_per_dvx)*R;
   
   // ***** Case Ax < 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvx*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE + Az*dt_per_dvx)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ax*Az*dt_per_dvx*(FOURTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*(ONE + Az*dt_per_dvx)*R;

   // *********************************
   // ***** Store calculated flux *****
   // *********************************
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] = myflux;
}

__global__ void cuda_acc_xfluxes3DcorrY(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   volatile Real Vx,Vy,Vz;
   Real Ax,Ay,Az,R,dt_per_dvx;
   Real myflux = 0.0f;
   cuint STATUS = 0;
   
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   Real xm2,xm1,xcc,xp1;
   
   dt_per_dvx = dt/blockParams[BlockParams::DVY];
   
   // Ax > 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-2,threadIdx.z-1);
   R = xcc-xm1;
   
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   // Ax > 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-2,threadIdx.z+1);
   R = xcc-xm1;
   
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   // Ax < 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-2,threadIdx.z-1);
   R = xcc-xm1;
   
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   // Ax < 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-2,threadIdx.z+1);
   R = xcc-xm1;
   
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvx*Ay)*R;
   
   // Ax > 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+2,threadIdx.z-1);
   R = xcc-xm1;
   
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + dt_per_dvx*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + dt_per_dvx*Ay)*R;
   
   // Ax > 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+1);
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+2,threadIdx.z+1);
   R = xcc-xm1;
   
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + dt_per_dvx*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + dt_per_dvx*Ay)*R;
   
   // Ax < 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z-1);
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+2,threadIdx.z-1);
   R = xcc-xm1;
   
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + Ay*dt_per_dvx)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + Ay*dt_per_dvx)*R;
   
   // Ax < 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z+1);
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+2,threadIdx.z+1);
   R = xcc-xm1;
   
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + Ay*dt_per_dvx)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVZ]*(ONE + Ay*dt_per_dvx)*R;

   // Store
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] += myflux;
}

__global__ void cuda_acc_xfluxes3DcorrZ(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   volatile Real Vx,Vy,Vz;
   Real Ax,Ay,Az,R,dt_per_dvx;
   Real myflux = 0.0f;
   cuint STATUS = 0;
   
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
      i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   Real xm2,xm1,xcc,xp1;
   
   dt_per_dvx = dt/blockParams[BlockParams::DVZ];
   
   // Ax > 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-2);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   // Ax > 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-2);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   // Ax < 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-2);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   // Ax < 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-1);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-2);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE - dt_per_dvx*Az)*R;
   
   // Ax > 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+2);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   // Ax > 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z  );
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z+2);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   // Ax < 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z  );
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+2);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   // Ax < 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z  );
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z+2);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*R;
   
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvx*dt/blockParams[BlockParams::DVY]*(ONE + dt_per_dvx*Az)*R;
   
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] += myflux;
}


__global__ void cuda_acc_yfluxes(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   Real Vx,Vy,Vz,Ax,Ay,Az,R,myflux,dt_per_dvy;
   cuint STATUS = tex1Dfetch(texRef_nbrsVel,offsetNbrs+MYBLOCK*28+27);
   
   // Fetch block parameters into shared mem array:
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   // Fetch cell parameters to shared mem array:
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   // Fetch the 4-point stencil of volume averages for this thread:
   Real xm2,xm1,xcc,xp1;
   
   // ******************************************************************
   // ***** Calculate flux at interface between (i,j,k), (i,j-1,k) *****
   // ******************************************************************
    
   // ***** Add contributions from vy-faces *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   dt_per_dvy = dt/blockParams[BlockParams::DVY];
   
   // Contribution from Ay > 0:
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z);
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux = Ay*xm1;

   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z);
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-2,threadIdx.z);
   R = (xcc - xm1)*limiter((xm1-xm2)/(xcc-xm1+EPSILON));
   myflux += R*HALF*Ay*(ONE - dt_per_dvy*Ay);

   // Contribution from Ay < 0:
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Ay*xcc;

   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+1,threadIdx.z);
   R = (xcc - xm1)*limiter((xp1-xcc)/(xcc-xm1+EPSILON));
   myflux -= R*HALF*Ay*(ONE + dt_per_dvy*Ay);

   // ***** Add Contributions from Vx-faces *****
   // correction wave = 0.5*|Ax|*Ay*(1-dt/dvx*|Ax|)*dt/dvx*R(limited)
   dt_per_dvy = dt/blockParams[BlockParams::DVX];
   
   // ***** Case Ax > 0, Ay > 0 *****
   dt_per_dvy = dt/blockParams[BlockParams::DVX];
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z);   
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvy*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y-1,threadIdx.z);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ax*Ay*(ONE-dt_per_dvy*Ax)*R*dt_per_dvy;

   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Ay*(ONE-dt_per_dvy*Ax)*R*dt_per_dvy;

   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   // ***** Case Ax > 0, Ay < 0 ***** 
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z);   
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvy*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // 2D Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y  ,threadIdx.z);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ax*Ay*(ONE-dt_per_dvy*Ax)*R*dt_per_dvy;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R; // 3D
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Ay*(ONE-dt_per_dvy*Ax)*R*dt_per_dvy;

   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R; // 3D
   
   // ***** Case Ax < 0, Ay > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z);
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvy*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y-1,threadIdx.z);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ax*Ay*(ONE+dt_per_dvy*Ax)*R*dt_per_dvy;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R; // 3D
      
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Ay*(ONE+dt_per_dvy*Ax)*R*dt_per_dvy;

   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R; // 3D
   
   // ***** Case Ax < 0, Ay < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z);
   R = xcc-xm1;
   myflux += Ax*Ay*dt_per_dvy*(SIXTH*Az*dt/blockParams[BlockParams::DVZ] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y  ,threadIdx.z);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ax*Ay*(ONE+dt_per_dvy*Ax)*R*dt_per_dvy;

   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R; // 3D
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Ay*(ONE+dt_per_dvy*Ax)*R*dt_per_dvy;
 
   Az = fabs(accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R; // 3D
   
   // ***** Add contributions from vz-faces *****
   // correction wave = 0.5*|Az|*Ay*(1-dt/dvz*|Az|)*dt/dvz*R(limited)
   dt_per_dvy = dt / blockParams[BlockParams::DVZ];
   
   // ***** Case Ay > 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvy*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction wave:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE-dt_per_dvy*Az)*dt_per_dvy*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R; // 3D
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE-dt_per_dvy*Az)*dt_per_dvy*R;

   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R; // 3D
   
   // ***** Case Ay < 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));

   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvy*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction wave:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE-dt_per_dvy*Az)*dt_per_dvy*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R; // 3D
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE-dt_per_dvy*Az)*dt_per_dvy*R;

   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R; // 3D
   
   // ***** Case Ay > 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvy*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE+dt_per_dvy*Az)*dt_per_dvy*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R; // 3D
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE+dt_per_dvy*Az)*dt_per_dvy*R;

   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R; // 3D
   
   // ***** Case Ay < 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvy*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE+dt_per_dvy*Az)*dt_per_dvy*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R; // 3D
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE+dt_per_dvy*Az)*dt_per_dvy*R;

   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R; // 3D

   // ***** Store calculated Flux *****
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] = myflux;
}

__global__ void cuda_acc_yfluxes3DcorrX(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   volatile Real Vx,Vy,Vz;
   Real Ax,Ay,Az,R,dt_per_dvy;
   Real myflux = 0.0f;
   cuint STATUS = 0;
   
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   Real xm2,xm1,xcc,xp1;
   
   // ***** 3D INCREMENT WAVES FROM VX-FACES *****
   dt_per_dvy = dt/blockParams[BlockParams::DVX];
   
   // 3D Increment wave, Ax > 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y-1,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   // 3D Increment wave, Ax > 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y-1,threadIdx.z+1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   // 3D Increment wave, Ax > 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y  ,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   // 3D Increment wave, Ax > 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y  ,threadIdx.z+1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE - dt_per_dvy*Ax)*R;
   
   // 3D Increment wave, Ax < 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;

   // ----- 3D Correction Waves -----
   
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y-1,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R; // ETUMERKKI ?
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R; // ETUMERKKI ?
   
   // 3D Increment wave, Ax < 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z+1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y-1,threadIdx.z+1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R;
   
   // 3D Increment wave, Ax < 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y  ,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R;
   
   // 3D Increment wave, Ax < 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y  ,threadIdx.z+1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y  ,threadIdx.z+1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVZ]*(ONE + Ax*dt_per_dvy)*R;
 
   // ***** Store Calculated Flux *****
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] += myflux;
}

__global__ void cuda_acc_yfluxes3DcorrZ(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   volatile Real Vx,Vy,Vz;
   Real Ax,Ay,Az,R,dt_per_dvy;
   Real myflux = 0.0f;
   cuint STATUS = 0;
   
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   Real xm2,xm1,xcc,xp1;
   
   // ***** 3D INCREMENT WAVES FROM VZ-FACES *****
   dt_per_dvy = dt/blockParams[BlockParams::DVZ];

   // 3D Increment wave, Ax > 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R;

   // 3D Increment wave, Ax < 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R; // ETUMERKKI?
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R; // ERUMERKKI?

   // 3D Increment wave, Ax > 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R;

   // 3D Increment wave, Ax < 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-2);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z+1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE - Az*dt_per_dvy)*R;

   // 3D Increment wave, Ax > 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;

   // 3D Increment wave, Ax < 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;

   // 3D Increment wave, Ax > 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z  );
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;

   // 3D Increment wave, Ax < 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+1)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z+1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z  );
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z+2);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;
   
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvy*dt/blockParams[BlockParams::DVX]*(ONE + Az*dt_per_dvy)*R;

   // *********************************
   // ***** Store calculated flux *****
   // *********************************
   //flux[offsetFlux+MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] = myflux;
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] += myflux;
}








__global__ void cuda_acc_zfluxes(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   Real Vx,Vy,Vz,Ax,Ay,Az,R,myflux,dt_per_dvz;
   cuint STATUS = tex1Dfetch(texRef_nbrsVel,offsetNbrs+MYBLOCK*28+27);
   
   // Fetch block parameters into shared mem array:
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   // Fetch cell parameters to shared mem array:
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   // Fetch the 4-point stencil of volume averages for this thread:
   Real xm2,xm1,xcc,xp1;
   
   // ******************************************************************
   // ***** Calculate flux at interface between (i,j,k), (i,j,k-1) *****
   // ******************************************************************
      
   // ***** Add contributions from vz-faces *****
   dt_per_dvz = dt/blockParams[BlockParams::DVZ];
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + threadIdx.z*blockParams[BlockParams::DVZ];
   
   // Contribution from Az > 0:
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y,threadIdx.z-1);
   myflux = Az*xm1;

   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y,threadIdx.z  );
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y,threadIdx.z-2);
   R = (xcc - xm1)*limiter((xm1-xm2)/(xcc-xm1+EPSILON));
   R = HALF*Az*(ONE - dt_per_dvz*Az)*R;
   myflux += R;

   // Contribution from Az < 0:
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += Az*xcc;

   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y,threadIdx.z+1);
   R = (xcc - xm1)*limiter((xp1-xcc)/(xcc-xm1+EPSILON));
   R = -HALF*Az*(ONE + dt_per_dvz*Az)*R;
   myflux += R;

   // ***** Add contributions from vx-faces *****
   dt_per_dvz = dt/blockParams[BlockParams::DVX];
   
   // ***** Case Ax > 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvz*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // Correction wave:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ax*Az*(ONE-dt_per_dvz*Ax)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE-Ax*dt_per_dvz)*R; // 3D
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Az*(ONE - dt_per_dvz*Ax)*dt_per_dvz*R;

   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R; // 3D
   
   // ***** Case Ax < 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));

   // Increment wave 2D+3D:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvz*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // Correction wave 2D + 3D:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ax*Az*(ONE + dt_per_dvz*Ax)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // 3D
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Az*(ONE+dt_per_dvz*Ax)*dt_per_dvz*R;

   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // 3D
   
   // ***** Case Ax > 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvz*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // 2D Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y,threadIdx.z  );
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ax*Az*(ONE - dt_per_dvz*Ax)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R; // 3D
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y,threadIdx.z  );
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(EPSILON+R));
   myflux -= HALF*Ax*Az*(ONE-dt_per_dvz*Ax)*dt_per_dvz*R;

   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R; // 3D
   
   // ***** Case Ax < 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ax*Az*dt_per_dvz*(SIXTH*Ay*dt/blockParams[BlockParams::DVY] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y,threadIdx.z  );
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ax*Az*(ONE + dt_per_dvz*Ax)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // 3D
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y,threadIdx.z  );
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ax*Az*(ONE+dt_per_dvz*Ax)*dt_per_dvz*R;

   Ay = fabs(accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // 3D
   
   // ***** Add contributions from vy-faces *****
   dt_per_dvz = dt/blockParams[BlockParams::DVY];
   
   // ***** Case Ay > 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvz*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-2,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE-dt_per_dvz*Ay)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - Ay*dt_per_dvz)*R; // 3D
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE-dt_per_dvz*Ay)*dt_per_dvz*R;

   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - Ay*dt_per_dvz)*R; // 3D
   
   // ***** Case Ay < 0, Az > 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvz*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+2,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE + dt_per_dvz*Ay)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + Ay*dt_per_dvz)*R; // 3D
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE + dt_per_dvz*Ay)*dt_per_dvz*R;

   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + Ay*dt_per_dvz)*R; // 3D
   
   // ***** Case Ay > 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvz*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction waves:
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-2,threadIdx.z  );
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE-dt_per_dvz*Ay)*dt_per_dvz*R;

   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - Ay*dt_per_dvz)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+1,threadIdx.z  );
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE-dt_per_dvz*Ay)*dt_per_dvz*R;
   
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - Ay*dt_per_dvz)*R;

   // ***** Case Ay < 0, Az < 0 *****
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   
   // Increment wave:
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z  );
   R = xcc-xm1;
   myflux += Ay*Az*dt_per_dvz*(SIXTH*Ax*dt/blockParams[BlockParams::DVX] - HALF)*R;

   // 2D Correction waves:
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+2,threadIdx.z  );
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += HALF*Ay*Az*(ONE + dt_per_dvz*Ay)*dt_per_dvz*R;
   
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + Ay*dt_per_dvz)*R; // 3D
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x,threadIdx.y-1,threadIdx.z  );
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= HALF*Ay*Az*(ONE + dt_per_dvz*Ay)*dt_per_dvz*R;
   
   Ax = fabs(accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + Ay*dt_per_dvz)*R; // 3D
   
   // ***** Store Calculated Flux *****
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] = myflux;
}

__global__ void cuda_acc_zfluxes3DcorrY(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   Real Vx,Vy,Vz,Ax,Ay,Az,R,myflux,dt_per_dvz;
   cuint STATUS = tex1Dfetch(texRef_nbrsVel,offsetNbrs+MYBLOCK*28+27);
   myflux = 0.0f;
   
   // Fetch block parameters into shared mem array:
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   // Fetch cell parameters to shared mem array:
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   // Fetch the 4-point stencil of volume averages for this thread:
   Real xm2,xm1,xcc,xp1;
   
   // **** 3D INCREMENT WAVES FROM VY-FACES *****
   dt_per_dvz = dt/blockParams[BlockParams::DVY];
   
   // 3D Increment wave, Ax > 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-2,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   // 3D Increment wave, Ax < 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;

   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-2,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   // 3D Increment wave, Ax > 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+2,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   // 3D Increment wave, Ax < 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z-1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+2,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   // 3D Increment wave, Ax > 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-2,threadIdx.z  );
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z  );
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   // 3D Increment wave, Ax < 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-2,threadIdx.z  );
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z  );
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE - dt_per_dvz*Ay)*R;
   
   // 3D Increment wave:  Ax > 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x-HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y  ,threadIdx.z  );
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+2,threadIdx.z  );
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   // 3D Increment wave, Ax < 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1+HALF)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y  ,threadIdx.z  );
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+2,threadIdx.z  );
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;
   
   Vy = blockParams[BlockParams::VYCRD] + threadIdx.y*blockParams[BlockParams::DVY];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z  );
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVX]*(ONE + dt_per_dvz*Ay)*R;

   // *********************************
   // ***** STORE CALCULATED FLUX *****
   // *********************************
   //flux[offsetFlux+MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] = myflux;
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] += myflux;
}

__global__ void cuda_acc_zfluxes3DcorrX(uint offsetAvgs,uint offsetFlux,uint offsetBlockParams,uint offsetCellParams,uint offsetNbrs,Real* flux,Real dt,Real q_per_m) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   
   int i;
   Real Vx,Vy,Vz,Ax,Ay,Az,R,myflux,dt_per_dvz;
   cuint STATUS = tex1Dfetch(texRef_nbrsVel,offsetNbrs+MYBLOCK*28+27);
   myflux = 0.0f;
   
   // Fetch block parameters into shared mem array:
   __shared__ Real blockParams[SIZE_BLOCKPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_BLOCKPARAMS;
   blockParams[i] = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+i);
   
   // Fetch cell parameters to shared mem array:
   __shared__ Real cellParams[SIZE_CELLPARAMS];
   i = index(threadIdx.x,threadIdx.y,threadIdx.z) % SIZE_CELLPARAMS;
   cellParams[i] = tex1Dfetch(texRef_cellParams,offsetCellParams+i);
   
   // Fetch the 4-point stencil of volume averages for this thread:
   Real xm2,xm1,xcc,xp1;
   
   // **** 3D INCREMENT WAVES FROM VX-FACES *****
   dt_per_dvz = dt/blockParams[BlockParams::DVX];
   
   // 3D Increment wave, Ax > 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // --- 3D Correction waves ---
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y-1,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*dt_per_dvz*dt/blockParams[BlockParams::DVY]*Ax*Ay*Az*(ONE - Ax*dt_per_dvz)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*dt_per_dvz*dt/blockParams[BlockParams::DVY]*Ax*Ay*Az*(ONE - Ax*dt_per_dvz)*R;
   
   // 3D Increment wave, Ax > 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y+1,threadIdx.z-1);
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*dt_per_dvz*dt/blockParams[BlockParams::DVY]*Ax*Ay*Az*(ONE - Ax*dt_per_dvz)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z-1);
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*dt_per_dvz*dt/blockParams[BlockParams::DVY]*Ax*Ay*Az*(ONE - Ax*dt_per_dvz)*R;
   
   // 3D Increment wave, Ax < 0, Ay > 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z-1);
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y-1,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // ETUMERRKI
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // ETUMERKKI
   
   // 3D Increment wave, Ax < 0, Ay < 0, Az > 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z-HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z-1);
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z-1);
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Wave -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y+1,threadIdx.z-1);
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // ETUMERKKI?
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmax(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z-1);
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R; // ETUMERKKI?
   
   // 3D Increment wave, Ax > 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y-1,threadIdx.z  );
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z  );
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R;
   
   // 3D Increment wave, Ax > 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z  );
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Waves -----
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-2,threadIdx.y+1,threadIdx.z  );
   R *= limiter((xm1-xm2)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Ax = fmax(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z  );
   R = xp1-xcc;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE - Ax*dt_per_dvz)*R;
   
   // 3D Increment wave, Ax < 0, Ay > 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y-HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y-1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y-1,threadIdx.z  );
   R = xcc-xm1;
   myflux -= SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Waves -----
   xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y-1,threadIdx.z  );
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmax(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y-1,threadIdx.z  );
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R;
   
   // 3D Increment wave, Ax < 0, Ay < 0, Az < 0:
   Vx = blockParams[BlockParams::VXCRD] + (threadIdx.x+1)*blockParams[BlockParams::DVX];
   Vy = blockParams[BlockParams::VYCRD] + (threadIdx.y+1+HALF)*blockParams[BlockParams::DVY];
   Vz = blockParams[BlockParams::VZCRD] + (threadIdx.z+HALF)*blockParams[BlockParams::DVZ];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xcc = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y+1,threadIdx.z  );
   xm1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y+1,threadIdx.z  );
   R = xcc-xm1;
   myflux += SIXTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*R;
   
   // ----- 3D Correction Waves -----
    xp1 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x+2,threadIdx.y+1,threadIdx.z  );
   R *= limiter((xp1-xcc)/(R+EPSILON));
   myflux -= FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R;
   
   Vx = blockParams[BlockParams::VXCRD] + threadIdx.x*blockParams[BlockParams::DVX];
   Ax = fmin(ZERO,accelerationX(Vx,Vy,Vz,cellParams,q_per_m));
   Ay = fmin(ZERO,accelerationY(Vx,Vy,Vz,cellParams,q_per_m));
   Az = fmin(ZERO,accelerationZ(Vx,Vy,Vz,cellParams,q_per_m));
   xm2 = fetchAvgs(MYBLOCK,offsetAvgs,offsetNbrs,STATUS,threadIdx.x-1,threadIdx.y+1,threadIdx.z  );
   R = xm1-xm2;
   R *= limiter((xcc-xm1)/(R+EPSILON));
   myflux += FOURTH*Ax*Ay*Az*dt_per_dvz*dt/blockParams[BlockParams::DVY]*(ONE + Ax*dt_per_dvz)*R;
   
   // ***** STORE CALCULATED FLUX *****
   flux[MYBLOCK*WID3 + index(threadIdx.x,threadIdx.y,threadIdx.z)] += myflux;
}

__global__ void cuda_acc_vx_changes(Real* dfdt,uint offsetDeriv,uint offsetFlux,uint offsetNbrs,uint offsetBlockParams,Real dt) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   cuint MYIND = index(threadIdx.x,threadIdx.y,threadIdx.z);

   cuint STATUS = tex1Dfetch(texRef_nbrsVel    ,offsetNbrs+MYBLOCK*28+27);   
   creal DVX    = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+BlockParams::DVX);
   
   creal Fx_neg = fetchFlux(MYBLOCK,offsetFlux,offsetNbrs,STATUS,threadIdx.x  ,threadIdx.y,threadIdx.z);
   creal Fx_pos = fetchFlux(MYBLOCK,offsetFlux,offsetNbrs,STATUS,threadIdx.x+1,threadIdx.y,threadIdx.z);
   
   dfdt[MYBLOCK*WID3 + MYIND] = dt/DVX*(Fx_neg-Fx_pos);
}

__global__ void cuda_acc_vy_changes(Real* dfdt,uint offsetDeriv,uint offsetFlux,uint offsetNbrs,uint offsetBlockParams,Real dt) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   cuint MYIND = index(threadIdx.x,threadIdx.y,threadIdx.z);
   
   cuint STATUS = tex1Dfetch(texRef_nbrsVel    ,offsetNbrs+MYBLOCK*28+27);
   creal DVY    = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+BlockParams::DVY);
   
   creal Fy_neg = fetchFlux(MYBLOCK,offsetFlux,offsetNbrs,STATUS,threadIdx.x,threadIdx.y  ,threadIdx.z);
   creal Fy_pos = fetchFlux(MYBLOCK,offsetFlux,offsetNbrs,STATUS,threadIdx.x,threadIdx.y+1,threadIdx.z);
   
   dfdt[MYBLOCK*WID3 + MYIND] += dt/DVY*(Fy_neg-Fy_pos);
}

__global__ void cuda_acc_vz_changes(Real* dfdt,uint offsetDeriv,uint offsetFlux,uint offsetNbrs,uint offsetBlockParams,Real dt) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   cuint MYIND = index(threadIdx.x,threadIdx.y,threadIdx.z);
   
   cuint STATUS = tex1Dfetch(texRef_nbrsVel    ,offsetNbrs+MYBLOCK*28+27);
   creal DVZ    = tex1Dfetch(texRef_blockParams,offsetBlockParams+MYBLOCK*SIZE_BLOCKPARAMS+BlockParams::DVZ);
   
   creal Fz_neg = fetchFlux(MYBLOCK,offsetFlux,offsetNbrs,STATUS,threadIdx.x,threadIdx.y,threadIdx.z  );
   creal Fz_pos = fetchFlux(MYBLOCK,offsetFlux,offsetNbrs,STATUS,threadIdx.x,threadIdx.y,threadIdx.z+1);
   
   dfdt[MYBLOCK*WID3 + MYIND] += dt/DVZ*(Fz_neg-Fz_pos);
}

__global__ void cuda_acc_propagate(Real* avgs,Real* dfdt,uint offsetAvgs,uint offsetDeriv) {
   cuint MYBLOCK = blockIndex(blockIdx.x,blockIdx.y);
   cuint MYIND = index(threadIdx.x,threadIdx.y,threadIdx.z);
   avgs[MYBLOCK*WID3 + MYIND] += dfdt[MYBLOCK*WID3 + MYIND];
}
