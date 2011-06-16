#ifndef CPU_TRANS_H
#define CPU_TRANS_H

#include <vector>
#include "../definitions.h"
#include "../common.h"
#include "cpu_common.h"

const uint I0_J0 = 0;
const uint I1_J0 = 1;
const uint I2_J0 = 2;
const uint I0_J1 = 3;
const uint I1_J1 = 4;
const uint I2_J1 = 5;
const uint I0_J2 = 6;
const uint I1_J2 = 7;
const uint I2_J2 = 8;

const Real EPSILON = 1.0e-25;
const Real ONE     = 1.0;
const Real ZERO    = 0.0;
const Real HALF    = 0.5;

#ifdef SPAT3D
   #warning IJ indices not implemented for 3D
#endif

template<typename REAL> REAL limiter(const REAL& theta_up,const REAL& theta_lo,const REAL& xcc) {
   return superbee(theta_up/theta_lo);
}

template<typename REAL> bool signs(const REAL& a,const REAL& b,REAL& sign) {
   const REAL PLUS_ONE = 1.0;
   const REAL MINUS_ONE = -1.0;
   if (a >= ZERO) {
      if (b >= ZERO) {
	 sign = PLUS_ONE;
	 return true;
      } else {
	 return false;
      }
   } else {
      if (b < ZERO) {
	 sign = MINUS_ONE;
	 return true;
      } else {
	 return false;
      }
   }
}

template<typename REAL> void cpu_blockVelocityMoments(const REAL* const avgs,const REAL* const blockParams,REAL* const cellParams) {
   const REAL HALF = 0.5;
   
   REAL n_sum = 0.0;
   REAL nvx_sum = 0.0;
   REAL nvy_sum = 0.0;
   REAL nvz_sum = 0.0;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      const REAL VX = blockParams[BlockParams::VXCRD] + HALF*blockParams[BlockParams::DVX];
      const REAL VY = blockParams[BlockParams::VYCRD] + HALF*blockParams[BlockParams::DVY];
      const REAL VZ = blockParams[BlockParams::VZCRD] + HALF*blockParams[BlockParams::DVZ];
      
      n_sum   += avgs[cellIndex(i,j,k)];
      nvx_sum += avgs[cellIndex(i,j,k)]*VX;
      nvy_sum += avgs[cellIndex(i,j,k)]*VY;
      nvz_sum += avgs[cellIndex(i,j,k)]*VZ;
   }
   
   // Accumulate contributions coming from this velocity block to the 
   // spatial cell velocity moments. If multithreading / OpenMP is used, 
   // these updates need to be atomic:
   const REAL DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   cellParams[CellParams::RHO  ] += n_sum * DV3;
   cellParams[CellParams::RHOVX] += nvx_sum * DV3;
   cellParams[CellParams::RHOVY] += nvy_sum * DV3;
   cellParams[CellParams::RHOVZ] += nvz_sum * DV3;
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcSpatDerivs(CELL& cell,const UINT& BLOCK,const std::vector<const CELL*>& nbrPtrs) { }

/**
 * @param AVGS Pointer to array containing volume averages.
 * @param CELL_PARAMS Pointer to array containing spatial cell parameters.
 * @param BLOCK_PARAMS Pointer to array containing velocity block parameters.
 * @param flux Array in which calculated df/dt values are copied.
 * @param nbrsSpa Array containing spatial space neighbour lists.
 * @param BLOCK Which velocity block is to be calculated, acts as an offset into spatial neighbour list.
 * @param DT Time step.
 */
template<typename REAL,typename UINT> void cpu_calcSpatDfdt(const REAL* const AVGS,const REAL* const CELL_PARAMS,const REAL* const BLOCK_PARAMS,REAL* const flux,
							    const UINT* const nbrsSpa,const UINT& BLOCK,const REAL& DT) {
   // NOTES ON OPTIMIZATION:
   // 
   // Skipping for-loops based on the signs of velocity components reduced 
   // time consumption by 37.6% on Arto's laptop and by 46.3% on meteo, 
   // compared to a version where all for loops are calculated.
   // Practically all performance metrics (flops, cache usage, etc.) got 
   // slightly worse, but the overall time spent in function was cut by a 
   // large margin.
   
   // Create a temporary buffer for storing df/dt updates and init to zero value:
   const UINT SIZE_FLUXBUFFER = 9*WID3;
   REAL dfdt[SIZE_FLUXBUFFER];
   for (uint i=0; i<SIZE_FLUXBUFFER; ++i) dfdt[i] = 0.0;
   #ifdef SPAT3D
      #warning Not enough dfdt buffers for 3D
   #endif
   
   // Pointer to velocity block whose df/dt contributions are calculated:
   const REAL* const blockAvgs   = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 13]*SIZE_VELBLOCK;
   const REAL* const blockParams = BLOCK_PARAMS + BLOCK*SIZE_BLOCKPARAMS;
   
   // ***** Consider the interface between (i-1,j,k) and (i,j,k): *****
   const REAL dt_per_dx = DT / CELL_PARAMS[CellParams::DX];
   const REAL dt_per_dy = DT / CELL_PARAMS[CellParams::DY];
   const REAL* nbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 12]*SIZE_VELBLOCK; //  -x nbr
   
   // Case Vx > 0:
   Real vx_sign,vy_sign;
   bool sign_vx_constant = signs(blockParams[BlockParams::VXCRD],blockParams[BlockParams::VXCRD] + (WID-1+HALF)*blockParams[BlockParams::DVX],vx_sign);
   bool sign_vy_constant = signs(blockParams[BlockParams::VYCRD],blockParams[BlockParams::VYCRD] + (WID-1+HALF)*blockParams[BlockParams::DVY],vy_sign);
   
   if (sign_vx_constant == false || vx_sign > ZERO) {
      const REAL* const nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 27]*SIZE_VELBLOCK; // --x nbr
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 const UINT MYIND = cellIndex(i,j,k);
	 const REAL Vx = max(ZERO,blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX]);
	 const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];

	 const REAL xcc = blockAvgs[MYIND];
	 const REAL xm1 = nbr_minus1[MYIND];
	 REAL R = xcc - xm1;

	 // Increment waves:
	 dfdt[I1_J1*WID3 + MYIND] += Vx*xm1*dt_per_dx;
	 dfdt[I0_J1*WID3 + MYIND] -= Vx*xm1*dt_per_dx;

	 // Transverse increment waves:
	 dfdt[I1_J2*WID3 + MYIND] -= HALF*dt_per_dx*Vx*max(ZERO,Vy)*R * dt_per_dy; // Vy > 0
	 dfdt[I1_J1*WID3 + MYIND] += HALF*dt_per_dx*Vx*max(ZERO,Vy)*R * dt_per_dy;
	 dfdt[I1_J0*WID3 + MYIND] += HALF*dt_per_dx*Vx*min(ZERO,Vy)*R * dt_per_dy; // Vy < 0
	 dfdt[I1_J1*WID3 + MYIND] -= HALF*dt_per_dx*Vx*min(ZERO,Vy)*R * dt_per_dy;
	 
	 // Correction waves:
	 const REAL xm2 = nbr_minus2[MYIND];
	 R *= limiter(xm1-xm2,R+EPSILON,xcc);
	 const REAL corr_wave = HALF*Vx*(ONE - dt_per_dx*Vx)*R;
	 dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dx;
	 dfdt[I0_J1*WID3 + MYIND] -= corr_wave * dt_per_dx;

	 // Transverse correction waves:
	 dfdt[I1_J2*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy > 0
	 dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J2*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J1*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;

	 dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy < 0
	 dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
      }
   }
   
   // Case Vx < 0:
   if (sign_vx_constant == false || vx_sign < ZERO) {
      const REAL* const nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 14]*SIZE_VELBLOCK; //  +x nbr
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 const UINT MYIND = cellIndex(i,j,k);
	 const REAL Vx = min(ZERO,blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX]);
	 const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	 
	 const REAL xcc = blockAvgs[cellIndex(i,j,k)];
	 const REAL xm1 = nbr_minus1[cellIndex(i,j,k)];
	 REAL R = xcc - xm1;
	 
	 // Increment waves:
	 dfdt[I1_J1*WID3 + MYIND] += Vx*xcc*dt_per_dx;
	 dfdt[I0_J1*WID3 + MYIND] -= Vx*xcc*dt_per_dx;
	 
	 // Transverse increment waves:
	 dfdt[I0_J2*WID3 + MYIND] -= HALF*dt_per_dx*Vx*fmax(ZERO,Vy)*R * dt_per_dy; // Vy > 0
	 dfdt[I0_J1*WID3 + MYIND] += HALF*dt_per_dx*Vx*fmax(ZERO,Vy)*R * dt_per_dy;
	 
	 dfdt[I0_J0*WID3 + MYIND] += HALF*dt_per_dx*Vx*fmin(ZERO,Vy)*R * dt_per_dy; // Vy < 0
	 dfdt[I0_J1*WID3 + MYIND] -= HALF*dt_per_dx*Vx*fmin(ZERO,Vy)*R * dt_per_dy;

	 // Correction waves:
	 const REAL xm2 = nbr_minus2[MYIND];
	 R *= limiter(xm2-xcc,R+EPSILON,xcc);
	 const REAL corr_wave = -HALF*Vx*(ONE + dt_per_dx*Vx)*R;
	 dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dx;
	 dfdt[I0_J1*WID3 + MYIND] -= corr_wave * dt_per_dx;

	 // Transverse correction waves:
	 dfdt[I1_J2*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy > 0
	 dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J2*WID3 + MYIND] -= fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J1*WID3 + MYIND] += fmax(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 
	 dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy; // Vy < 0
	 dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
	 dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vy)*dt_per_dx*corr_wave * dt_per_dy;
      }
   }
   
   // ***** Consider the interface between (i,j-1,k) and (i,j,k): *****
   nbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 10]*SIZE_VELBLOCK;
   
   // Case Vy > 0:
   if (sign_vy_constant == false || vy_sign > ZERO) {
      const REAL* const nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 28]*SIZE_VELBLOCK; // --y nbr
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 const UINT MYIND = cellIndex(i,j,k);
	 const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	 const REAL Vy = max(ZERO,blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY]);
	 
	 const REAL xcc = blockAvgs[cellIndex(i,j,k)];
	 const REAL xm1 = nbr_minus1[cellIndex(i,j,k)];
	 REAL R = xcc - xm1;
      
	 // Increment waves:
	 dfdt[I1_J1*WID3 + MYIND] += Vy*xm1*dt_per_dy;
	 dfdt[I1_J0*WID3 + MYIND] -= Vy*xm1*dt_per_dy;
	 
	 // Transverse increment waves:
	 dfdt[I2_J1*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx; // Vx > 0
	 dfdt[I1_J1*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx;
	 
	 dfdt[I0_J1*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx; // Vx < 0
	 dfdt[I1_J1*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx;
	 
	 // Correction waves:
	 const REAL xm2 = nbr_minus2[MYIND];
	 R *= limiter(xm1-xm2,R+EPSILON,xcc);
	 const REAL corr_wave = HALF*Vy*(ONE - dt_per_dy*Vy)*R;
	 dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dy;
	 dfdt[I1_J0*WID3 + MYIND] -= corr_wave * dt_per_dy;
	 
	 // Transverse correction waves:
	 dfdt[I2_J1*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx > 0
	 dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I2_J0*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I1_J0*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 
	 dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx < 0
	 dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
      }
   }
   
   // Case Vy < 0:
   if (sign_vy_constant == false || vy_sign < ZERO) {
      const REAL* const nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 16]*SIZE_VELBLOCK; //  +y nbr
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 const UINT MYIND = cellIndex(i,j,k);
	 const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	 const REAL Vy = min(ZERO,blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY]);
	 
	 const REAL xcc = blockAvgs[cellIndex(i,j,k)];
	 const REAL xm1 = nbr_minus1[cellIndex(i,j,k)];
	 REAL R = xcc - xm1;
	 
	 // Increment waves:
	 dfdt[I1_J1*WID3 + MYIND] += Vy*xcc*dt_per_dy;
	 dfdt[I1_J0*WID3 + MYIND] -= Vy*xcc*dt_per_dy;
	 
	 // Transverse increment waves:
	 dfdt[I2_J0*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx; // Vx > 0
	 dfdt[I1_J0*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmax(ZERO,Vx)*R * dt_per_dx;
	 
	 dfdt[I0_J0*WID3 + MYIND] += HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx; // Vx < 0
	 dfdt[I1_J0*WID3 + MYIND] -= HALF*dt_per_dy*Vy*fmin(ZERO,Vx)*R * dt_per_dx;
	 
	 // Correction waves:
	 const REAL xm2 = nbr_minus2[MYIND];
	 R *= limiter(xm2-xcc,R+EPSILON,xcc);
	 const REAL corr_wave = -HALF*Vy*(ONE + dt_per_dy*Vy)*R;
	 dfdt[I1_J1*WID3 + MYIND] += corr_wave * dt_per_dy;
	 dfdt[I1_J0*WID3 + MYIND] -= corr_wave * dt_per_dy;
	 
	 // Transverse correction waves:
	 dfdt[I2_J1*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx > 0
	 dfdt[I1_J1*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I2_J0*WID3 + MYIND] -= fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I1_J0*WID3 + MYIND] += fmax(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;

	 dfdt[I0_J1*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx; // Vx < 0
	 dfdt[I1_J1*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I0_J0*WID3 + MYIND] += fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
	 dfdt[I1_J0*WID3 + MYIND] -= fmin(ZERO,Vx)*dt_per_dy*corr_wave * dt_per_dx;
      }
   }

   // Accumulate calculated df/dt values from temporary buffer to 
   // main memory. If multithreading is used, these updates need 
   // to be atomistic:
   const UINT boundaryFlags = nbrsSpa[BLOCK*SIZE_NBRS_SPA + 30];
   for (uint nbr=0; nbr<9; ++nbr) {
      // If the neighbour does not exist, do not copy data:
      if (((boundaryFlags >> (9+nbr)) & 1) == 0) continue;
      
      const UINT nbrBlock = nbrsSpa[BLOCK*SIZE_NBRS_SPA + 9 + nbr];
      for (uint i=0; i<SIZE_VELBLOCK; ++i) flux[nbrBlock*WID3 + i] += dfdt[nbr*WID3 + i];
   }

   #ifdef SPAT3D
      #warning dfdt copies not implemented for 3D in cpu_calcSpatDfdt
   #endif
}

template<typename REAL,typename UINT> void cpu_propagateSpat(REAL* const avgs,const REAL* const flux,const REAL* const nbrFluxes,
							     const REAL* const blockParams,const REAL* const cellParams,const UINT& BLOCK) {
   // Propagate distribution function:
   if (nbrFluxes == NULL) {
      // No remote neighbour contributions to df/dt 
      for (UINT i=0; i<WID3; ++i) {
	 avgs[BLOCK*WID3 + i] += flux[BLOCK*WID3 + i];
      }
   } else {
      // Cell has remote neighbour contributions to df/dt
      for (UINT i=0; i<WID3; ++i) {
	 avgs[BLOCK*WID3 + i] += flux[BLOCK*WID3 + i] + nbrFluxes[BLOCK*WID3 + i];
      }
   }
}

/** Propagate distribution function in spatial space and calculate 
 * velocity moments rho, rhovx, rhovy, and rhovz simultaneously.
 * @param avgs Array containing volume averages of distribution function for this 
 * spatial cell.
 * @param flux Array containing the local contributions to df/dt.
 * @param nbrFlux Array containing remote neighbour contributions to df/dt.
 * @param blockParams Array containing velocity block parameters for this spatial cell.
 * @param cellParams Array containing spatial cell parameters.
 * @param BLOCK Which velocity block is being propagated, acts as an index into arrays.
 */
template<typename REAL,typename UINT> void cpu_propagateSpatWithMoments(REAL* const avgs,const REAL* const flux,const REAL* const nbrFluxes,
								       const REAL* const blockParams,REAL* const cellParams,const UINT& BLOCK) {
   // Propagate distribution function:
   if (nbrFluxes == NULL) {
      // No remote neighbour contributions to df/dt
      for (UINT i=0; i<WID3; ++i) {
	 avgs[BLOCK*WID3 + i] += flux[BLOCK*WID3 + i];
      }
   } else {
      // Cell has remote neighbour contributions to df/dt
      for (UINT i=0; i<WID3; ++i) {
	 avgs[BLOCK*WID3 + i] += flux[BLOCK*WID3 + i] + nbrFluxes[BLOCK*WID3 + i];
      }
   }
   
   // Calculate velocity moments:
   cpu_blockVelocityMoments(avgs+BLOCK*SIZE_VELBLOCK,blockParams+BLOCK*SIZE_BLOCKPARAMS,cellParams);
}

template<typename REAL,typename UINT> void cpu_calcVelocityMoments(const REAL* const avgs,const REAL* const blockParams,REAL* const cellParams,const UINT& BLOCK) {
   // Calculate velocity moments:
   cpu_blockVelocityMoments(avgs+BLOCK*SIZE_VELBLOCK,blockParams+BLOCK*SIZE_BLOCKPARAMS,cellParams);
}

#endif
