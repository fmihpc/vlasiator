#ifndef CPU_TRANS_H
#define CPU_TRANS_H

#include <vector>
#include "../definitions.h"
#include "../common.h"
#include "cpu_common.h"

const uint I0_J0_K0 = 0*WID3;
const uint I1_J0_K0 = 1*WID3;
const uint I2_J0_K0 = 2*WID3;
const uint I0_J1_K0 = 3*WID3;
const uint I1_J1_K0 = 4*WID3;
const uint I2_J1_K0 = 5*WID3;
const uint I0_J2_K0 = 6*WID3;
const uint I1_J2_K0 = 7*WID3;
const uint I2_J2_K0 = 8*WID3;

const uint I0_J0_K1 = 9*WID3;
const uint I1_J0_K1 = 10*WID3;
const uint I2_J0_K1 = 11*WID3;
const uint I0_J1_K1 = 12*WID3;
const uint I1_J1_K1 = 13*WID3;
const uint I2_J1_K1 = 14*WID3;
const uint I0_J2_K1 = 15*WID3;
const uint I1_J2_K1 = 16*WID3;
const uint I2_J2_K1 = 17*WID3;

const uint I0_J0_K2 = 18*WID3;
const uint I1_J0_K2 = 19*WID3;
const uint I2_J0_K2 = 20*WID3;
const uint I0_J1_K2 = 21*WID3;
const uint I1_J1_K2 = 22*WID3;
const uint I2_J1_K2 = 23*WID3;
const uint I0_J2_K2 = 24*WID3;
const uint I1_J2_K2 = 25*WID3;
const uint I2_J2_K2 = 26*WID3;

const Real EPSILON = 1.0e-25;
const Real ZERO    = 0.0;
const Real SIXTH   = 1.0/6.0;
const Real FOURTH  = 0.25;
const Real HALF    = 0.5;
const Real ONE     = 1.0;
creal TWO    = 2.0;

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
   
   // Create a temporary buffer for storing df/dt updates and init to zero value:
   //const UINT SIZE_FLUXBUFFER = 9*WID3;
   const UINT SIZE_FLUXBUFFER = 27*WID3;
   REAL dfdt[SIZE_FLUXBUFFER];
   for (uint i=0; i<SIZE_FLUXBUFFER; ++i) dfdt[i] = 0.0;
   #ifdef SPAT3D
      #warning Not enough dfdt buffers for 3D
   #endif
   
   // Pointer to velocity block whose df/dt contributions are calculated:
   const REAL* const blockAvgs   = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 13]*SIZE_VELBLOCK;
   const REAL* const blockParams = BLOCK_PARAMS + BLOCK*SIZE_BLOCKPARAMS;
   
   Real vx_sign,vy_sign,vz_sign;
   bool sign_vx_constant = signs(blockParams[BlockParams::VXCRD],blockParams[BlockParams::VXCRD] + (WID-1+HALF)*blockParams[BlockParams::DVX],vx_sign);
   bool sign_vy_constant = signs(blockParams[BlockParams::VYCRD],blockParams[BlockParams::VYCRD] + (WID-1+HALF)*blockParams[BlockParams::DVY],vy_sign);
   bool sign_vz_constant = signs(blockParams[BlockParams::VZCRD],blockParams[BlockParams::VZCRD] + (WID-1+HALF)*blockParams[BlockParams::DVZ],vz_sign);

   // ***** Consider the interface between (i-1,j,k) and (i,j,k): *****
   const REAL dt_per_dx = DT / CELL_PARAMS[CellParams::DX];
   const REAL dt_per_dy = DT / CELL_PARAMS[CellParams::DY];
   const REAL dt_per_dz = DT / CELL_PARAMS[CellParams::DZ];

   const REAL* const xnbr_plus1  = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 14]*SIZE_VELBLOCK; //  +x nbr
   const REAL* const xnbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 12]*SIZE_VELBLOCK; //  -x nbr
   const REAL* const xnbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 27]*SIZE_VELBLOCK; // --x nbr
   const REAL* const ynbr_plus1  = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 16]*SIZE_VELBLOCK; //  +y nbr
   const REAL* const ynbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 10]*SIZE_VELBLOCK;
   const REAL* const ynbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 28]*SIZE_VELBLOCK; // --y nbr
   const REAL* const znbr_plus1  = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 22]*SIZE_VELBLOCK; // +z nbr
   const REAL* const znbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 4 ]*SIZE_VELBLOCK;
   const REAL* const znbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 29]*SIZE_VELBLOCK; // --z nbr
   
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];

      // Fetch required values from global arrays:
      const UINT MYIND = cellIndex(i,j,k);
      const REAL xcc = blockAvgs[MYIND];
      const REAL xp1 = xnbr_plus1[MYIND];
      const REAL xm1 = xnbr_minus1[MYIND];
      const REAL xm2 = xnbr_minus2[MYIND];
      const REAL yp1 = ynbr_plus1[MYIND];
      const REAL ym1 = ynbr_minus1[MYIND];
      const REAL ym2 = ynbr_minus2[MYIND];
      const REAL zp1 = znbr_plus1[MYIND];
      const REAL zm1 = znbr_minus1[MYIND];
      const REAL zm2 = znbr_minus2[MYIND];
      
      if (Vx >= ZERO) {
	 // 1D increment and correction waves:
	 const REAL RX = xcc - xm1;
	 const REAL incrWaveX = Vx*xm1*dt_per_dx;
	 dfdt[I1_J1_K1 + MYIND] += incrWaveX;
	 dfdt[I0_J1_K1 + MYIND] -= incrWaveX;
	 #ifndef LEV_1ST_ORDER
	 const REAL corrWaveX = HALF*Vx*(ONE - dt_per_dx*Vx)*RX*limiter(xm1-xm2,RX+EPSILON,xcc)*dt_per_dx;
	 dfdt[I1_J1_K1 + MYIND] += corrWaveX;
	 dfdt[I0_J1_K1 + MYIND] -= corrWaveX;
         #endif
	 
	 const REAL RY = xcc - ym1;
	 const REAL RZ = xcc - zm1;
	 if (Vy >= ZERO) { // Vx > 0 Vy > 0
	    const REAL transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J2_K1 + MYIND] -= transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYIND] += transIncrWave*RX;
	    dfdt[I2_J1_K1 + MYIND] -= transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYIND] += transIncrWave*RY;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*RX*limiter(xm1-xm2,RX+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I1_J2_K1 + MYIND] += transCorrWaveXY;
	    dfdt[I0_J2_K1 + MYIND] -= transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYIND] += transCorrWaveXY;
	    const REAL transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*RY*limiter(ym1-ym2,RY+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I2_J1_K1 + MYIND] += transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveYX;
	    dfdt[I2_J0_K1 + MYIND] -= transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYIND] += transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx > 0 Vy > 0 Vz > 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J2_K2 + MYIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J1_K2 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I1_J2_K1 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I2_J1_K2 + MYIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K2 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J1_K1 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J2_K1 + MYIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J2_K1 + MYIND] += doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K1 + MYIND] += doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE - dt_per_dx*Vx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xm1-xm2,RX+EPSILON,xcc);
	       dfdt[I1_J2_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(ym1-ym2,RY+EPSILON,xcc);
	       dfdt[I2_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K2 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc);
	       dfdt[I2_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYIND] += doubleTransCorrWaveZXY;
	       #endif
	    } else { // Vx > 0 Vy > 0 Vz < 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J2_K0 + MYIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I1_J2_K1 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I2_J1_K0 + MYIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J1_K1 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J2_K0 + MYIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J2_K0 + MYIND] += doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K0 + MYIND] += doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE - Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xm1-xm2,RX+EPSILON,xcc);
	       dfdt[I1_J2_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(ym1-ym2,RY+EPSILON,xcc);
	       dfdt[I2_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc);
	       dfdt[I2_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYIND] += doubleTransCorrWaveZXY;
	       #endif
	    }
	 } else { // Vx > 0 Vy < 0
	    const REAL transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J0_K1 + MYIND] += transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYIND] -= transIncrWave*RX;
	    dfdt[I2_J0_K1 + MYIND] -= transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYIND] += transIncrWave*RY;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*RX*limiter(xm1-xm2,RX+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveXY;
	    dfdt[I0_J0_K1 + MYIND] += transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveXY;
	    const REAL transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*RY*limiter(yp1-xcc,RY+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I2_J1_K1 + MYIND] += transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveYX;
	    dfdt[I2_J0_K1 + MYIND] -= transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYIND] += transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx > 0 Vy < 0 Vz > 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J0_K2 + MYIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J1_K2 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I2_J0_K2 + MYIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J0_K2 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J0_K1 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J0_K1 + MYIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K1 + MYIND] -= doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xm1-xm2,RX+EPSILON,xcc);
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(yp1-xcc,RY+EPSILON,xcc);
	       dfdt[I2_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K2 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc);
	       dfdt[I2_J0_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J0_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       #endif
	    } else { // Vx > 0 Vy < 0 Vz < 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J0_K0 + MYIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J0_K1 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I2_J0_K0 + MYIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J0_K0 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J0_K1 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J0_K0 + MYIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J0_K0 + MYIND] -= doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K0 + MYIND] -= doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xm1-xm2,RX+EPSILON,xcc);
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(yp1-xcc,RY+EPSILON,xcc);
	       dfdt[I2_J1_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc);
	       dfdt[I2_J0_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J0_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       #endif
	    }
	 }
	 
	 if (Vz >= ZERO) {
	    const REAL transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYIND] -= transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYIND] += transIncrWave*RX;
	    dfdt[I2_J1_K1 + MYIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYIND] += transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*RX*limiter(xm1-xm2,RX+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I1_J1_K2 + MYIND] += transCorrWaveXZ;
	    dfdt[I0_J1_K2 + MYIND] -= transCorrWaveXZ;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYIND] += transCorrWaveXZ;
	    const REAL transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I2_J1_K1 + MYIND] += transCorrWaveZX;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveZX;
	    dfdt[I2_J1_K0 + MYIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYIND] += transCorrWaveZX;
	    #endif
	 } else {
	    const REAL transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K0 + MYIND] += transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYIND] -= transIncrWave*RX;
	    dfdt[I2_J1_K0 + MYIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYIND] += transIncrWave*RZ;

	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*RX*limiter(xm1-xm2,RX+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveXZ;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K0 + MYIND] += transCorrWaveXZ;
	    const REAL transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I2_J1_K1 + MYIND] += transCorrWaveZX;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveZX;
	    dfdt[I2_J1_K0 + MYIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYIND] += transCorrWaveZX;
            #endif
	 }
      } else { // Vx < 0
	 const REAL RX = xcc - xm1;
	 const REAL incrWaveX = Vx*xcc*dt_per_dx;
	 dfdt[I1_J1_K1 + MYIND] += incrWaveX;
	 dfdt[I0_J1_K1 + MYIND] -= incrWaveX;
	 #ifndef LEV_1ST_ORDER
	 const REAL corrWaveX = -HALF*Vx*(ONE + dt_per_dx*Vx)*RX*limiter(xp1-xcc,RX+EPSILON,xcc)*dt_per_dx;
	 dfdt[I1_J1_K1 + MYIND] += corrWaveX;
	 dfdt[I0_J1_K1 + MYIND] -= corrWaveX;
	 #endif
  
	 const REAL RY = xcc - ym1;
	 const REAL RZ = xcc - zm1;
	 if (Vy >= ZERO) { // Vx < 0 Vy > 0
	    const REAL transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J2_K1 + MYIND] -= transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYIND] += transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYIND] -= transIncrWave*RY;

	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*RX*limiter(xp1-xcc,RX+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I1_J2_K1 + MYIND] += transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveXY;
	    dfdt[I0_J2_K1 + MYIND] -= transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYIND] += transCorrWaveXY;
	    const REAL transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*RY*limiter(ym1-ym2,RY+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveYX;
	    dfdt[I0_J0_K1 + MYIND] += transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx < 0 Vy > 0 Vz > 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J2_K2 + MYIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J2_K1 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J1_K2 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J1_K2 + MYIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K2 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J2_K1 + MYIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransIncrWave*RZ;
	       dfdt[I1_J2_K1 + MYIND] -= doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xp1-xcc,RX+EPSILON,xcc);
	       dfdt[I1_J2_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(ym1-ym2,RY+EPSILON,xcc);
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K2 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc);
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K0 + MYIND] += doubleTransCorrWaveZXY;
               #endif
	    } else { // Vx < 0 Vy > 0 Vz < 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J2_K0 + MYIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J1_K0 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J2_K1 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J1_K0 + MYIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K0 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I0_J2_K0 + MYIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J2_K0 + MYIND] -= doubleTransIncrWave*RZ;
	       dfdt[I0_J1_K0 + MYIND] -= doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xp1-xcc,RX+EPSILON,xcc);
	       dfdt[I1_J2_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(ym1-ym2,RY+EPSILON,xcc);
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc);
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K0 + MYIND] += doubleTransCorrWaveZXY;
               #endif
	    }
	 } else { // Vx < 0 Vy < 0
	    const REAL transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J0_K1 + MYIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYIND] -= transIncrWave*RX;
	    dfdt[I0_J0_K1 + MYIND] += transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYIND] -= transIncrWave*RY;

	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*RX*limiter(xp1-xcc,RX+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveXY;
	    dfdt[I0_J0_K1 + MYIND] += transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveXY;
	    const REAL transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*RY*limiter(yp1-xcc,RY+EPSILON,xcc)*dt_per_dx*dt_per_dy;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveYX;
	    dfdt[I0_J0_K1 + MYIND] += transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx < 0 Vy < 0 Vz > 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J0_K2 + MYIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J1_K2 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J0_K1 + MYIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J0_K2 + MYIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J0_K2 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J0_K1 + MYIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J0_K1 + MYIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J0_K1 + MYIND] += doubleTransIncrWave*RZ;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xp1-xcc,RX+EPSILON,xcc);
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K2 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(yp1-xcc,RY+EPSILON,xcc);
	       dfdt[I1_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K2 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc);
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
               #endif
	    } else { // Vx < 0 Vy < 0 Vz < 0
	       const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J0_K0 + MYIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J1_K0 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J0_K0 + MYIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I0_J0_K1 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransIncrWave*RY;
	       dfdt[I0_J0_K0 + MYIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I0_J1_K0 + MYIND] += doubleTransIncrWave*RZ;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransIncrWave*RZ;
	       
	       #ifndef LEV_1ST_ORDER
	       const REAL doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*RX*limiter(xp1-xcc,RX+EPSILON,xcc);
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K0 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransCorrWaveXYZ;
	       const REAL doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*RY*limiter(yp1-xcc,RY+EPSILON,xcc);
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K0 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransCorrWaveYZX;
	       const REAL doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc);
	       dfdt[I1_J0_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYIND] += doubleTransCorrWaveZXY;
	       #endif
	    }
	 }
	 
	 if (Vz >= ZERO) {
	    const REAL transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J1_K2 + MYIND] -= transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYIND] -= transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*RX*limiter(xp1-xcc,RX+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I1_J1_K2 + MYIND] += transCorrWaveXZ;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K2 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYIND] += transCorrWaveXZ;
	    const REAL transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveZX;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveZX;
	    dfdt[I0_J1_K0 + MYIND] += transCorrWaveZX;
	    #endif	    
	 } else {
	    const REAL transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J1_K0 + MYIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYIND] -= transIncrWave*RX;
	    dfdt[I0_J1_K0 + MYIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYIND] -= transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*RX*limiter(xp1-xcc,RX+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveXZ;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K0 + MYIND] += transCorrWaveXZ;
	    const REAL transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc)*dt_per_dx*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveZX;
	    dfdt[I0_J1_K1 + MYIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveZX;
	    dfdt[I0_J1_K0 + MYIND] += transCorrWaveZX;
	    #endif
	 }
      }
      
      if (Vy >= ZERO) {
	 const REAL RY = xcc - ym1;
	 const REAL incrWaveY = Vy*ym1*dt_per_dy;
	 dfdt[I1_J1_K1 + MYIND] += incrWaveY;
	 dfdt[I1_J0_K1 + MYIND] -= incrWaveY;
	 #ifndef LEV_1ST_ORDER
	 const REAL corrWaveY = HALF*Vy*(ONE - dt_per_dy*Vy)*RY*limiter(ym1-ym2,RY+EPSILON,xcc)*dt_per_dy;
	 dfdt[I1_J1_K1 + MYIND] += corrWaveY;
	 dfdt[I1_J0_K1 + MYIND] -= corrWaveY;
         #endif
	 
	 const REAL RZ = xcc - zm1;
	 if (Vz >= ZERO) {
	    const REAL transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYIND] -= transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYIND] += transIncrWave*RY;
	    dfdt[I1_J2_K1 + MYIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYIND] += transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*RY*limiter(ym1-ym2,RY+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K2 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYIND] += transCorrWaveYZ;
	    const REAL transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J2_K1 + MYIND] += transCorrWaveZY;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J2_K0 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYIND] += transCorrWaveZY;
	    #endif
	 } else {
	    const REAL transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K0 + MYIND] += transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYIND] -= transIncrWave*RY;
	    dfdt[I1_J2_K0 + MYIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYIND] += transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*RY*limiter(ym1-ym2,RY+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K0 + MYIND] += transCorrWaveYZ;
	    const REAL transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J2_K1 + MYIND] += transCorrWaveZY;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J2_K0 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYIND] += transCorrWaveZY;
	    #endif
	 }
      } else { // Vy < 0
	 const REAL RY = xcc - ym1;
	 const REAL incrWaveY = Vy*xcc*dt_per_dy;
	 dfdt[I1_J1_K1 + MYIND] += incrWaveY;
	 dfdt[I1_J0_K1 + MYIND] -= incrWaveY;
	 #ifndef LEV_1ST_ORDER
	 const REAL corrWaveY = -HALF*Vy*(ONE + dt_per_dy*Vy)*RY*limiter(yp1-xcc,RY+EPSILON,xcc)*dt_per_dy;
	 dfdt[I1_J1_K1 + MYIND] += corrWaveY;
	 dfdt[I1_J0_K1 + MYIND] -= corrWaveY;
         #endif

	 const REAL RZ = xcc - zm1;
	 if (Vz >= ZERO) {
	    const REAL transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J0_K2 + MYIND] -= transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYIND] += transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYIND] -= transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*RY*limiter(yp1-xcc,RY+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K1 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K2 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYIND] += transCorrWaveYZ;
	    const REAL transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*RZ*limiter(zm1-zm2,RZ+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveZY;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J0_K0 + MYIND] += transCorrWaveZY;
	    #endif
	 } else {
	    const REAL transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J0_K0 + MYIND] += transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYIND] -= transIncrWave*RY;
	    dfdt[I1_J0_K0 + MYIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYIND] -= transIncrWave*RZ;
	    
	    #ifndef LEV_1ST_ORDER
	    const REAL transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*RY*limiter(yp1-xcc,RY+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K0 + MYIND] += transCorrWaveYZ;
	    const REAL transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc)*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K1 + MYIND] += transCorrWaveZY;
	    dfdt[I1_J0_K1 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYIND] -= transCorrWaveZY;
	    dfdt[I1_J0_K0 + MYIND] += transCorrWaveZY;
	    #endif
	 }
      }
      
      if (Vz >= ZERO) {
	 const REAL R = xcc - zm1;
	 const REAL incrWaveZ = Vz*zm1*dt_per_dz;
	 dfdt[I1_J1_K1 + MYIND] += incrWaveZ;
	 dfdt[I1_J1_K0 + MYIND] -= incrWaveZ;
	 #ifndef LEV_1ST_ORDER
	 const REAL corrWaveZ = HALF*Vz*(ONE - dt_per_dz*Vz)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dz;
	 dfdt[I1_J1_K1 + MYIND] += corrWaveZ;
	 dfdt[I1_J1_K0 + MYIND] -= corrWaveZ;
         #endif
      } else {
	 const REAL RZ = xcc - zm1;
	 const REAL incrWaveZ = Vz*xcc*dt_per_dz;
	 dfdt[I1_J1_K1 + MYIND] += incrWaveZ;
	 dfdt[I1_J1_K0 + MYIND] -= incrWaveZ;
	 #ifndef LEV_1ST_ORDER
	 const REAL corrWaveZ = -HALF*Vz*(ONE + dt_per_dz*Vz)*RZ*limiter(zp1-xcc,RZ+EPSILON,xcc)*dt_per_dz;
	 dfdt[I1_J1_K1 + MYIND] += corrWaveZ;
	 dfdt[I1_J1_K0 + MYIND] -= corrWaveZ;
         #endif
      }      
   }

   // Accumulate calculated df/dt values from temporary buffer to 
   // main memory. If multithreading is used, these updates need 
   // to be atomistic:
   const UINT boundaryFlags = nbrsSpa[BLOCK*SIZE_NBRS_SPA + 30];
   for (uint nbr=0; nbr<27; ++nbr) {
      // If the neighbour does not exist, do not copy data:
      if (((boundaryFlags >> nbr) & 1) == 0) continue;
      
      const UINT nbrBlock = nbrsSpa[BLOCK*SIZE_NBRS_SPA + nbr];
      for (uint i=0; i<SIZE_VELBLOCK; ++i) flux[nbrBlock*WID3 + i] += dfdt[nbr*WID3 + i];
   }
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
