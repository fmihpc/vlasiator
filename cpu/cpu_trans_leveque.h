#ifndef CPU_TRANS_H
#define CPU_TRANS_H

#include <vector>
#include "../definitions.h"
#include "../common.h"
#include "cpu_common.h"

const uint I0_J0_K0 = 0;
const uint I1_J0_K0 = 1;
const uint I2_J0_K0 = 2;
const uint I0_J1_K0 = 3;
const uint I1_J1_K0 = 4;
const uint I2_J1_K0 = 5;
const uint I0_J2_K0 = 6;
const uint I1_J2_K0 = 7;
const uint I2_J2_K0 = 8;

const uint I0_J0_K1 = 9;
const uint I1_J0_K1 = 10;
const uint I2_J0_K1 = 11;
const uint I0_J1_K1 = 12;
const uint I1_J1_K1 = 13;
const uint I2_J1_K1 = 14;
const uint I0_J2_K1 = 15;
const uint I1_J2_K1 = 16;
const uint I2_J2_K1 = 17;

const uint I0_J0_K2 = 18;
const uint I1_J0_K2 = 19;
const uint I2_J0_K2 = 20;
const uint I0_J1_K2 = 21;
const uint I1_J1_K2 = 22;
const uint I2_J1_K2 = 23;
const uint I0_J2_K2 = 24;
const uint I1_J2_K2 = 25;
const uint I2_J2_K2 = 26;

const Real EPSILON = 1.0e-25;
const Real ZERO    = 0.0;
const Real SIXTH   = 1.0/6.0;
const Real FOURTH  = 0.25;
const Real HALF    = 0.5;
const Real ONE     = 1.0;

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
   const REAL* nbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 12]*SIZE_VELBLOCK; //  -x nbr
   
   // Case Vx > 0:
   const REAL* nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 27]*SIZE_VELBLOCK; // --x nbr
   if (sign_vx_constant == false || vx_sign > ZERO) {
      for (UINT k=0; k<WID; ++k) {
	 const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
	 //const REAL Vz = ZERO;
	 for (UINT j=0; j<WID; ++j) {
	    const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	    //const REAL Vy = ZERO;
	    for (UINT i=0; i<WID; ++i) {
	       const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	       //const REAL Vx = ZERO;
	 
	       const UINT MYIND = cellIndex(i,j,k);
	       const REAL xcc = blockAvgs[MYIND];
	       const REAL xm1 = nbr_minus1[MYIND];
	       const REAL xm2 = nbr_minus2[MYIND];
	       
	       const REAL R = xcc - xm1;

	       // Increment waves:
	       const REAL incrWaveX = Vx*xm1*dt_per_dx;
	       dfdt[I1_J1_K1*WID3 + MYIND] += incrWaveX;
	       dfdt[I0_J1_K1*WID3 + MYIND] -= incrWaveX;
	       #ifndef LEV_1ST_ORDER
	       // Correction waves:
	       const REAL corrWaveX = HALF*Vx*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx;
	       dfdt[I1_J1_K1*WID3 + MYIND] += corrWaveX;
	       dfdt[I0_J1_K1*WID3 + MYIND] -= corrWaveX;
	       #endif
	       if (Vy >= ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= transIncrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= transIncrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  const REAL transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K2*WID3 + MYIND] -= doubleTransIncrWave; // Fy
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I1_J2_K2*WID3 + MYIND] -= doubleTransIncrWave; // Fz
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz 
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  #endif
	       } else if (Vy >= ZERO && Vz < ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= transIncrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  const REAL transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransIncrWave; // Fy
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransIncrWave; // Fz
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = +FOURTH*Vx*Vy*Vz*(ONE - Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave; // Fy
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;

		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  #endif
	       } else if (Vy < ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= transIncrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  const REAL transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransIncrWave; // Fy
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransIncrWave; // Fz
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  const REAL transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransIncrWave; // Fy
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransIncrWave; // Fz
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       }
	    }
	 }
      }
   }

   // Case Vx < 0:
   const REAL* nbr_plus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 14]*SIZE_VELBLOCK; //  +x nbr
   if (sign_vx_constant == false || vx_sign < ZERO) {
      for (UINT k=0; k<WID; ++k) {
	 const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
	 //const REAL Vz = ZERO;
	 for (UINT j=0; j<WID; ++j) {
	    const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	    //const REAL Vy = ZERO;
	    for (UINT i=0; i<WID; ++i) {
	       const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	       //const REAL Vx = ZERO;
	    
	       const UINT MYIND = cellIndex(i,j,k);
	       const REAL xcc = blockAvgs[MYIND];
	       const REAL xm1 = nbr_minus1[MYIND];
	       const REAL xm2 = nbr_plus1[MYIND];
	       
	       const REAL R = xcc - xm1;
	       //const REAL corr_wave = -HALF*Vx*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx;
	       
	       // Increment waves:
	       const REAL incrWaveX = Vx*xcc*dt_per_dx;
	       dfdt[I1_J1_K1*WID3 + MYIND] += incrWaveX;
	       dfdt[I0_J1_K1*WID3 + MYIND] -= incrWaveX;
	       #ifndef LEV_1ST_ORDER
	       // Correction waves:
	       const REAL corrWaveX = -HALF*Vx*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx;
	       dfdt[I1_J1_K1*WID3 + MYIND] += corrWaveX;
	       dfdt[I0_J1_K1*WID3 + MYIND] -= corrWaveX;
	       #endif
	       if (Vy >= ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= transIncrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= transIncrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  const REAL transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J2_K2*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J2_K2*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransCorrWave;

		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else if (Vy >= ZERO && Vz < ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= transIncrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transIncrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  const REAL transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransIncrWave; // Fy
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransIncrWave; // Fz
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else if (Vy < ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transIncrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= transIncrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  const REAL transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransIncrWave; // Fy
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransIncrWave; // Fz
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else {
		  // Transverse increment waves:
		  const REAL transIncrWaveXY = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transIncrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transIncrWaveXY;
		  const REAL transIncrWaveXZ = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transIncrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transIncrWaveXZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveXY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveXY;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXY;
		  const REAL transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveXZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveXZ;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveXZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransIncrWave; // Fy
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransIncrWave; // Fz
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       }
	    }
	 }
      }
   }

   // ***** Consider the interface between (i,j-1,k) and (i,j,k): *****
   nbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 10]*SIZE_VELBLOCK;

   // Case Vy > 0:
   nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 28]*SIZE_VELBLOCK; // --y nbr
   if (sign_vy_constant == false || vy_sign > ZERO) {
      for (UINT k=0; k<WID; ++k) {
	 const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
	 //const REAL Vz = ZERO;
	 for (UINT i=0; i<WID; ++i) {
	    const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	    //const REAL Vx = ZERO;
	    for (UINT j=0; j<WID; ++j) {
	       const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	       //const REAL Vy = ZERO;
	    
	       const UINT MYIND = cellIndex(i,j,k);
	       const REAL xcc = blockAvgs[MYIND];
	       const REAL xm1 = nbr_minus1[MYIND];
	       const REAL xm2 = nbr_minus2[MYIND];
	    
	       const REAL R = xcc - xm1;
	       //const REAL corr_wave = HALF*Vy*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy;
	    
	       // Increment waves:
	       const REAL incrWaveY = Vy*xm1*dt_per_dy;
	       dfdt[I1_J1_K1*WID3 + MYIND] += incrWaveY;
	       dfdt[I1_J0_K1*WID3 + MYIND] -= incrWaveY;
	       #ifndef LEV_1ST_ORDER
	       // Correction waves:
	       const REAL corrWaveY = HALF*Vy*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy;
	       dfdt[I1_J1_K1*WID3 + MYIND] += corrWaveY;
	       dfdt[I1_J0_K1*WID3 + MYIND] -= corrWaveY;
	       #endif
	       if (Vx >= ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= transIncrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= transIncrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  const REAL transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K2*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I2_J1_K2*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  #endif
	       } else if (Vx >= ZERO && Vz < ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= transIncrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  const REAL transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave; // Fx
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  #endif
	       } else if (Vx < ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= transIncrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  const REAL transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransIncrWave;
                  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  const REAL transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       }
	    }
	 }
      }
   }

   // Case Vy < 0:
   nbr_plus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 16]*SIZE_VELBLOCK; //  +y nbr
   if (sign_vy_constant == false || vy_sign < ZERO) {
      for (UINT k=0; k<WID; ++k) {
	 const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
	 //const REAL Vz = ZERO;
	 for (UINT i=0; i<WID; ++i) {
	    const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	    //const REAL Vx = ZERO;
	    for (UINT j=0; j<WID; ++j) {
	       const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	       //const REAL Vy = ZERO;
	    
	       const UINT MYIND = cellIndex(i,j,k);
	       const REAL xcc = blockAvgs[MYIND];
	       const REAL xm1 = nbr_minus1[MYIND];
	       const REAL xm2 = nbr_plus1[MYIND];
	       
	       const REAL R = xcc - xm1;
	       //const REAL corr_wave = -HALF*Vy*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy;

	       // Increment waves:
	       const REAL incrWaveY = Vy*xcc*dt_per_dy;
	       dfdt[I1_J1_K1*WID3 + MYIND] += incrWaveY;
	       dfdt[I1_J0_K1*WID3 + MYIND] -= incrWaveY;
	       #ifndef LEV_1ST_ORDER
	       // Correction waves:
	       const REAL corrWaveY = -HALF*Vy*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy;
	       dfdt[I1_J1_K1*WID3 + MYIND] += corrWaveY;
	       dfdt[I1_J0_K1*WID3 + MYIND] -= corrWaveY;
	       #endif
	       if (Vx >= ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= transIncrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= transIncrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  const REAL transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J0_K2*WID3 + MYIND] -= doubleTransIncrWave; // Fx
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I2_J0_K2*WID3 + MYIND] -= doubleTransIncrWave; // FZ
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else if (Vx >= ZERO && Vz < ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= transIncrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transIncrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  const REAL transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransIncrWave; // Fx
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransIncrWave; // Fz
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else if (Vx < ZERO && Vz >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transIncrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= transIncrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  const REAL transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K2*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransIncrWave; // Fx
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransIncrWave; // Fz
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J1_K2*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K2*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K2*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else {
		  // Transverse increment waves:
		  const REAL transIncrWaveYX = HALF*Vx*Vy*dt_per_dx*dt_per_dy*R - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transIncrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transIncrWaveYX;
		  const REAL transIncrWaveYZ = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transIncrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transIncrWaveYZ;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dy;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveYX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I0_J0_K1*WID3 + MYIND] += transCorrWaveYX;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYX;
		  const REAL transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveYZ;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveYZ;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveYZ;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransIncrWave; // Fx
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransIncrWave; // Fz
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fz
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       }
	    }
	 }
      }
   }

   // ***** Consider the interface between (i,j,k-1) and (i,j,k): *****
   nbr_minus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 4]*SIZE_VELBLOCK;

   // Case Vz > 0:
   if (sign_vz_constant == false || vz_sign > ZERO) {
      nbr_minus2 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 29]*SIZE_VELBLOCK; // --z nbr
      for (UINT j=0; j<WID; ++j) {
	 const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	 //const REAL Vy = ZERO;
	 for (UINT i=0; i<WID; ++i) {
	    const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	    //const REAL Vx = ZERO;
	    for (UINT k=0; k<WID; ++k) {
	       const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
	       //const REAL Vz = ZERO;
	       
	       const UINT MYIND = cellIndex(i,j,k);
	       const REAL xcc = blockAvgs[MYIND];
	       const REAL xm1 = nbr_minus1[MYIND];
	       const REAL xm2 = nbr_minus2[MYIND];

	       const REAL R = xcc - xm1;
	       
	       // Increment waves:
	       const REAL incrWaveZ = Vz*xm1*dt_per_dz;
	       dfdt[I1_J1_K1*WID3 + MYIND] += incrWaveZ;
	       dfdt[I1_J1_K0*WID3 + MYIND] -= incrWaveZ;
	       #ifndef LEV_1ST_ORDER
	       // Correction waves:
	       const REAL corrWaveZ = HALF*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dz;
	       dfdt[I1_J1_K1*WID3 + MYIND] += corrWaveZ;
	       dfdt[I1_J1_K0*WID3 + MYIND] -= corrWaveZ;
	       #endif
	       if (Vx >= ZERO && Vy >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= transIncrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= transIncrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double Transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J2_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I2_J2_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;		  
		  dfdt[I2_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  #endif
	       } else if (Vx >= ZERO && Vy < ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= transIncrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double Transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else if (Vx < ZERO && Vy >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= transIncrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double Transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;

		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K1*WID3 + MYIND] += transIncrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K1*WID3 + MYIND] += transIncrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double Transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J0_K1*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = FOURTH*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm1-xm2,R+EPSILON,xcc);
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       }
	    }
	 }
      }
   }
   
   // Case Vz < 0:
   if (sign_vz_constant == false || vz_sign < ZERO) {
      nbr_plus1 = AVGS + nbrsSpa[BLOCK*SIZE_NBRS_SPA + 22]*SIZE_VELBLOCK; // +z nbr
      for (UINT j=0; j<WID; ++j) {
	 const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
	 //const REAL Vy = ZERO;
	 for (UINT i=0; i<WID; ++i) {
	    const REAL Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
	    //const REAL Vx = ZERO;
	    for (UINT k=0; k<WID; ++k) {
	       const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
	       //const REAL Vz = ZERO;
	       
	       const UINT MYIND = cellIndex(i,j,k);
	       const REAL xcc = blockAvgs[MYIND];
	       const REAL xm1 = nbr_minus1[MYIND];
	       const REAL xm2 = nbr_plus1[MYIND];
	       
	       const REAL R = xcc - xm1;
	       
	       // Increment waves:
	       const REAL incrWaveZ = Vz*xcc*dt_per_dz;
	       dfdt[I1_J1_K1*WID3 + MYIND] += incrWaveZ;
	       dfdt[I1_J1_K0*WID3 + MYIND] -= incrWaveZ;
	       #ifndef LEV_1ST_ORDER
	       const REAL corrWaveZ = -HALF*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dz;
	       dfdt[I1_J1_K1*WID3 + MYIND] += corrWaveZ;
	       dfdt[I1_J1_K0*WID3 + MYIND] -= corrWaveZ;
               #endif
	       
	       if (Vx >= ZERO && Vy >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= transIncrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= transIncrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZY;
                  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J2_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I2_J2_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransCorrWave;

		  dfdt[I1_J2_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  #endif
	       } else if (Vx >= ZERO && Vy < ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= transIncrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transIncrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I2_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I1_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I2_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I2_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else if (Vx < ZERO && Vy >= ZERO) {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transIncrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= transIncrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J2_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J2_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I0_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J2_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J2_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       } else {
		  // Transverse increment waves:
		  const REAL transIncrWaveZX = HALF*Vx*Vz*dt_per_dx*dt_per_dz*R - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transIncrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transIncrWaveZX;
		  const REAL transIncrWaveZY = HALF*Vy*Vz*dt_per_dy*dt_per_dz*R - SIXTH*fabs(Vx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transIncrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transIncrWaveZY;
		  #ifndef LEV_1ST_ORDER
		  // Transverse correction waves:
		  const REAL transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dx*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZX;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZX;
		  dfdt[I0_J1_K0*WID3 + MYIND] += transCorrWaveZX;
		  const REAL transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*R*limiter(xm2-xcc,R+EPSILON,xcc)*dt_per_dy*dt_per_dz;
		  dfdt[I1_J1_K1*WID3 + MYIND] += transCorrWaveZY;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= transCorrWaveZY;
		  dfdt[I1_J0_K0*WID3 + MYIND] += transCorrWaveZY;
		  #endif
		  // Double transverse increment waves:
		  const REAL doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz*R;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransIncrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransIncrWave;
		  #ifndef LEV_1ST_ORDER
		  // Double transverse correction waves:
		  const REAL doubleTransCorrWave = -FOURTH*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*dt_per_dz*R*limiter(xm2-xcc,R+EPSILON,xcc);
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fx
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  
		  dfdt[I0_J1_K1*WID3 + MYIND] -= doubleTransCorrWave; // Fy
		  dfdt[I0_J0_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J1_K0*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I0_J0_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K1*WID3 + MYIND] += doubleTransCorrWave;
		  dfdt[I1_J0_K1*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J1_K0*WID3 + MYIND] -= doubleTransCorrWave;
		  dfdt[I1_J0_K0*WID3 + MYIND] += doubleTransCorrWave;
                  #endif
	       }
	    }
	 }
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
