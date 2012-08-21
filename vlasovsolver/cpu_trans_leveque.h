/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CPU_TRANS_H
#define CPU_TRANS_H

#include <vector>
#include <limits>
#include <dccrg.hpp>

#include "../definitions.h"
#include "../common.h"
#include "limiters_vec4.h"
#include "leveque_common.h"
#include "spatial_cell.hpp"
#include "vec4.h"



const uint I0_J0_K0 = 0*WID2;
const uint I1_J0_K0 = 1*WID2;
const uint I2_J0_K0 = 2*WID2;
const uint I0_J1_K0 = 3*WID2;
const uint I1_J1_K0 = 4*WID2;
const uint I2_J1_K0 = 5*WID2;
const uint I0_J2_K0 = 6*WID2;
const uint I1_J2_K0 = 7*WID2;
const uint I2_J2_K0 = 8*WID2;

const uint I0_J0_K1 = 9*WID2;
const uint I1_J0_K1 = 10*WID2;
const uint I2_J0_K1 = 11*WID2;
const uint I0_J1_K1 = 12*WID2;
const uint I1_J1_K1 = 13*WID2;
const uint I2_J1_K1 = 14*WID2;
const uint I0_J2_K1 = 15*WID2;
const uint I1_J2_K1 = 16*WID2;
const uint I2_J2_K1 = 17*WID2;

const uint I0_J0_K2 = 18*WID2;
const uint I1_J0_K2 = 19*WID2;
const uint I2_J0_K2 = 20*WID2;
const uint I0_J1_K2 = 21*WID2;
const uint I1_J1_K2 = 22*WID2;
const uint I2_J1_K2 = 23*WID2;
const uint I0_J2_K2 = 24*WID2;
const uint I1_J2_K2 = 25*WID2;
const uint I2_J2_K2 = 26*WID2;



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
      const REAL VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      const REAL VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const REAL VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
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
template<typename REAL,typename UINT> void cpu_calcSpatDfdt(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid, spatial_cell::SpatialCell* cell,const UINT& blockId,const REAL& dt) {
//const REAL* const AVGS,cons        t REAL* const CELL_PARAMS,const REAL* const BLOCK_PARAMS,REAL* const flux,
//					        		    const UINT* const nbrsSpa

   spatial_cell::Velocity_Block* block=cell->at(blockId); //returns a reference to block
// Create a temporary buffer for storing df/dt updates and init to zero value:

   if(WID!=4){
      //hard-limit from vectorization
      std::cerr<<"ABORT, cpu_calcSpatDfdt assumes blockwidth=4" << std::endl;
      exit(1);
   }
   
   const UINT SIZE_FLUXBUFFER = 27*WID2;
   Vec4 dfdt[SIZE_FLUXBUFFER] = {};
   
   for(UINT i=0;i<SIZE_FLUXBUFFER;i++){
      dfdt[i]=Vec4(0);
   }
//  Pointer to velocity block whose df/dt contribution        s are calculated:
   const REAL*__restrict__ const blockAvgs   = block->data;
   const REAL*__restrict__ const blockParams = block->parameters;
   
   // ***** Consider the interface between (i-1,j,k) and (i,j,k): *****
   const REAL dt_per_dx = dt / cell->parameters[CellParams::DX];
   const REAL dt_per_dy = dt / cell->parameters[CellParams::DY];
   const REAL dt_per_dz = dt / cell->parameters[CellParams::DZ];

//FIXME, the indices could have a enum or namespace
   const REAL*__restrict__ const xnbr_plus1  = mpiGrid[cell->neighbors[14]]->at(blockId)->data;
   const REAL*__restrict__ const xnbr_minus1 = mpiGrid[cell->neighbors[12]]->at(blockId)->data;
   const REAL*__restrict__ const xnbr_minus2 = mpiGrid[cell->neighbors[27]]->at(blockId)->data; // --x nbr
   const REAL*__restrict__ const ynbr_plus1  = mpiGrid[cell->neighbors[16]]->at(blockId)->data; //  +y nbr
   const REAL*__restrict__ const ynbr_minus1 = mpiGrid[cell->neighbors[10]]->at(blockId)->data;
   const REAL*__restrict__ const ynbr_minus2 = mpiGrid[cell->neighbors[28]]->at(blockId)->data; // --y nbr
   const REAL*__restrict__ const znbr_plus1  = mpiGrid[cell->neighbors[22]]->at(blockId)->data; // +z nbr
   const REAL*__restrict__ const znbr_minus1 = mpiGrid[cell->neighbors[4]]->at(blockId)->data;
   const REAL*__restrict__ const znbr_minus2 = mpiGrid[cell->neighbors[29]]->at(blockId)->data; // --z nbr


   Vec4 allVx(0.0);
   Vec4 posVx(0.0);
   Vec4 negVx(0.0);
   bool hasPosVx=false;
   bool hasNegVx=false;   
   for(UINT i=0;i<WID;++i){
      const REAL VxVal=blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      allVx.insert(i,VxVal);
      
      if(VxVal>=0){
         posVx.insert(i,VxVal);
         hasPosVx=true;
      }
      else{
         negVx.insert(i,VxVal);
         hasNegVx=true;
      }
   }
   
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
      const UINT MYVECIND =k*WID + j;
      const UINT MYIND = MYVECIND*WID;
      const REAL Vz = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      const REAL Vy = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      
      // Fetch required values from global arrays:


      Vec4 xcc;
      Vec4 xp1,xm1,xm2;
      Vec4 yp1,ym1,ym2;
      Vec4 zp1,zm1,zm2;
      xcc.load(blockAvgs+MYIND);
      xp1.load(xnbr_plus1+MYIND);
      xm1.load(xnbr_minus1+MYIND);
      xm2.load(xnbr_minus2+MYIND);
      yp1.load(ynbr_plus1+MYIND);
      ym1.load(ynbr_minus1+MYIND);
      ym2.load(ynbr_minus2+MYIND);
      zp1.load(znbr_plus1+MYIND);
      zm1.load(znbr_minus1+MYIND);
      zm2.load(znbr_minus2+MYIND);

      const Vec4 RX = xcc - xm1;
      const Vec4 RY = xcc - ym1;
      const Vec4 RZ = xcc - zm1;
      const Vec4 RX_limiter_xm1xm2_dtdx=RX*limiter_vec4(xm1-xm2,RX+EPSILON,xcc)*dt_per_dx;
      const Vec4 RY_limiter_ym1ym2_dtdy=RY*limiter_vec4(ym1-ym2,RY+EPSILON,xcc)*dt_per_dy;
      const Vec4 RZ_limiter_zm1zm2_dtdz=RZ*limiter_vec4(zm1-zm2,RZ+EPSILON,xcc)*dt_per_dz;
      const Vec4 RX_limiter_xp1xcc_dtdx=RX*limiter_vec4(xp1-xcc,RX+EPSILON,xcc)*dt_per_dx;
      const Vec4 RY_limiter_yp1xcc_dtdy=RY*limiter_vec4(yp1-xcc,RY+EPSILON,xcc)*dt_per_dy;
      const Vec4 RZ_limiter_zp1xcc_dtdz=RZ*limiter_vec4(zp1-xcc,RZ+EPSILON,xcc)*dt_per_dz;

         
      if (hasPosVx) {
         // 1D increment and correction waves:
         const Vec4 Vx=posVx;

         const Vec4 incrWaveX = Vx*xm1*dt_per_dx;
	 dfdt[I1_J1_K1 + MYVECIND] += incrWaveX;
	 dfdt[I0_J1_K1 + MYVECIND] -= incrWaveX;
	 #if MOVER_VLASOV_ORDER > 1
	 const Vec4 corrWaveX = HALF*Vx*(ONE - dt_per_dx*Vx)*RX_limiter_xm1xm2_dtdx;
	 dfdt[I1_J1_K1 + MYVECIND] += corrWaveX;
	 dfdt[I0_J1_K1 + MYVECIND] -= corrWaveX;
         #endif
	 

	 if (Vy >= ZERO) { // Vx > 0 Vy > 0
	    const Vec4 transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J2_K1 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYVECIND] += transIncrWave*RX;
	    dfdt[I2_J1_K1 + MYVECIND] -= transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYVECIND] += transIncrWave*RY;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*RX_limiter_xm1xm2_dtdx*dt_per_dy;
	    dfdt[I1_J2_K1 + MYVECIND] += transCorrWaveXY;
	    dfdt[I0_J2_K1 + MYVECIND] -= transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYVECIND] += transCorrWaveXY;
	    const Vec4 transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*RY_limiter_ym1ym2_dtdy*dt_per_dx;
	    dfdt[I2_J1_K1 + MYVECIND] += transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveYX;
	    dfdt[I2_J0_K1 + MYVECIND] -= transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYVECIND] += transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx > 0 Vy > 0 Vz > 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J2_K2 + MYVECIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J1_K2 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I1_J2_K1 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I2_J1_K2 + MYVECIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K2 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J1_K1 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J2_K1 + MYVECIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J2_K1 + MYVECIND] += doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K1 + MYVECIND] += doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE - dt_per_dx*Vx)*dt_per_dy*dt_per_dz*RX_limiter_xm1xm2_dtdx;
	       dfdt[I1_J2_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dz*RY_limiter_ym1ym2_dtdy;
	       dfdt[I2_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zm1zm2_dtdz;
	       dfdt[I2_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       #endif
	    } else { // Vx > 0 Vy > 0 Vz < 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J2_K0 + MYVECIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I1_J2_K1 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I2_J1_K0 + MYVECIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J1_K1 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J2_K0 + MYVECIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J2_K0 + MYVECIND] += doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K0 + MYVECIND] += doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE - Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xm1xm2_dtdx;
	       dfdt[I1_J2_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dz*RY_limiter_ym1ym2_dtdy;
	       dfdt[I2_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zp1xcc_dtdz;
	       dfdt[I2_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J2_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       #endif
	    }
	 } else { // Vx > 0 Vy < 0
	    const Vec4 transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J0_K1 + MYVECIND] += transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I2_J0_K1 + MYVECIND] -= transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYVECIND] += transIncrWave*RY;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXY = HALF*Vx*Vy*(ONE - dt_per_dx*Vx)*RX_limiter_xm1xm2_dtdx*dt_per_dy;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveXY;
	    dfdt[I0_J0_K1 + MYVECIND] += transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveXY;
	    const Vec4 transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*RY_limiter_yp1xcc_dtdy*dt_per_dx;
	    dfdt[I2_J1_K1 + MYVECIND] += transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveYX;
	    dfdt[I2_J0_K1 + MYVECIND] -= transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYVECIND] += transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx > 0 Vy < 0 Vz > 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J0_K2 + MYVECIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J1_K2 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I2_J0_K2 + MYVECIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J0_K2 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J0_K1 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I2_J0_K1 + MYVECIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K1 + MYVECIND] -= doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xm1xm2_dtdx;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*RY_limiter_yp1xcc_dtdy*dt_per_dz;
	       dfdt[I2_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zm1zm2_dtdz;
	       dfdt[I2_J0_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J0_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       #endif
	    } else { // Vx > 0 Vy < 0 Vz < 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I1_J0_K0 + MYVECIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I2_J0_K0 + MYVECIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J0_K0 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J0_K1 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I2_J0_K0 + MYVECIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J0_K0 + MYVECIND] -= doubleTransIncrWave*RZ;
	       dfdt[I2_J1_K0 + MYVECIND] -= doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = HALF*Vx*Vy*Vz*(ONE-Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xm1xm2_dtdx;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*RY_limiter_yp1xcc_dtdy*dt_per_dz;
	       dfdt[I2_J1_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I2_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zp1xcc_dtdz;
	       dfdt[I2_J0_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I2_J0_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I2_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       #endif
	    }
	 }
	 
	 if (Vz >= ZERO) {
	    const Vec4 transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYVECIND] += transIncrWave*RX;
	    dfdt[I2_J1_K1 + MYVECIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYVECIND] += transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*RX_limiter_xm1xm2_dtdx*dt_per_dz;
	    dfdt[I1_J1_K2 + MYVECIND] += transCorrWaveXZ;
	    dfdt[I0_J1_K2 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYVECIND] += transCorrWaveXZ;
	    const Vec4 transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*RZ_limiter_zm1zm2_dtdz*dt_per_dx;
	    dfdt[I2_J1_K1 + MYVECIND] += transCorrWaveZX;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I2_J1_K0 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYVECIND] += transCorrWaveZX;
	    #endif
	 } else {
	    const Vec4 transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K0 + MYVECIND] += transIncrWave*RX;
	    dfdt[I1_J1_K1 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I2_J1_K0 + MYVECIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYVECIND] += transIncrWave*RZ;

	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXZ = HALF*Vx*Vz*(ONE - dt_per_dx*Vx)*RX_limiter_xm1xm2_dtdx*dt_per_dz;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveXZ;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K0 + MYVECIND] += transCorrWaveXZ;
	    const Vec4 transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*RZ_limiter_zp1xcc_dtdz*dt_per_dx;
	    dfdt[I2_J1_K1 + MYVECIND] += transCorrWaveZX;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I2_J1_K0 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYVECIND] += transCorrWaveZX;
            #endif
	 }
         
         
      }

      if (hasNegVx) {
            // 1D increment and correction waves:
         const Vec4 Vx=negVx;

         const Vec4 incrWaveX = Vx*xcc*dt_per_dx;
	 dfdt[I1_J1_K1 + MYVECIND] += incrWaveX;
	 dfdt[I0_J1_K1 + MYVECIND] -= incrWaveX;
	 #if MOVER_VLASOV_ORDER > 1
	 const Vec4 corrWaveX = -HALF*Vx*(ONE + dt_per_dx*Vx)*RX_limiter_xp1xcc_dtdx;
	 dfdt[I1_J1_K1 + MYVECIND] += corrWaveX;
	 dfdt[I0_J1_K1 + MYVECIND] -= corrWaveX;
	 #endif
  
	 if (Vy >= ZERO) { // Vx < 0 Vy > 0
	    const Vec4 transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J2_K1 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYVECIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYVECIND] += transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYVECIND] -= transIncrWave*RY;

	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*RX_limiter_xp1xcc_dtdx*dt_per_dy;
	    dfdt[I1_J2_K1 + MYVECIND] += transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveXY;
	    dfdt[I0_J2_K1 + MYVECIND] -= transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYVECIND] += transCorrWaveXY;
	    const Vec4 transCorrWaveYX = HALF*Vx*Vy*(ONE - dt_per_dy*Vy)*RY_limiter_ym1ym2_dtdy*dt_per_dx;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveYX;
	    dfdt[I0_J0_K1 + MYVECIND] += transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx < 0 Vy > 0 Vz > 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J2_K2 + MYVECIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J2_K1 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J1_K2 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J1_K2 + MYVECIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K2 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J2_K1 + MYVECIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransIncrWave*RZ;
	       dfdt[I1_J2_K1 + MYVECIND] -= doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xp1xcc_dtdx;
	       dfdt[I1_J2_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dz*RY_limiter_ym1ym2_dtdy;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zm1zm2_dtdz;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K0 + MYVECIND] += doubleTransCorrWaveZXY;
               #endif
	    } else { // Vx < 0 Vy > 0 Vz < 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J2_K0 + MYVECIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J1_K0 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J2_K1 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J1_K0 + MYVECIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I0_J2_K0 + MYVECIND] += TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J2_K0 + MYVECIND] -= doubleTransIncrWave*RZ;
	       dfdt[I0_J1_K0 + MYVECIND] -= doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xp1xcc_dtdx;
	       dfdt[I1_J2_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J2_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J2_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = HALF*Vx*Vy*Vz*(ONE - dt_per_dy*Vy)*dt_per_dx*dt_per_dz*RY_limiter_ym1ym2_dtdy;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zp1xcc_dtdz;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J2_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J2_K0 + MYVECIND] += doubleTransCorrWaveZXY;
               #endif
	    }
	 } else { // Vx < 0 Vy < 0
	    const Vec4 transIncrWave = HALF*Vx*Vy*dt_per_dx*dt_per_dy - SIXTH*Vx*Vy*fabs(Vz)*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J0_K1 + MYVECIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I0_J0_K1 + MYVECIND] += transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYVECIND] -= transIncrWave*RY;

	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXY = -HALF*Vx*Vy*(ONE + dt_per_dx*Vx)*RX_limiter_xp1xcc_dtdx*dt_per_dy;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveXY;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveXY;
	    dfdt[I0_J0_K1 + MYVECIND] += transCorrWaveXY;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveXY;
	    const Vec4 transCorrWaveYX = -HALF*Vx*Vy*(ONE + dt_per_dy*Vy)*RY_limiter_yp1xcc_dtdy*dt_per_dx;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveYX;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveYX;
	    dfdt[I0_J0_K1 + MYVECIND] += transCorrWaveYX;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveYX;
	    #endif
	    
	    if (Vz >= ZERO) { // Vx < 0 Vy < 0 Vz > 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J0_K2 + MYVECIND] += TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J1_K2 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J0_K1 + MYVECIND] -= doubleTransIncrWave*RX;
	       dfdt[I0_J0_K2 + MYVECIND] += TWO*doubleTransIncrWave*RY; // y
	       dfdt[I1_J0_K2 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J0_K1 + MYVECIND] -= doubleTransIncrWave*RY;
	       dfdt[I0_J0_K1 + MYVECIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransIncrWave*RZ;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xp1xcc_dtdx;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K2 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K2 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*RY_limiter_yp1xcc_dtdy*dt_per_dz;
	       dfdt[I1_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K2 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K2 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = HALF*Vx*Vy*Vz*(ONE - dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zm1zm2_dtdz;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
               #endif
	    } else { // Vx < 0 Vy < 0 Vz < 0
	       const Vec4 doubleTransIncrWave = SIXTH*Vx*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	       dfdt[I0_J0_K0 + MYVECIND] -= TWO*doubleTransIncrWave*RX; // x
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransIncrWave*RX;
	       dfdt[I0_J0_K0 + MYVECIND] -= TWO*doubleTransIncrWave*RY; // y
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransIncrWave*RY;
	       dfdt[I0_J0_K0 + MYVECIND] -= TWO*doubleTransIncrWave*RZ; // z
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransIncrWave*RZ;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransIncrWave*RZ;
	       
	       #if MOVER_VLASOV_ORDER > 1
	       const Vec4 doubleTransCorrWaveXYZ = -HALF*Vx*Vy*Vz*(ONE+Vx*dt_per_dx)*dt_per_dy*dt_per_dz*RX_limiter_xp1xcc_dtdx;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K0 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveXYZ;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransCorrWaveXYZ;
	       const Vec4 doubleTransCorrWaveYZX = -HALF*Vx*Vy*Vz*(ONE + dt_per_dy*Vy)*dt_per_dx*RY_limiter_yp1xcc_dtdy*dt_per_dz;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K0 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveYZX;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransCorrWaveYZX;
	       const Vec4 doubleTransCorrWaveZXY = -HALF*Vx*Vy*Vz*(ONE + dt_per_dz*Vz)*dt_per_dx*dt_per_dy*RZ_limiter_zp1xcc_dtdz;
	       dfdt[I1_J0_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I1_J0_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J0_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K1 + MYVECIND] += doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K1 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I1_J1_K0 + MYVECIND] -= doubleTransCorrWaveZXY;
	       dfdt[I0_J1_K0 + MYVECIND] += doubleTransCorrWaveZXY;
	       #endif
	    }
	 }
	 
	 if (Vz >= ZERO) {
	    const Vec4 transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J1_K2 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYVECIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYVECIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYVECIND] -= transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*RX_limiter_xp1xcc_dtdx*dt_per_dz;
	    dfdt[I1_J1_K2 + MYVECIND] += transCorrWaveXZ;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K2 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYVECIND] += transCorrWaveXZ;
	    const Vec4 transCorrWaveZX = HALF*Vx*Vz*(ONE - dt_per_dz*Vz)*RZ_limiter_zm1zm2_dtdz*dt_per_dx;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveZX;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I0_J1_K0 + MYVECIND] += transCorrWaveZX;
	    #endif	    
	 } else {
	    const Vec4 transIncrWave = HALF*Vx*Vz*dt_per_dx*dt_per_dz - SIXTH*Vx*fabs(Vy)*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I0_J1_K0 + MYVECIND] += transIncrWave*RX;
	    dfdt[I0_J1_K1 + MYVECIND] -= transIncrWave*RX;
	    dfdt[I0_J1_K0 + MYVECIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYVECIND] -= transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveXZ = -HALF*Vx*Vz*(ONE + dt_per_dx*Vx)*RX_limiter_xp1xcc_dtdx*dt_per_dz;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveXZ;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveXZ;
	    dfdt[I0_J1_K0 + MYVECIND] += transCorrWaveXZ;
	    const Vec4 transCorrWaveZX = -HALF*Vx*Vz*(ONE + dt_per_dz*Vz)*RZ_limiter_zp1xcc_dtdz*dt_per_dx;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveZX;
	    dfdt[I0_J1_K1 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveZX;
	    dfdt[I0_J1_K0 + MYVECIND] += transCorrWaveZX;
	    #endif
	 }

      }
      
      if (Vy >= ZERO) {
	 const Vec4 incrWaveY = Vy*ym1*dt_per_dy;
	 dfdt[I1_J1_K1 + MYVECIND] += incrWaveY;
	 dfdt[I1_J0_K1 + MYVECIND] -= incrWaveY;
	 #if MOVER_VLASOV_ORDER > 1
	 const Vec4 corrWaveY = HALF*Vy*(ONE - dt_per_dy*Vy)*RY_limiter_ym1ym2_dtdy;
	 dfdt[I1_J1_K1 + MYVECIND] += corrWaveY;
	 dfdt[I1_J0_K1 + MYVECIND] -= corrWaveY;
         #endif
	 
	 if (Vz >= ZERO) {
	    const Vec4 transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*abs(allVx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYVECIND] -= transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYVECIND] += transIncrWave*RY;
	    dfdt[I1_J2_K1 + MYVECIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYVECIND] += transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*RY_limiter_ym1ym2_dtdy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYVECIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K2 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYVECIND] += transCorrWaveYZ;
	    const Vec4 transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*RZ_limiter_zm1zm2_dtdz*dt_per_dy;
	    dfdt[I1_J2_K1 + MYVECIND] += transCorrWaveZY;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J2_K0 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYVECIND] += transCorrWaveZY;
	    #endif
	 } else {
	    const Vec4 transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*abs(allVx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J1_K0 + MYVECIND] += transIncrWave*RY;
	    dfdt[I1_J1_K1 + MYVECIND] -= transIncrWave*RY;
	    dfdt[I1_J2_K0 + MYVECIND] -= transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYVECIND] += transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveYZ = HALF*Vy*Vz*(ONE - dt_per_dy*Vy)*RY_limiter_ym1ym2_dtdy*dt_per_dz;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K0 + MYVECIND] += transCorrWaveYZ;
	    const Vec4 transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*RZ_limiter_zp1xcc_dtdz*dt_per_dy;
	    dfdt[I1_J2_K1 + MYVECIND] += transCorrWaveZY;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J2_K0 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYVECIND] += transCorrWaveZY;
	    #endif
	 }
      } else { // Vy < 0
	 const Vec4 incrWaveY = Vy*xcc*dt_per_dy;
	 dfdt[I1_J1_K1 + MYVECIND] += incrWaveY;
	 dfdt[I1_J0_K1 + MYVECIND] -= incrWaveY;
	 #if MOVER_VLASOV_ORDER > 1
	 const Vec4 corrWaveY = -HALF*Vy*(ONE + dt_per_dy*Vy)*RY_limiter_yp1xcc_dtdy;
	 dfdt[I1_J1_K1 + MYVECIND] += corrWaveY;
	 dfdt[I1_J0_K1 + MYVECIND] -= corrWaveY;
         #endif

	 if (Vz >= ZERO) {
	    const Vec4 transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*abs(allVx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J0_K2 + MYVECIND] -= transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYVECIND] += transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYVECIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K1 + MYVECIND] -= transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*RY_limiter_yp1xcc_dtdy*dt_per_dz;
	    dfdt[I1_J1_K2 + MYVECIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K1 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K2 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYVECIND] += transCorrWaveYZ;
	    const Vec4 transCorrWaveZY = HALF*Vy*Vz*(ONE - dt_per_dz*Vz)*RZ_limiter_zm1zm2_dtdz*dt_per_dy;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveZY;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J0_K0 + MYVECIND] += transCorrWaveZY;
	    #endif
	 } else {
	    const Vec4 transIncrWave = HALF*Vy*Vz*dt_per_dy*dt_per_dz - SIXTH*abs(allVx)*Vy*Vz*dt_per_dx*dt_per_dy*dt_per_dz;
	    dfdt[I1_J0_K0 + MYVECIND] += transIncrWave*RY;
	    dfdt[I1_J0_K1 + MYVECIND] -= transIncrWave*RY;
	    dfdt[I1_J0_K0 + MYVECIND] += transIncrWave*RZ;
	    dfdt[I1_J1_K0 + MYVECIND] -= transIncrWave*RZ;
	    
	    #if MOVER_VLASOV_ORDER > 1
	    const Vec4 transCorrWaveYZ = -HALF*Vy*Vz*(ONE + dt_per_dy*Vy)*RY_limiter_yp1xcc_dtdy*dt_per_dz;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveYZ;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveYZ;
	    dfdt[I1_J0_K0 + MYVECIND] += transCorrWaveYZ;
	    const Vec4 transCorrWaveZY = -HALF*Vy*Vz*(ONE + dt_per_dz*Vz)*RZ_limiter_zp1xcc_dtdz*dt_per_dy;
	    dfdt[I1_J1_K1 + MYVECIND] += transCorrWaveZY;
	    dfdt[I1_J0_K1 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J1_K0 + MYVECIND] -= transCorrWaveZY;
	    dfdt[I1_J0_K0 + MYVECIND] += transCorrWaveZY;
	    #endif
	 }
      }
      
      if (Vz >= ZERO) {
	 const Vec4 incrWaveZ = Vz*zm1*dt_per_dz;
	 dfdt[I1_J1_K1 + MYVECIND] += incrWaveZ;
	 dfdt[I1_J1_K0 + MYVECIND] -= incrWaveZ;
	 #if MOVER_VLASOV_ORDER > 1
	 const Vec4 corrWaveZ = HALF*Vz*(ONE - dt_per_dz*Vz)*RZ_limiter_zm1zm2_dtdz;
	 dfdt[I1_J1_K1 + MYVECIND] += corrWaveZ;
	 dfdt[I1_J1_K0 + MYVECIND] -= corrWaveZ;
         #endif
      } else {
	 const Vec4 incrWaveZ = Vz*xcc*dt_per_dz;
	 dfdt[I1_J1_K1 + MYVECIND] += incrWaveZ;
	 dfdt[I1_J1_K0 + MYVECIND] -= incrWaveZ;
	 #if MOVER_VLASOV_ORDER > 1
	 const Vec4 corrWaveZ = -HALF*Vz*(ONE + dt_per_dz*Vz)*RZ_limiter_zp1xcc_dtdz;
	 dfdt[I1_J1_K1 + MYVECIND] += corrWaveZ;
	 dfdt[I1_J1_K0 + MYVECIND] -= corrWaveZ;
         #endif
      }      
   }
   

   // Accumulate calculated df/dt values from temporary buffer to 
   // main memory. If multithreading is used over blocks, these updates do
   //not need to be atomistic as long as all threads work on the same
   //cell

   //FIXME, this extra copying is perhaps not needed, why not use block->fx's directly?  check
   const UINT procBoundaryFlags = cell->procBoundaryFlag;
   for (uint nbr=0; nbr<27; ++nbr) {
      // If the neighbour does not exist, do not copy data:     
      if (((procBoundaryFlags >> nbr) & 1) == 0) continue;
      spatial_cell::Velocity_Block* nbrblock = mpiGrid[cell->neighbors[nbr]]->at(blockId);
       //if neighbour is null_block, do not copy fluxes to it
      if(mpiGrid[cell->neighbors[nbr]]->is_null_block(nbrblock)) continue;
      
      //the flux array is not zeroed, rather a marker has been put to mark if it is non-initialized
      //check for that, and initialize by using = instead of += if that is the case

      if( nbrblock->fx[0] == std::numeric_limits<REAL>::max()) {
         for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
            const UINT MYVECIND =k*WID + j;
            const UINT MYIND = MYVECIND*WID;
            dfdt[nbr*WID2 + MYVECIND].store((nbrblock->fx)+MYIND);  
         }
      }
      else{
         for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
            const UINT MYVECIND =k*WID + j;
            const UINT MYIND = MYVECIND*WID;
            Vec4 fx;
            fx.load((nbrblock->fx)+MYIND);
            fx+=dfdt[nbr*WID2 + MYVECIND];
            fx.store((nbrblock->fx)+MYIND);
         }
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
template<typename REAL,typename UINT> void cpu_propagateSpatWithMoments(
	REAL * const nbrFluxes,
	spatial_cell::SpatialCell *cell,
	const UINT blockId,
	const uint block_i
) {
   spatial_cell::Velocity_Block* block=cell->at(blockId);

   // Propagate distribution function:
   if (nbrFluxes == NULL) {
      // No remote neighbour contributions to df/dt
      for (UINT i=0; i<WID3; ++i) {
         block->data[i]+=block->fx[i];
      }
   } else {
      // Cell has remote neighbour contributions to df/dt
      for (UINT i=0; i<WID3; ++i) {
         block->data[i]+=block->fx[i]+nbrFluxes[block_i*WID3 + i];
      }
   }
   
   // Calculate velocity moments:
   cpu_blockVelocityMoments(block->data,block->parameters,cell->parameters);

}

template<typename UINT> void cpu_calcVelocityMoments(spatial_cell::SpatialCell *cell,const UINT blockId){
   spatial_cell::Velocity_Block* block=cell->at(blockId); //returns a reference to block            
   // Calculate velocity moments:
   cpu_blockVelocityMoments(block->data,block->parameters,cell->parameters);

}

#endif
