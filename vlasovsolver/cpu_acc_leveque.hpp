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

#ifndef CPU_ACC_H
#define CPU_ACC_H

#include "definitions.h"
#include "common.h"
#include "spatial_cell.hpp"
#include "project.h"
#include "leveque_common.h"

using namespace spatial_cell;
// Constant for switch statement in solver to decide which of the eight 
// possibilities should be calculated. 
const uint AXP_AYP_AZP = 0; // Ax > 0, Ay > 0, Az > 0
const uint AXP_AYP_AZN = 1; //                 Az < 0
const uint AXP_AYN_AZP = 2; //         Ay < 0, Az > 0
const uint AXP_AYN_AZN = 3; //                 Az < 0
const uint AXN_AYP_AZP = 4; // Ax < 0, Ay > 0, Az > 0
const uint AXN_AYP_AZN = 5; //                 Az < 0
const uint AXN_AYN_AZP = 6; //         Ay < 0, Az > 0
const uint AXN_AYN_AZN = 7; //                 Az < 0

template<typename T> inline T accIndex(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}
template<typename T> inline T fullInd(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> inline T isBoundary(const T& STATE,const T& BND) {return STATE & BND;}
template<typename T> inline T findex(const T& i,const T& j,const T& k) {return k*36+j*6+i;}

void accumulateChanges(Velocity_Block *block,const Real* const dF){
   Real* nbrFlux;
   //typedef Parameters P;
   
   // NOTE: velocity block can have up to 26 neighbours, and we need to copy changes 
   // to each existing neighbour here.
   // NOTE2: dF is a (6,6,6) sized block, nbrFlux is (4,4,4).
   // NOTE3: This function copies values correctly
   const unsigned int MIN = 0;
   const unsigned int MAX = 5;
   const unsigned int ACCMIN = 0;
   const unsigned int ACCMAX = 3;
//   typedef velocity_block_neighbors NbrsVel;
   
   // ***** NEIGHBOURS TO NEGATIVE VZ DIRECTION *****           
   // Neighbour (iv-1,jv-1,kv-1):
   
   if (block->neighbors[velocity_neighbor::XM1_YM1_ZM1] != NULL){
      nbrFlux = block->neighbors[velocity_neighbor::XM1_YM1_ZM1]->fx;
      nbrFlux[accIndex(ACCMAX,ACCMAX,ACCMAX)] += dF[findex(MIN,MIN,MIN)];
   }
   // Neighbour (iv  ,jv-1,kv-1):
   if (block->neighbors[velocity_neighbor::XCC_YM1_ZM1] != NULL) {
      nbrFlux = block->neighbors[velocity_neighbor::XCC_YM1_ZM1]->fx;
      for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,ACCMAX)] += dF[findex(i+1,MIN,MIN)];
   }
   // Neighbour (iv+1,jv-1,kv-1):
   if (block->neighbors[velocity_neighbor::XP1_YM1_ZM1] != NULL) {
      nbrFlux = block->neighbors[velocity_neighbor::XP1_YM1_ZM1]->fx;
      nbrFlux[accIndex(ACCMIN,ACCMAX,ACCMAX)] += dF[findex(MAX,MIN,MIN)];
   }
   // Neighbour (iv-1,jv  ,kv-1):
   if (block->neighbors[velocity_neighbor::XM1_YCC_ZM1] != NULL) {
      nbrFlux = block->neighbors[velocity_neighbor::XM1_YCC_ZM1]->fx;
      for (unsigned int j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,ACCMAX)] += dF[findex(MIN,j+1,MIN)];
   }
   // Neighbour (iv  ,jv  ,kv-1):
   if (block->neighbors[velocity_neighbor::XCC_YCC_ZM1] != NULL) {
      nbrFlux = block->neighbors[velocity_neighbor::XCC_YCC_ZM1]->fx;
      for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,j,ACCMAX)] += dF[findex(i+1,j+1,MIN)];
   }
   // Neighbour (iv+1,jv  ,kv-1):
   if (block->neighbors[velocity_neighbor::XP1_YCC_ZM1] != NULL) {
      nbrFlux  = block->neighbors[velocity_neighbor::XP1_YCC_ZM1]->fx;
      for (unsigned int j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,ACCMAX)] += dF[findex(MAX,j+1,MIN)];
   }
   // Neighbour (iv-1,jv+1,kv-1):
   if (block->neighbors[velocity_neighbor::XM1_YP1_ZM1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YP1_ZM1]->fx;
      nbrFlux[accIndex(ACCMAX,ACCMIN,ACCMAX)] += dF[findex(MIN,MAX,MIN)];
   }
   // Neighbour (iv  ,jv+1,kv-1):
   if (block->neighbors[velocity_neighbor::XCC_YP1_ZM1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XCC_YP1_ZM1]->fx;
      for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,ACCMAX)] += dF[findex(i+1,MAX,MIN)];
   }
   // Neighbour (iv+1,jv+1,kv-1):
   if (block->neighbors[velocity_neighbor::XP1_YP1_ZM1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XP1_YP1_ZM1]->fx;
      nbrFlux[accIndex(ACCMIN,ACCMIN,ACCMAX)] += dF[findex(MAX,MAX,MIN)];
   }
   
   // ***** NEIGHBOURS IN SAME VZ PLANE *****
   // Neighbour (iv-1,jv-1,kv  ):
   if (block->neighbors[velocity_neighbor::XM1_YM1_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YM1_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) nbrFlux[accIndex(ACCMAX,ACCMAX,k)] += dF[findex(MIN,MIN,k+1)];
   }
   // Neighbour (iv  ,jv-1,kv  ):
   if (block->neighbors[velocity_neighbor::XCC_YM1_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XCC_YM1_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,k)] += dF[findex(i+1,MIN,k+1)];
   }
   // Neighbour (iv+1,jv-1,kv  ):
   if (block->neighbors[velocity_neighbor::XP1_YM1_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XP1_YM1_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) nbrFlux[accIndex(ACCMIN,ACCMAX,k)] += dF[findex(MAX,MIN,k+1)];
   }
   // Neighbour (iv-1,jv  ,kv  ):
   if (block->neighbors[velocity_neighbor::XM1_YCC_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YCC_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,k)] += dF[findex(MIN,j+1,k+1)];
   }
   
   // This block (iv  ,jv  ,kv  ):
   for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) {
      block->fx[accIndex(i,j,k)] += dF[findex(i+1,j+1,k+1)];
   }
   // Neighbour (iv+1,jv  ,kv  ):
   if (block->neighbors[velocity_neighbor::XP1_YCC_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XP1_YCC_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,k)] += dF[findex(MAX,j+1,k+1)];
   }
   // Neighbour (iv-1,jv+1,kv  ):
   if (block->neighbors[velocity_neighbor::XM1_YP1_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YP1_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) nbrFlux[accIndex(ACCMAX,ACCMIN,k)] += dF[findex(MIN,MAX,k+1)];
   }
   // Neighbour (iv  ,jv+1,kv  ):
   if (block->neighbors[velocity_neighbor::XCC_YP1_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XCC_YP1_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,k)] += dF[findex(i+1,MAX,k+1)];
   }   
   // Neighbour (iv+1,jv+1,kv  ):
   if (block->neighbors[velocity_neighbor::XP1_YP1_ZCC] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XP1_YP1_ZCC]->fx;
      for (unsigned int k=0; k<WID; ++k) nbrFlux[accIndex(ACCMIN,ACCMIN,k)] += dF[findex(MAX,MAX,k+1)];
   }

   // ***** NEIGHBOURS TO POSITIVE VZ DIRECTION *****
   // Neighbour (iv-1,jv-1,kv+1):
   if (block->neighbors[velocity_neighbor::XM1_YM1_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YM1_ZP1]->fx;
      nbrFlux[accIndex(ACCMAX,ACCMAX,ACCMIN)] += dF[findex(MIN,MIN,MAX)];
   }   
   // Neighbour (iv  ,jv-1,kv+1):
   if (block->neighbors[velocity_neighbor::XCC_YM1_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XCC_YM1_ZP1]->fx;
      for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,ACCMIN)] += dF[findex(i+1,MIN,MAX)];
   }
   // Neighbour (iv+1,jv-1,kv+1):
   if (block->neighbors[velocity_neighbor::XP1_YM1_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XP1_YM1_ZP1]->fx;
      nbrFlux[accIndex(ACCMIN,ACCMAX,ACCMIN)] += dF[findex(MAX,MIN,MAX)];
   }
   // Neighbour (iv-1,jv  ,kv+1):
   if (block->neighbors[velocity_neighbor::XM1_YCC_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YCC_ZP1]->fx;
      for (unsigned int j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,ACCMIN)] += dF[findex(MIN,j+1,MAX)];
   }
   // Neighbour (iv  ,jv  ,kv+1):
   if (block->neighbors[velocity_neighbor::XCC_YCC_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XCC_YCC_ZP1]->fx;
      for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,j,ACCMIN)] += dF[findex(i+1,j+1,MAX)];
   }
   // Neighbour (iv+1,jv  ,kv+1):
   if (block->neighbors[velocity_neighbor::XP1_YCC_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XP1_YCC_ZP1]->fx;
      for (unsigned int j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,ACCMIN)] += dF[findex(MAX,j+1,MAX)];
   }
   // Neighbour (iv-1,jv+1,kv+1):
   if (block->neighbors[velocity_neighbor::XM1_YP1_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XM1_YP1_ZP1]->fx;
      nbrFlux[accIndex(ACCMAX,ACCMIN,ACCMIN)] += dF[findex(MIN,MAX,MAX)];
   }
   // Neighbour (iv  ,jv+1,kv+1):
   if (block->neighbors[velocity_neighbor::XCC_YP1_ZP1] != NULL) {
      nbrFlux =  block->neighbors[velocity_neighbor::XCC_YP1_ZP1]->fx;
      for (unsigned int i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,ACCMIN)] += dF[findex(i+1,MAX,MAX)];
   }
   // Neighbour (iv+1,jv+1,kv+1):
   if (block->neighbors[velocity_neighbor::XP1_YP1_ZP1] != NULL) {
      nbrFlux = block->neighbors[velocity_neighbor::XP1_YP1_ZP1]->fx;
      nbrFlux[accIndex(ACCMIN,ACCMIN,ACCMIN)] += dF[findex(MAX,MAX,MAX)];
   }
}


void fetchAllAverages(Velocity_Block* block,Real* const avgs) {

   Velocity_Block *nbrBlock;
   for (unsigned int i=0; i<8*WID3; ++i) avgs[i] = 0.0;
   
   
   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (block->neighbors[velocity_neighbor::XM1_YCC_ZCC] == NULL) {
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<2; ++i) {
	 avgs[fullInd(i  ,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = block->neighbors[velocity_neighbor::XM1_YCC_ZCC];
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (unsigned int i=0; i<2; ++i) {
            avgs[fullInd(i  ,j+2,k+2)] = nbrBlock->data[accIndex(i+2,j,k)];
	 }
      }
   }
   
   // Copy averages from +x neighbour, or calculate using a boundary function:
   if (block->neighbors[velocity_neighbor::XP1_YCC_ZCC] == NULL) {
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<2; ++i) {
	 avgs[fullInd(i+6,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = block->neighbors[velocity_neighbor::XP1_YCC_ZCC];
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (unsigned int i=0; i<2; ++i) {
	    avgs[fullInd(i+6,j+2,k+2)] = nbrBlock->data[accIndex(i,j,k)];
	 }
      }
   }
   
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (block->neighbors[velocity_neighbor::XCC_YM1_ZCC] == NULL) {
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<2; ++j) for (unsigned int i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j  ,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = block->neighbors[velocity_neighbor::XCC_YM1_ZCC];
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<2; ++j) {
	 #pragma ivdep
	 for (unsigned int i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j  ,k+2)] = nbrBlock->data[accIndex(i,j+2,k)];
	 }
      }
   }
   // Copy averages from +y neighbour, or calculate using a boundary function:
   if (block->neighbors[velocity_neighbor::XCC_YP1_ZCC] == NULL) {
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<2; ++j) for (unsigned int i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+6,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = block->neighbors[velocity_neighbor::XCC_YP1_ZCC];
      for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<2; ++j) {
	 #pragma ivdep
	 for (unsigned int i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j+6,k+2)] = nbrBlock->data[accIndex(i,j,k)];
	 }
      }
   }
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (block->neighbors[velocity_neighbor::XCC_YCC_ZM1] == NULL) {
      for (unsigned int k=0; k<2; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k  )] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = block->neighbors[velocity_neighbor::XCC_YCC_ZM1];
      for (unsigned int k=0; k<2; ++k) for (unsigned int j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (unsigned int i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j+2,k  )] = nbrBlock->data[accIndex(i,j,k+2)];
	 }
      }
   }
   // Copy averages from +z neighbour, or calculate using a boundary function:
   if (block->neighbors[velocity_neighbor::XCC_YCC_ZP1] == NULL) {
      for (unsigned int k=0; k<2; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k+6)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = block->neighbors[velocity_neighbor::XCC_YCC_ZP1];
      for (unsigned int k=0; k<2; ++k) for (unsigned int j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (unsigned int i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j+2,k+6)] =  nbrBlock->data[accIndex(i,j,k)];
	 }
      }
   }

   // Copy volume averages of this block:
   for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) {
      avgs[fullInd(i+2,j+2,k+2)] = block->data[accIndex(i,j,k)];
   }
}

//zero fluxes
void cpu_clearVelFluxes(SpatialCell *cell,const unsigned int& BLOCK) {
   Velocity_Block* block=cell->at(BLOCK);
   for (unsigned int i=0; i<SIZE_FLUXS; ++i)  block->fx[i]= 0.0;
}

void cpu_calcVelFluxes(SpatialCell *cell,const unsigned int& BLOCK,const Real& DT,creal* const accmat) {
   // Creation of temporary calculation block dfdt and avgs + 
   // value fetching and initializations seem to take about
   // ~4% of time used by calcVelFluxes
   
   // Allocate temporary array in which local dF changes are calculated:
   // F is the distribution function, dFdt is its change over timestep DT
   Real dfdt[216];
   for (unsigned int i=0; i<216; ++i) dfdt[i] = 0.0;
   
   Velocity_Block* block=cell->at(BLOCK);
   Real avgs[8*WID3];
   fetchAllAverages(block,avgs);

   const Real dt_per_dvx = DT / block->parameters[BlockParams::DVX];
   const Real dt_per_dvy = DT / block->parameters[BlockParams::DVY];
   const Real dt_per_dvz = DT / block->parameters[BlockParams::DVZ];
   
   Real Ax,Ay,Az;
   for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) {
      Real R,incrWave,transIncrWave,doubleTransIncrWave;
      Real corrWave,transCorrWave,doubleTransCorrWave;
      unsigned int solverFlags;
      
      // ***********************************
      // ***** INTERFACE BETWEEN I,I-1 *****
      // ***********************************
      const Real xcc = avgs[fullInd(i+2,j+2,k+2)];
      const Real xp1 = avgs[fullInd(i+3,j+2,k+2)];
      const Real xm1 = avgs[fullInd(i+1,j+2,k+2)];
      const Real xm2 = avgs[fullInd(i  ,j+2,k+2)];
      calcAccFaceX(Ax,Ay,Az,i,j,k,cell->parameters,block->parameters);

      solverFlags = 0;
      if (Az < ZERO) solverFlags = (solverFlags | (1 << 0));
      if (Ay < ZERO) solverFlags = (solverFlags | (1 << 1));
      if (Ax < ZERO) solverFlags = (solverFlags | (1 << 2));

      switch (solverFlags) {
       case AXP_AYP_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
	 #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;

	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+2)] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvx*Ax)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYP_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;

	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k  )] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 #endif
	 break;
       case AXP_AYN_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;

	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+2)] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE-Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
         #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE-Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
	 #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+2)] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+2,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+2)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] += doubleTransCorrWave;
         #endif
	 break;	 
       case AXN_AYP_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
         dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k  )] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
         #endif
	 break;	 
       case AXN_AYN_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+2)] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k  )] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+1,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
      }
      
      // ***********************************
      // ***** INTERFACE BETWEEN J,J-1 *****
      // ***********************************
      const Real yp1 = avgs[fullInd(i+2,j+3,k+2)];
      const Real ym1 = avgs[fullInd(i+2,j+1,k+2)];
      const Real ym2 = avgs[fullInd(i+2,j  ,k+2)];
      calcAccFaceY(Ax,Ay,Az,i,j,k,cell->parameters,block->parameters);
      
      solverFlags = 0;
      if (Az < ZERO) solverFlags = (solverFlags | (1 << 0));
      if (Ay < ZERO) solverFlags = (solverFlags | (1 << 1));
      if (Ax < ZERO) solverFlags = (solverFlags | (1 << 2));

      switch (solverFlags) {
       case AXP_AYP_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
	 #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+2)] -= TWO*doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] += doubleTransCorrWave;
	 #endif
	 break;
       case AXP_AYP_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k  )] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
	 #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+2)] -= TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j  ,k+2)] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j  ,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] += doubleTransCorrWave;
	 #endif	 
	 break;
       case AXP_AYN_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k  )] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+2)] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] -= TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j+1,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+2)] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k  )] -= TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
      }

      // ***********************************
      // ***** INTERFACE BETWEEN K,K-1 *****
      // ***********************************
      const Real zp1 = avgs[fullInd(i+2,j+2,k+3)];
      const Real zm1 = avgs[fullInd(i+2,j+2,k+1)];
      const Real zm2 = avgs[fullInd(i+2,j+2,k  )];
      calcAccFaceZ(Ax,Ay,Az,i,j,k,cell->parameters,block->parameters);
      
      solverFlags = 0;
      if (Az < ZERO) solverFlags = (solverFlags | (1 << 0));
      if (Ay < ZERO) solverFlags = (solverFlags | (1 << 1));
      if (Ax < ZERO) solverFlags = (solverFlags | (1 << 2));
      
      switch (solverFlags) {
       case AXP_AYP_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
	 #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+2,k+1)] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] += doubleTransCorrWave;
	 #endif
	 break;
       case AXP_AYP_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+2,k  )] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j+2,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+1)] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k  )] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+1)] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k  )] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
         #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k  )] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i  ,j+1,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
      }
   }
   
   // If multithreading is used (OpenMP/pthreads), then some of the updates 
   // in accumulateChanges may need to be atomic. dfdt values are copied to 
   // neighbouring blocks within the same spatial cell.
   // 
   // accumulateChanges seems to take ~3% of time used by calcVelFluxes
   accumulateChanges(block,dfdt);

}

void cpu_propagateVel(SpatialCell *cell,const unsigned int& BLOCK,const Real& DT) {
   Velocity_Block* block=cell->at(BLOCK);
   for (unsigned int i=0; i<WID3; ++i)  block->data[i] += block->fx[i];



   
}
#endif
