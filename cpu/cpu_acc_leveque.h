#ifndef CPU_ACC_H
#define CPU_ACC_H

#include <map>
#include "definitions.h"
#include "common.h"
#include "cell_spatial.h"
#include "cpu_common.h"
#include "project.h"

template<typename T> inline T accIndex(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}
template<typename T> inline T fullInd(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> inline T isBoundary(const T& STATE,const T& BND) {return STATE & BND;}

template<typename T> inline T findex(const T& i,const T& j,const T& k) {return k*36+j*6+i;}

template<typename T> inline T limiter(const T& THETA) {
   //return vanLeer(THETA);
   //return MClimiter(THETA);
   return superbee(THETA);
   //return modbee2(THETA);
}

template<typename T> T velDerivs1(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2) {
   return superbee(xl1,xcc,xr1);
}

template<typename T> T velDerivs2(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2) {
   return convert<T>(0.0);
}

template<typename REAL,typename UINT> void accumulateChanges(const UINT& BLOCK,const REAL* const dF,REAL* const flux,const UINT* const nbrsVel) {
   UINT nbrBlock;
   REAL* nbrFlux;
   const UINT STATE = nbrsVel[NbrsVel::STATE];
   typedef Parameters P;

   // NOTE: velocity block can have up to 26 neighbours, and we need to copy changes 
   // to each existing neighbour here.
   // NOTE2: dF is a (6,6,6) sized block, nbrFlux is (4,4,4).
   // NOTE3: This function copies values correctly
   const UINT MIN = 0;
   const UINT MAX = 5;
   const UINT ACCMIN = 0;
   const UINT ACCMAX = 3;
   
   // Accumulate changes to this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) 
     flux[BLOCK*SIZE_FLUXS+accIndex(i,j,k)] += dF[findex(i+1,j+1,k+1)];
   
   // Accumulate changes to (iv-1,jv,kv) neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VXNEG];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,k)] += dF[findex(MIN,j+1,k+1)];
      // Accumulate changes to (iv-1,jv+-1,kv) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMAX,ACCMAX,k)] += dF[findex(MIN,MIN,k+1)];
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMAX,ACCMIN,k)] += dF[findex(MIN,MAX,k+1)];
      }
   }
   // Accumulate changes to (iv+1,jv,kv) neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VXPOS];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,k)] += dF[findex(MAX,j+1,k+1)];
      // Accumulate changes to (iv+1,jv+-1,kv) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMIN,ACCMAX,k)] += dF[findex(MAX,MIN,k+1)];
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMIN,ACCMIN,k)] += dF[findex(MAX,MAX,k+1)];
      }
   }
   // Accumulate changes to -vy neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VYNEG];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,k)] += dF[findex(i+1,MIN,k+1)];
      // Accumulate changes to (iv,jv-1,kv+-1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,ACCMAX)] += dF[findex(i+1,MIN,MIN)];
      }
      if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,ACCMIN)] += dF[findex(i+1,MIN,MAX)];
      }
   }
   // Accumulate changes to +vy neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VYPOS];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,k)] += dF[findex(i+1,MAX,k+1)];
      // Accumulate changes to (iv,jv+1,kv+-1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini; // temp solution
      	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,ACCMAX)] += dF[findex(i+1,MAX,MIN)];
      }
      if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,ACCMIN)] += dF[findex(i+1,MAX,MAX)];
      }
   }
   // Accumulate changes to -vz neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VZNEG];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,j,ACCMAX)] += dF[findex(i+1,j+1,MIN)];
      // Accumulate changes to (iv+-1,jv,kv-1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,ACCMAX)] += dF[findex(MIN,j+1,MIN)];
      }
      if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,ACCMAX)] += dF[findex(MAX,j+1,MIN)];
      }
   }
   // Accumulate changes to +vz neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VZPOS];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,j,ACCMIN)] += dF[findex(i+1,j+1,MAX)];
      // Accumulate changes to (iv+-1,jv,kv+1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,ACCMIN)] += dF[findex(MIN,j+1,MAX)];
      }
      if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,ACCMIN)] += dF[findex(MAX,j+1,MAX)];
      }
   }
   
   // Accumulate changes to 8 corner neighbours:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv-1,jv-1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMAX,ACCMAX)] += dF[findex(MIN,MIN,MIN)];
	 }
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv-1,jv-1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMAX,ACCMIN)] += dF[findex(MIN,MIN,MAX)];
	 }
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv-1,jv+1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMIN,ACCMAX)] += dF[findex(MIN,MAX,MIN)];
	 } 
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv-1,jv+1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMIN,ACCMIN)] += dF[findex(MIN,MAX,MAX)];
	 }
      }
   }
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv+1,jv-1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMAX,ACCMAX)] += dF[findex(MAX,MIN,MIN)];
	 }
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv+1,jv-1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMAX,ACCMIN)] += dF[findex(MAX,MIN,MAX)];
	 }
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) { 
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv+1,jv+1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMIN,ACCMAX)] += dF[findex(MAX,MAX,MIN)];
	 }
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv+1,jv+1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMIN,ACCMIN)] += dF[findex(MAX,MAX,MAX)];
	 }
      }
   }
}

template<typename REAL,typename UINT> void fetchAllAverages(const UINT& BLOCK,REAL* const avgs,const REAL* const cpu_avgs,const UINT* const nbrsVel) {
   UINT nbrBlock;
   for (UINT i=0; i<8*WID3; ++i) avgs[i] = 0.0;
   const UINT STATE = nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE];
   
   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<2; ++i) {
	 avgs[fullInd(i  ,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK; 
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<2; ++i) {
	 //avgs[fullInd(i  ,j+2,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i+2,j,k)];
	 avgs[fullInd(i  ,j+2,k+2)] = tmp[accIndex(i+2,j,k)];
	 }
      }
   }
   // Copy averages from +x neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<2; ++i) {
	 avgs[fullInd(i+6,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<2; ++i) {
	    //avgs[fullInd(i+6,j+2,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i  ,j,k)];
	    avgs[fullInd(i+6,j+2,k+2)] = tmp[accIndex(i  ,j,k)];
	 }
      }
   }
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VY_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j  ,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j  ,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j+2,k)];
	    avgs[fullInd(i+2,j  ,k+2)] = tmp[accIndex(i,j+2,k)];
	 }
      }
   }
   // Copy averages from +y neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VY_POS_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+6,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j+6,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j  ,k)];
	    avgs[fullInd(i+2,j+6,k+2)] = tmp[accIndex(i,j  ,k)];
	 }
      }
   }
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) > 0) {
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k  )] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j+2,k  )] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,k+2)];
	    avgs[fullInd(i+2,j+2,k  )] = tmp[accIndex(i,j,k+2)];
	 }
      }
   }
   // Copy averages from +z neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VZ_POS_BND) > 0) {
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k+6)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j+2,k+6)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,k)];
	    avgs[fullInd(i+2,j+2,k+6)] = tmp[accIndex(i,j,k)];
	 }
      }
   }
   // Copy volume averages of this block:
   creal* const tmp = cpu_avgs + BLOCK*SIZE_VELBLOCK;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
      #pragma ivdep
      for (UINT i=0; i<WID; ++i) {
	 //avgs[fullInd(i+2,j+2,k+2)] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
	 avgs[fullInd(i+2,j+2,k+2)] = tmp[accIndex(i,j,k)];
      }
   }
}
   
template<typename REAL,typename UINT> void fetchAveragesX(const UINT& BLOCK,REAL* const avgs,const REAL* const cpu_avgs,const UINT* const nbrsVel) {
   // The size of array avgs is (5,4,4).
   const UINT XS=5;
   const UINT YS=4;
   const UINT YSXS = YS*XS;
   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VX_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 avgs[k*YSXS+j*XS] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 avgs[k*YSXS+j*XS] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(convert<UINT>(3),j,k)];
      }
   }
   // Copy volume averages of this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[k*YSXS+j*XS+(i+1)] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT> void fetchAveragesY(const UINT& BLOCK,REAL* const avgs,const REAL* const cpu_avgs,const UINT* const nbrsVel) {
   // The size of array avgs (4,5,4).
   cuint XS=4;
   cuint YS=5;
   cuint YSXS = YS*XS;
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VY_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) {
	 avgs[k*YSXS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) {
	 avgs[k*YSXS+i] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,convert<UINT>(3),k)];
      }
   }
   // Copy volume averages of this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[k*YSXS+(j+1)*XS+i] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
   }   
}
   
template<typename REAL,typename UINT> void fetchAveragesZ(const UINT& BLOCK,REAL* const avgs,const REAL* const cpu_avgs,const UINT* const nbrsVel) {   
   // The size of array avgs (4,4,5).
   const UINT XS=4;
   const UINT YS=4;
   const UINT YSXS=YS*XS;
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VZ_NEG_BND) > 0) {
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[j*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[j*XS+i] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,convert<UINT>(3))];
      }
   }
   // Copy volume averages of this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[(k+1)*YSXS+j*XS+i] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT> void fetchDerivsX(const UINT& BLOCK,REAL* const d1x,const REAL* const cpu_d1x,const UINT* const nbrsVel) {
   // The size of array avgs (5,4,4).
   const UINT XS=5;
   const UINT YS=4;
   const UINT YSXS=YS*XS;
   // Copy derivatives from -x neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VX_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 d1x[k*YSXS+j*XS] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 d1x[k*YSXS+j*XS] = cpu_d1x[nbrBlock*SIZE_DERIV + accIndex(convert<UINT>(3),j,k)];
      }
   }
   // Copy derivatives for this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      d1x[k*YSXS+j*XS+(i+1)] = cpu_d1x[BLOCK*SIZE_DERIV + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT> void fetchDerivsY(const UINT& BLOCK,REAL* const d1y,const REAL* const cpu_d1y,const UINT* const nbrsVel) {
   // The size of array avgs (4,5,4).
   const UINT XS=4;
   const UINT YS=5;
   const UINT YSXS=YS*XS;
   // Copy derivatives from -y neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VY_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) {
	 d1y[k*YSXS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) {
	 d1y[k*YSXS+i] = cpu_d1y[nbrBlock*SIZE_DERIV + accIndex(i,convert<UINT>(3),k)];
      }
   }
   // Copy derivatives for this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      d1y[k*YSXS+(j+1)*XS+i] = cpu_d1y[BLOCK*SIZE_DERIV + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT> void fetchDerivsZ(const UINT& BLOCK,REAL* const d1z,const REAL* const cpu_d1z,const UINT* const nbrsVel) {
   // The size of array avgs (4,4,5)
   const UINT XS=4;
   const UINT YS=4;
   const UINT YSXS=YS*XS;
   // Copy derivatives from -z neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VZ_NEG_BND) > 0) {
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 d1z[j*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 d1z[j*XS+i] = cpu_d1z[nbrBlock*SIZE_DERIV + accIndex(i,j,convert<UINT>(3))];
      }
   }
   // Copy derivatives for this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      d1z[(k+1)*YSXS+j*XS+i] = cpu_d1z[BLOCK*SIZE_DERIV + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT> void fetchFluxesX(const UINT& BLOCK,REAL* const flux,const REAL* const cpu_fx,const UINT* const nbrsVel) {
   // The size of array fx is (5,4,4):
   const UINT XS = 5;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   // Fetch fluxes from +x neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VX_POS_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 flux[k*YSXS+j*XS+4] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXPOS];
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 flux[k*YSXS+j*XS+4] = cpu_fx[nbrBlock*SIZE_FLUXS + accIndex(convert<UINT>(0),j,k)];
      }
   }
   // Fetch the fluxes of this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      flux[k*YSXS+j*XS+i] = cpu_fx[BLOCK*SIZE_FLUXS + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT> void fetchFluxesY(const UINT& BLOCK,REAL* const flux,const REAL* const cpu_fy,const UINT* const nbrsVel) {
   // The size of array fy is (4,5,4):
   const UINT XS = 4;
   const UINT YS = 5;
   const UINT YSXS = YS*XS;
   // Fetch fluxes from +y neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VY_POS_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) {
	 flux[k*YSXS+4*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYPOS];
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) {
	 flux[k*YSXS+4*XS+i] = cpu_fy[nbrBlock*SIZE_FLUXS + accIndex(i,convert<uint>(0),k)];
      }
   }
   // Fetch the fluxes of this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      flux[k*YSXS+j*XS+i] = cpu_fy[BLOCK*SIZE_FLUXS + accIndex(i,j,k)];
   }
}
   
template<typename REAL,typename UINT> void fetchFluxesZ(const UINT& BLOCK,REAL* const flux,const REAL* const cpu_fz,const UINT* const nbrsVel) {
   // The size of array fz is (4,4,5):
   const UINT XS = 4;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   // Fetch fluxes from +z neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE],NbrsVel::VZ_POS_BND) > 0) {
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 flux[4*YSXS+j*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      const UINT nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZPOS];
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 flux[4*YSXS+j*XS+i] = cpu_fz[nbrBlock*SIZE_FLUXS + accIndex(i,j,convert<UINT>(0))];
      }
   }
   // Fetch the fluxes of this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      flux[k*YSXS+j*XS+i] = cpu_fz[BLOCK*SIZE_FLUXS + accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_clearVelFluxes(CELL& cell,const UINT& BLOCK) {
   for (UINT i=0; i<SIZE_FLUXS; ++i) cell.cpu_fx[BLOCK*SIZE_FLUXS + i] = 0.0;
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxes(CELL& cell,const UINT& BLOCK,const REAL& DT,creal* const accmat) {
   
   const REAL EPSILON = 1.0e-15;

   // Allocate temporary array in which local dF changes are calculated:
   // F is the distribution function, dFdt is its change over timestep DT
   Real dFdt[216];
   for (UINT i=0; i<216; ++i) dFdt[i] = 0.0;
   
   const REAL* const cellParams = cell.cpu_cellParams;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS; 
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   
   REAL AX,AY,AZ,VX,VY,VZ;
   REAL dF,theta;
   REAL AX_NEG,AX_POS,AY_NEG,AY_POS,AZ_NEG,AZ_POS;
   REAL INCR_WAVE,CORR_WAVE;
   REAL theta_up,theta_lo;
   const REAL DVX = blockParams[BlockParams::DVX];
   const REAL DVY = blockParams[BlockParams::DVY];
   const REAL DVZ = blockParams[BlockParams::DVZ];
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      // ***********************************
      // ***** INTERFACE BETWEEN I,I-1 *****
      // ***********************************
      
      // Calculate acceleration at face (i-1/2,j,k):
      VX = blockParams[BlockParams::VXCRD];
      VY = blockParams[BlockParams::VYCRD] + (j+0.5)*DVY;
      VZ = blockParams[BlockParams::VZCRD] + (k+0.5)*DVZ;
      AX = Parameters::q_per_m * (cellParams[CellParams::EX] + VY*cellParams[CellParams::BZ] - VZ*cellParams[CellParams::BY]);
      AY = Parameters::q_per_m * (cellParams[CellParams::EY] + VZ*cellParams[CellParams::BX] - VX*cellParams[CellParams::BZ]);
      AZ = Parameters::q_per_m * (cellParams[CellParams::EZ] + VX*cellParams[CellParams::BY] - VY*cellParams[CellParams::BX]);
      
      AX_NEG = std::min(convert<REAL>(0.0),AX);
      AX_POS = std::max(convert<REAL>(0.0),AX);
      AY_NEG = std::min(convert<REAL>(0.0),AY);
      AY_POS = std::max(convert<REAL>(0.0),AY);
      AZ_NEG = std::min(convert<REAL>(0.0),AZ);
      AZ_POS = std::max(convert<REAL>(0.0),AZ);

      // Calculate slope-limited x-derivative of F:
      if (AX > 0.0) theta_up = avgs[fullInd(i+1,j+2,k+2)] - avgs[fullInd(i  ,j+2,k+2)];
      else theta_up = avgs[fullInd(i+3,j+2,k+2)] - avgs[fullInd(i+2,j+2,k+2)];
      theta_lo = avgs[fullInd(i+2,j+2,k+2)] - avgs[fullInd(i+1,j+2,k+2)] + EPSILON;
      theta = limiter(theta_up/theta_lo);
      
      // Donor cell upwind method, Fx updates:
      dF = avgs[fullInd(i+2,j+2,k+2)] - avgs[fullInd(i+1,j+2,k+2)]; // F(i,j,k) - F(i-1,j,k), jump in F
      INCR_WAVE = AX_POS*avgs[fullInd(i+1,j+2,k+2)] + AX_NEG*avgs[fullInd(i+2,j+2,k+2)];// + 0.5*fabs(AX)*(1.0-fabs(AX)*DT/DVX)*dF*theta;
      CORR_WAVE = 0.5*fabs(AX)*(1.0-fabs(AX)*DT/DVX)*dF*theta;

      cell.cpu_d1x[BLOCK*SIZE_DERIV+accIndex(i,j,k)] = theta;
      cell.cpu_d2x[BLOCK*SIZE_DERIV+accIndex(i,j,k)] = dF*theta;

      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE+CORR_WAVE)*DT/DVX;               // Positive change to i
      dFdt[findex(i  ,j+1,k+1)] -= (INCR_WAVE+CORR_WAVE)*DT/DVX;               // Negative change to i-1

      INCR_WAVE = -0.5*DT/DVX*dF;
      CORR_WAVE *= DT/DVX;

      // Transverse propagation, Fy updates:
      dFdt[findex(i+1,j+2,k+1)] += (INCR_WAVE*AX_POS*AY_POS + CORR_WAVE*AY_POS)*DT/DVY; // Fy: i,j+1
      dFdt[findex(i+1,j+1,k+1)] -= (INCR_WAVE*AX_POS*AY_POS + CORR_WAVE*AY_POS)*DT/DVY; 
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE*AX_POS*AY_NEG + CORR_WAVE*AY_NEG)*DT/DVY; // Fy: i,j
      dFdt[findex(i+1,j  ,k+1)] -= (INCR_WAVE*AX_POS*AY_NEG + CORR_WAVE*AY_NEG)*DT/DVY;
      dFdt[findex(i  ,j+2,k+1)] += (INCR_WAVE*AX_NEG*AY_POS - CORR_WAVE*AY_POS)*DT/DVY; // Fy: i-1,j+1
      dFdt[findex(i  ,j+1,k+1)] -= (INCR_WAVE*AX_NEG*AY_POS - CORR_WAVE*AY_POS)*DT/DVY;
      dFdt[findex(i  ,j+1,k+1)] += (INCR_WAVE*AX_NEG*AY_NEG - CORR_WAVE*AY_NEG)*DT/DVY; // Fy: i,j
      dFdt[findex(i  ,j  ,k+1)] -= (INCR_WAVE*AX_NEG*AY_NEG - CORR_WAVE*AY_NEG)*DT/DVY;

      // Transverse propagation, Fz updates:
      dFdt[findex(i+1,j+1,k+2)] += (INCR_WAVE*AX_POS*AZ_POS + CORR_WAVE*AZ_POS)*DT/DVZ; // Ax > 0, Az > 0
      dFdt[findex(i+1,j+1,k+1)] -= (INCR_WAVE*AX_POS*AZ_POS + CORR_WAVE*AZ_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE*AX_POS*AZ_NEG + CORR_WAVE*AZ_NEG)*DT/DVZ; // Ax > 0, Az < 0
      dFdt[findex(i+1,j+1,k  )] -= (INCR_WAVE*AX_POS*AZ_NEG + CORR_WAVE*AZ_NEG)*DT/DVZ;
      dFdt[findex(i  ,j+1,k+2)] += (INCR_WAVE*AX_NEG*AZ_POS - CORR_WAVE*AZ_POS)*DT/DVZ; // Ax < 0, Az > 0
      dFdt[findex(i  ,j+1,k+1)] -= (INCR_WAVE*AX_NEG*AZ_POS - CORR_WAVE*AZ_POS)*DT/DVZ;
      dFdt[findex(i  ,j+1,k+1)] += (INCR_WAVE*AX_NEG*AZ_NEG - CORR_WAVE*AZ_NEG)*DT/DVZ; // Ax < 0, Az < 0
      dFdt[findex(i  ,j+1,k  )] -= (INCR_WAVE*AX_NEG*AZ_NEG - CORR_WAVE*AZ_NEG)*DT/DVZ;

      INCR_WAVE = DT*DT/DVX/DVZ*dF/6.0;

      // Double transverse propagation, Fy updates:
      dFdt[findex(i+1,j+2,k+1)] += INCR_WAVE*AX_POS*AY_POS*fabs(AZ)*DT/DVY; // Ax > 0, Ay > 0
      dFdt[findex(i+1,j+1,k+1)] += INCR_WAVE*AX_POS*AY_NEG*fabs(AZ)*DT/DVY; //         Ay < 0
      dFdt[findex(i  ,j+2,k+1)] += INCR_WAVE*AX_NEG*AY_POS*fabs(AZ)*DT/DVY; // Ax < 0, Ay > 0
      dFdt[findex(i  ,j+1,k+1)] += INCR_WAVE*AX_NEG*AY_NEG*fabs(AZ)*DT/DVY; //         Ay < 0
      
      dFdt[findex(i+1,j+2,k+2)] -= INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVY; // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i+1,j+1,k+2)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVY; //         Ay < 0
      dFdt[findex(i+1,j+2,k  )] += INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVY; //         Ay > 0, Az < 0
      dFdt[findex(i+1,j+1,k  )] += INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVY; //         Ay < 0
      dFdt[findex(i  ,j+2,k+2)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVY; // Ax < 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+1,k+2)] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVY; //         Ay < 0
      dFdt[findex(i  ,j+2,k  )] += INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVY; //         Ay > 0, Az < 0
      dFdt[findex(i  ,j+1,k  )] += INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVY; //         Ay < 0
      
      dFdt[findex(i+1,j+1,k+1)] -= INCR_WAVE*AX_POS*AY_POS*fabs(AZ)*DT/DVY; // Ax > 0, Ay > 0
      dFdt[findex(i+1,j  ,k+1)] -= INCR_WAVE*AX_POS*AY_NEG*fabs(AZ)*DT/DVY; //         Ay < 0
      dFdt[findex(i  ,j+1,k+1)] -= INCR_WAVE*AX_NEG*AY_POS*fabs(AZ)*DT/DVY; // Ax < 0, Ay > 0
      dFdt[findex(i  ,j  ,k+1)] -= INCR_WAVE*AX_NEG*AY_NEG*fabs(AZ)*DT/DVY; //         Ay < 0
      
      dFdt[findex(i+1,j+1,k+2)] += INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVY; // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i+1,j  ,k+2)] += INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVY; //         Ay < 0
      dFdt[findex(i+1,j+1,k  )] -= INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVY; //         Ay > 0, Az < 0
      dFdt[findex(i+1,j  ,k  )] -= INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVY; //         Ay < 0
      dFdt[findex(i  ,j+1,k+2)] += INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVY; // Ax < 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j  ,k+2)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVY; //         Ay < 0
      dFdt[findex(i  ,j+1,k  )] -= INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVY; //         Ay > 0, Az < 0
      dFdt[findex(i  ,j  ,k  )] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVY; //         Ay < 0
      
      INCR_WAVE = DT*DT/DVX/DVY*dF/6.0;
      
      // Double transverse propagation, Fz updates:
      dFdt[findex(i+1,j+1,k+2)] += INCR_WAVE*AX_POS*fabs(AY)*AZ_POS*DT/DVZ; // Ax > 0, Az > 0
      dFdt[findex(i+1,j+1,k+1)] += INCR_WAVE*AX_POS*fabs(AY)*AZ_NEG*DT/DVZ;
      dFdt[findex(i  ,j+1,k+2)] += INCR_WAVE*AX_NEG*fabs(AY)*AZ_POS*DT/DVZ; // Ax < 0
      dFdt[findex(i  ,j+1,k+1)] += INCR_WAVE*AX_NEG*fabs(AY)*AZ_NEG*DT/DVZ;
      
      dFdt[findex(i+1,j+2,k+2)] -= INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVZ; // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i+1,j+2,k+1)] -= INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVZ; //                 Az < 0
      dFdt[findex(i+1,j  ,k+2)] += INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVZ; //         Ay < 0, Az > 0
      dFdt[findex(i+1,j  ,k+1)] += INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVZ; //                 Az < 0
      dFdt[findex(i  ,j+2,k+2)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVZ; // Ax < 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+2,k+1)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVZ; //                 Az < 0
      dFdt[findex(i  ,j  ,k+2)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVZ; //         Ay < 0, Az > 0
      dFdt[findex(i  ,j  ,k+1)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVZ; //                 Az < 0
      
      dFdt[findex(i+1,j+1,k+1)] -= INCR_WAVE*AX_POS*fabs(AY)*AZ_POS*DT/DVZ; // Ax > 0, Az > 0
      dFdt[findex(i+1,j+1,k  )] -= INCR_WAVE*AX_POS*fabs(AY)*AZ_NEG*DT/DVZ;
      dFdt[findex(i  ,j+1,k+1)] -= INCR_WAVE*AX_NEG*fabs(AY)*AZ_POS*DT/DVZ; // Ax < 0
      dFdt[findex(i  ,j+1,k  )] -= INCR_WAVE*AX_NEG*fabs(AY)*AZ_NEG*DT/DVZ;
      
      dFdt[findex(i+1,j+2,k+1)] += INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVZ; // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i+1,j+2,k  )] += INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVZ; //                 Az < 0
      dFdt[findex(i+1,j  ,k+1)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVZ; //         Ay < 0, Az > 0
      dFdt[findex(i+1,j  ,k  )] -= INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVZ; //                 Az < 0
      dFdt[findex(i  ,j+2,k+1)] += INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVZ; // Ax < 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+2,k  )] += INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVZ; //                 Az < 0
      dFdt[findex(i  ,j  ,k+1)] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVZ; //         Ay < 0, Az > 0
      dFdt[findex(i  ,j  ,k  )] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVZ; //                 Az < 0
   
      // ***********************************
      // ***** INTERFACE BETWEEN J,J-1 *****
      // ***********************************

      // Calculate acceleration at face (i,j-1/2,k):
      VX = blockParams[BlockParams::VXCRD] + (i+0.5)*DVX;
      VY = blockParams[BlockParams::VYCRD];
      VZ = blockParams[BlockParams::VZCRD] + (k+0.5)*DVZ;
      AX = Parameters::q_per_m * (VY*cellParams[CellParams::BZ] - VZ*cellParams[CellParams::BY]);
      AY = Parameters::q_per_m * (VZ*cellParams[CellParams::BX] - VX*cellParams[CellParams::BZ]);
      AZ = Parameters::q_per_m * (VX*cellParams[CellParams::BY] - VY*cellParams[CellParams::BX]);

      AX_NEG = std::min(convert<REAL>(0.0),AX);
      AX_POS = std::max(convert<REAL>(0.0),AX);
      AY_NEG = std::min(convert<REAL>(0.0),AY);
      AY_POS = std::max(convert<REAL>(0.0),AY);
      AZ_NEG = std::min(convert<REAL>(0.0),AZ);
      AZ_POS = std::max(convert<REAL>(0.0),AZ);

      // Calculate slope-limited y-derivative of F:
      if (AY > 0.0) theta_up = avgs[fullInd(i+2,j+1,k+2)] - avgs[fullInd(i+2,j  ,k+2)];
      else theta_up = avgs[fullInd(i+2,j+3,k+2)] - avgs[fullInd(i+2,j+2,k+2)];
      theta_lo = avgs[fullInd(i+2,j+2,k+2)] - avgs[fullInd(i+2,j+1,k+2)] + EPSILON;
      theta = limiter(theta_up/theta_lo);
 
      // Donor cell upwind method, Fy updates:
      dF = avgs[fullInd(i+2,j+2,k+2)] - avgs[fullInd(i+2,j+1,k+2)]; // F(i,j,k) - F(i-1,j,k), jump in F
      INCR_WAVE = AY_POS*avgs[fullInd(i+2,j+1,k+2)] + AY_NEG*avgs[fullInd(i+2,j+2,k+2)];// + 0.5*fabs(AY)*(1.0-fabs(AY)*DT/DVY)*dF*theta;
      CORR_WAVE = 0.5*fabs(AY)*(1.0-fabs(AY)*DT/DVY)*dF*theta;
      cell.cpu_d1y[BLOCK*SIZE_DERIV+accIndex(i,j,k)] = theta;
      cell.cpu_d2y[BLOCK*SIZE_DERIV+accIndex(i,j,k)] = dF*theta;
      //CORR_WAVE = 0.0;
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE+CORR_WAVE)*DT/DVY;               // Positive change to j
      dFdt[findex(i+1,j  ,k+1)] -= (INCR_WAVE+CORR_WAVE)*DT/DVY;               // Negative change to j-1

      // Transverse propagation, Fx updates:
      INCR_WAVE = -0.5*DT/DVY*dF;
      CORR_WAVE *= DT/DVY;
      cell.cpu_fy[BLOCK*SIZE_FLUXS+accIndex(i,j,k)] = INCR_WAVE*fabs(AX)*fabs(AY) + CORR_WAVE*fabs(AX);
      //CORR_WAVE = 0.0;
      dFdt[findex(i+2,j+1,k+1)] += (INCR_WAVE*AX_POS*AY_POS + CORR_WAVE*AX_POS)*DT/DVX; // Positive change to i
      dFdt[findex(i+1,j+1,k+1)] -= (INCR_WAVE*AX_POS*AY_POS + CORR_WAVE*AX_POS)*DT/DVX; // Negative change to i-1
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE*AX_NEG*AY_POS + CORR_WAVE*AX_NEG)*DT/DVX;
      dFdt[findex(i  ,j+1,k+1)] -= (INCR_WAVE*AX_NEG*AY_POS + CORR_WAVE*AX_NEG)*DT/DVX;
      dFdt[findex(i+2,j  ,k+1)] += (INCR_WAVE*AX_POS*AY_NEG - CORR_WAVE*AX_POS)*DT/DVX;
      dFdt[findex(i+1,j  ,k+1)] -= (INCR_WAVE*AX_POS*AY_NEG - CORR_WAVE*AX_POS)*DT/DVX;
      dFdt[findex(i+1,j  ,k+1)] += (INCR_WAVE*AX_NEG*AY_NEG - CORR_WAVE*AX_NEG)*DT/DVX;
      dFdt[findex(i  ,j  ,k+1)] -= (INCR_WAVE*AX_NEG*AY_NEG - CORR_WAVE*AX_NEG)*DT/DVX;

      // Transverse propagation, Fz updates:
      dFdt[findex(i+1,j+1,k+2)] += (INCR_WAVE*AY_POS*AZ_POS + CORR_WAVE*AZ_POS)*DT/DVZ; // Positive chnage to k
      dFdt[findex(i+1,j+1,k+1)] -= (INCR_WAVE*AY_POS*AZ_POS + CORR_WAVE*AZ_POS)*DT/DVZ; // Negative change to k-1
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE*AY_POS*AZ_NEG + CORR_WAVE*AZ_NEG)*DT/DVZ;
      dFdt[findex(i+1,j+1,k  )] -= (INCR_WAVE*AY_POS*AZ_NEG + CORR_WAVE*AZ_NEG)*DT/DVZ;
      dFdt[findex(i+1,j  ,k+2)] += (INCR_WAVE*AY_NEG*AZ_POS - CORR_WAVE*AZ_POS)*DT/DVZ;
      dFdt[findex(i+1,j  ,k+1)] -= (INCR_WAVE*AY_NEG*AZ_POS - CORR_WAVE*AZ_POS)*DT/DVZ;
      dFdt[findex(i+1,j  ,k+1)] += (INCR_WAVE*AY_NEG*AZ_NEG - CORR_WAVE*AZ_NEG)*DT/DVZ;
      dFdt[findex(i+1,j  ,k  )] -= (INCR_WAVE*AY_NEG*AZ_NEG - CORR_WAVE*AZ_NEG)*DT/DVZ;

      INCR_WAVE = DT*DT/DVY/DVZ*dF/6.0;
      
      // Double transverse propagation, Fx updates:
      dFdt[findex(i+2,j+1,k+1)] += INCR_WAVE*AX_POS*AY_POS*fabs(AZ)*DT/DVX; // Ax > 0, Ay > 0
      dFdt[findex(i+1,j+1,k+1)] += INCR_WAVE*AX_NEG*AY_POS*fabs(AZ)*DT/DVX;
      dFdt[findex(i+2,j  ,k+1)] += INCR_WAVE*AX_POS*AY_NEG*fabs(AZ)*DT/DVX;
      dFdt[findex(i+1,j  ,k+1)] += INCR_WAVE*AX_NEG*AY_NEG*fabs(AZ)*DT/DVX;
      
      dFdt[findex(i+2,j+1,k+2)] -= INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVX;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i+1,j+1,k+2)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+2,j  ,k+2)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i+1,j  ,k+2)] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+2,j+1,k  )] += INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVX;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i+1,j+1,k  )] += INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVX;   // Ax < 0
      dFdt[findex(i+2,j  ,k  )] += INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i+1,j  ,k  )] += INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVX;   // Ax < 0
      
      dFdt[findex(i+1,j+1,k+1)] -= INCR_WAVE*AX_POS*AY_POS*fabs(AZ)*DT/DVX; // Ax > 0, Ay > 0
      dFdt[findex(i  ,j+1,k+1)] -= INCR_WAVE*AX_NEG*AY_POS*fabs(AZ)*DT/DVX;
      dFdt[findex(i+1,j  ,k+1)] -= INCR_WAVE*AX_POS*AY_NEG*fabs(AZ)*DT/DVX;
      dFdt[findex(i  ,j  ,k+1)] -= INCR_WAVE*AX_NEG*AY_NEG*fabs(AZ)*DT/DVX;
      
      dFdt[findex(i+1,j+1,k+2)] += INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVX;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+1,k+2)] += INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+1,j  ,k+2)] += INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k+2)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+1,j+1,k  )] -= INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVX;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i  ,j+1,k  )] -= INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVX;   // Ax < 0
      dFdt[findex(i+1,j  ,k  )] -= INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVX;   // Ax > 0, Ay < 0, Az < 0
      dFdt[findex(i  ,j  ,k  )] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVX;   // Ax < 0
      
      INCR_WAVE = DT*DT/DVX/DVY*dF/6.0;
      
      // Double transverse propagation, Fz updates:
      dFdt[findex(i+1,j+1,k+2)] += INCR_WAVE*fabs(AX)*AY_POS*AZ_POS*DT/DVZ; // Ay > 0, Az > 0
      dFdt[findex(i+1,j  ,k+2)] += INCR_WAVE*fabs(AX)*AY_NEG*AZ_POS*DT/DVZ;
      dFdt[findex(i+1,j+1,k+1)] += INCR_WAVE*fabs(AX)*AY_POS*AZ_NEG*DT/DVZ;
      dFdt[findex(i+1,j  ,k+1)] += INCR_WAVE*fabs(AX)*AY_NEG*AZ_NEG*DT/DVZ;
      
      dFdt[findex(i+2,j+1,k+2)] -= INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVZ;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+1,k+2)] += INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVZ;   // Ax < 0
      dFdt[findex(i+2,j  ,k+2)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVZ;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k+2)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVZ;   // Ax < 0
      dFdt[findex(i+2,j+1,k+1)] -= INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVZ;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i  ,j+1,k+1)] += INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVZ;   // Ax < 0
      dFdt[findex(i+2,j  ,k+1)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVZ;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k+1)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVZ;   // Ax < 0
      
      dFdt[findex(i+1,j+1,k+1)] -= INCR_WAVE*fabs(AX)*AY_POS*AZ_POS*DT/DVZ; // Ay > 0, Az > 0
      dFdt[findex(i+1,j  ,k+1)] -= INCR_WAVE*fabs(AX)*AY_NEG*AZ_POS*DT/DVZ;
      dFdt[findex(i+1,j+1,k  )] -= INCR_WAVE*fabs(AX)*AY_POS*AZ_NEG*DT/DVZ;
      dFdt[findex(i+1,j  ,k  )] -= INCR_WAVE*fabs(AX)*AY_NEG*AZ_NEG*DT/DVZ;
      
      dFdt[findex(i+2,j+1,k+1)] += INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVZ;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+1,k+1)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVZ;   // Ax < 0
      dFdt[findex(i+2,j  ,k+1)] += INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVZ;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k+1)] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVZ;   // Ax < 0
      dFdt[findex(i+2,j+1,k  )] += INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVZ;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i  ,j+1,k  )] -= INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVZ;   // Ax < 0
      dFdt[findex(i+2,j  ,k  )] += INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVZ;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k  )] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVZ;
   
      // ***********************************
      // ***** INTERFACE BETWEEN K,K-1 *****
      // ***********************************

      // Calculate acceleration at face (i,j,k-1/2):
      VX = blockParams[BlockParams::VXCRD] + (i+0.5)*blockParams[BlockParams::DVX];
      VY = blockParams[BlockParams::VYCRD] + (j+0.5)*blockParams[BlockParams::DVY];
      VZ = blockParams[BlockParams::VZCRD];
      AX = Parameters::q_per_m * (VY*cellParams[CellParams::BZ] - VZ*cellParams[CellParams::BY]);
      AY = Parameters::q_per_m * (VZ*cellParams[CellParams::BX] - VX*cellParams[CellParams::BZ]);
      AZ = Parameters::q_per_m * (VX*cellParams[CellParams::BY] - VY*cellParams[CellParams::BX]);

      AX_NEG = std::min(convert<REAL>(0.0),AX);
      AX_POS = std::max(convert<REAL>(0.0),AX);
      AY_NEG = std::min(convert<REAL>(0.0),AY);
      AY_POS = std::max(convert<REAL>(0.0),AY);
      AZ_NEG = std::min(convert<REAL>(0.0),AZ);
      AZ_POS = std::max(convert<REAL>(0.0),AZ);
      
      // Calculate slope-limited y-derivative of F:
      if (AZ > 0.0) theta_up = avgs[fullInd(i+2,j+2,k+1)] - avgs[fullInd(i+2,j+2,k  )];
      else theta_up = avgs[fullInd(i+2,j+2,k+3)] - avgs[fullInd(i+2,j+2,k+2)];
      theta_lo = avgs[fullInd(i+2,j+2,k+2)] - avgs[fullInd(i+2,j+2,k+1)] + EPSILON;
      theta = limiter(theta_up/theta_lo);

      // Donor Cell Upwind methods, Fz updates:
      dF = avgs[fullInd(i+2,j+2,k+2)] - avgs[fullInd(i+2,j+2,k+1)]; // F(i,j,k) - F(i,j,k-1), jump in F
      INCR_WAVE = AZ_POS*avgs[fullInd(i+2,j+2,k+1)] + AZ_NEG*avgs[fullInd(i+2,j+2,k+2)];
      CORR_WAVE = 0.5*fabs(AZ)*(1.0-fabs(AZ)*DT/DVZ)*dF*theta;
      //CORR_WAVE = 0.0;
      cell.cpu_d1z[BLOCK*SIZE_DERIV+accIndex(i,j,k)] = theta;
      cell.cpu_d2z[BLOCK*SIZE_DERIV+accIndex(i,j,k)] = dF*theta;
      
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE+CORR_WAVE)*DT/DVZ;               // Positive change to k
      dFdt[findex(i+1,j+1,k  )] -= (INCR_WAVE+CORR_WAVE)*DT/DVZ;               // Negative change to k-1

      // Transverse propagation, Fx updates:
      INCR_WAVE = -0.5*DT/DVZ*dF;
      CORR_WAVE *= DT/DVZ;
      cell.cpu_fz[BLOCK*SIZE_FLUXS+accIndex(i,j,k)] = dF;//INCR_WAVE*fabs(AX)*fabs(AZ) + CORR_WAVE*fabs(AX);
      //CORR_WAVE = 0.0;
      dFdt[findex(i+2,j+1,k+1)] += (INCR_WAVE*AX_POS*AZ_POS + CORR_WAVE*AX_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k+1)] -= (INCR_WAVE*AX_POS*AZ_POS + CORR_WAVE*AX_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE*AX_NEG*AZ_POS + CORR_WAVE*AX_NEG)*DT/DVZ;
      dFdt[findex(i  ,j+1,k+1)] -= (INCR_WAVE*AX_NEG*AZ_POS + CORR_WAVE*AX_NEG)*DT/DVZ;
      dFdt[findex(i+2,j+1,k  )] += (INCR_WAVE*AX_POS*AZ_NEG - CORR_WAVE*AX_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k  )] -= (INCR_WAVE*AX_POS*AZ_NEG - CORR_WAVE*AX_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k  )] += (INCR_WAVE*AX_NEG*AZ_NEG - CORR_WAVE*AX_NEG)*DT/DVZ;
      dFdt[findex(i  ,j+1,k  )] -= (INCR_WAVE*AX_NEG*AZ_NEG - CORR_WAVE*AX_NEG)*DT/DVZ;

      dFdt[findex(i+1,j+2,k+1)] += (INCR_WAVE*AY_POS*AZ_POS + CORR_WAVE*AY_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k+1)] -= (INCR_WAVE*AY_POS*AZ_POS + CORR_WAVE*AY_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k+1)] += (INCR_WAVE*AY_NEG*AZ_POS + CORR_WAVE*AY_NEG)*DT/DVZ;
      dFdt[findex(i+1,j+0,k+1)] -= (INCR_WAVE*AY_NEG*AZ_POS + CORR_WAVE*AY_NEG)*DT/DVZ;
      dFdt[findex(i+1,j+2,k  )] += (INCR_WAVE*AY_POS*AZ_NEG - CORR_WAVE*AY_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k  )] -= (INCR_WAVE*AY_POS*AZ_NEG - CORR_WAVE*AY_POS)*DT/DVZ;
      dFdt[findex(i+1,j+1,k  )] += (INCR_WAVE*AY_NEG*AZ_NEG - CORR_WAVE*AY_NEG)*DT/DVZ;
      dFdt[findex(i+1,j+0,k  )] -= (INCR_WAVE*AY_NEG*AZ_NEG - CORR_WAVE*AY_NEG)*DT/DVZ;

      INCR_WAVE = DT*DT/DVY/DVZ*dF/6.0;
      
      // Double transverse propagation, Fx updates:
      dFdt[findex(i+2,j+1,k+1)] += INCR_WAVE*AX_POS*fabs(AY)*AZ_POS*DT/DVX; // Ax > 0, Az > 0
      dFdt[findex(i+1,j+1,k+1)] += INCR_WAVE*AX_NEG*fabs(AY)*AZ_POS*DT/DVX;
      dFdt[findex(i+2,j+1,k  )] += INCR_WAVE*AX_POS*fabs(AY)*AZ_NEG*DT/DVX;
      dFdt[findex(i+1,j+1,k  )] += INCR_WAVE*AX_NEG*fabs(AY)*AZ_NEG*DT/DVX;
      
      dFdt[findex(i+2,j+2,k+1)] -= INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVX;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i+1,j+2,k+1)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+2,j  ,k+1)] += INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i+1,j  ,k+1)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+2,j+2,k  )] -= INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVX;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i+1,j+2,k  )] -= INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVX;   // Ax < 0
      dFdt[findex(i+2,j  ,k  )] += INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i+1,j  ,k  )] += INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVX;   // Ax < 0
      
      dFdt[findex(i+1,j+1,k+1)] -= INCR_WAVE*AX_POS*fabs(AY)*AZ_POS*DT/DVX; // Ax > 0, Az > 0
      dFdt[findex(i  ,j+1,k+1)] -= INCR_WAVE*AX_NEG*fabs(AY)*AZ_POS*DT/DVX;
      dFdt[findex(i+1,j+1,k  )] -= INCR_WAVE*AX_POS*fabs(AY)*AZ_NEG*DT/DVX;
      dFdt[findex(i  ,j+1,k  )] -= INCR_WAVE*AX_NEG*fabs(AY)*AZ_NEG*DT/DVX;
      
      dFdt[findex(i+1,j+2,k+1)] += INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVX;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+2,k+1)] += INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+1,j  ,k+1)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k+1)] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVX;   // Ax < 0
      dFdt[findex(i+1,j+2,k  )] += INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVX;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i  ,j+2,k  )] += INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVX;   // Ax < 0
      dFdt[findex(i+1,j  ,k  )] -= INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVX;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k  )] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVX;   // Ax < 0
      
      INCR_WAVE = DT*DT/DVX/DVZ*dF/6.0;
      
      // Double transverse propagation, Fy updates:
      dFdt[findex(i+1,j+2,k+1)] += INCR_WAVE*fabs(AX)*AY_POS*AZ_POS*DT/DVY; // Ay > 0, Az > 0
      dFdt[findex(i+1,j+1,k+1)] += INCR_WAVE*fabs(AX)*AY_NEG*AZ_POS*DT/DVY;
      dFdt[findex(i+1,j+2,k  )] += INCR_WAVE*fabs(AX)*AY_POS*AZ_NEG*DT/DVY;
      dFdt[findex(i+1,j+1,k  )] += INCR_WAVE*fabs(AX)*AY_NEG*AZ_NEG*DT/DVY;
      
      dFdt[findex(i+2,j+2,k+1)] -= INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVY;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+2,k+1)] += INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVY;   // Ax < 0
      dFdt[findex(i+2,j+1,k+1)] -= INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVY;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j+1,k+1)] += INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVY;   // Ax < 0
      dFdt[findex(i+2,j+2,k  )] -= INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVY;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i  ,j+2,k  )] += INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVY;   // Ax < 0
      dFdt[findex(i+2,j+1,k  )] -= INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVY;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j+1,k  )] += INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVY;   // Ax < 0
      
      dFdt[findex(i+1,j+1,k+1)] -= INCR_WAVE*fabs(AX)*AY_POS*AZ_POS*DT/DVY; // Ay > 0, Az > 0
      dFdt[findex(i+1,j  ,k+1)] -= INCR_WAVE*fabs(AX)*AY_NEG*AZ_POS*DT/DVY;
      dFdt[findex(i+1,j+1,k  )] -= INCR_WAVE*fabs(AX)*AY_POS*AZ_NEG*DT/DVY;
      dFdt[findex(i+1,j  ,k  )] -= INCR_WAVE*fabs(AX)*AY_NEG*AZ_NEG*DT/DVY;
      
      dFdt[findex(i+2,j+1,k+1)] += INCR_WAVE*AX_POS*AY_POS*AZ_POS*DT/DVY;   // Ax > 0, Ay > 0, Az > 0
      dFdt[findex(i  ,j+1,k+1)] -= INCR_WAVE*AX_NEG*AY_POS*AZ_POS*DT/DVY;   // Ax < 0
      dFdt[findex(i+2,j  ,k+1)] += INCR_WAVE*AX_POS*AY_NEG*AZ_POS*DT/DVY;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k+1)] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_POS*DT/DVY;   // Ax < 0
      dFdt[findex(i+2,j+1,k  )] += INCR_WAVE*AX_POS*AY_POS*AZ_NEG*DT/DVY;   // Ax > 0, Ay > 0, Az < 0
      dFdt[findex(i  ,j+1,k  )] -= INCR_WAVE*AX_NEG*AY_POS*AZ_NEG*DT/DVY;   // Ax < 0
      dFdt[findex(i+2,j  ,k  )] += INCR_WAVE*AX_POS*AY_NEG*AZ_NEG*DT/DVY;   // Ax > 0, Ay < 0
      dFdt[findex(i  ,j  ,k  )] -= INCR_WAVE*AX_NEG*AY_NEG*AZ_NEG*DT/DVY;   // Ax < 0
   }
   accumulateChanges(BLOCK,dFdt,cell.cpu_fx,cell.cpu_nbrsVel+BLOCK*SIZE_NBRS_VEL);
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVel(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   for (UINT i=0; i<WID3; ++i) cell.cpu_avgs[BLOCK*SIZE_VELBLOCK+i] += cell.cpu_fx[BLOCK*SIZE_VELBLOCK+i];
}
#endif
