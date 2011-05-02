#ifndef CPU_ACC_H
#define CPU_ACC_H

#include "definitions.h"
#include "common.h"
#include "cell_spatial.h"
#include "cpu_common.h"
#include "project.h"

template<typename T> T accIndex(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}
template<typename T> T fullInd(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> T isBoundary(const T& STATE,const T& BND) {return STATE & BND;}

template<typename REAL,typename UINT> REAL acceleration_x(const UINT& j,const UINT& k,const REAL* const cellParams,const REAL* const blockParams) {
   const REAL BY = cellParams[CellParams::BY];
   const REAL BZ = cellParams[CellParams::BZ];
   const REAL VY = blockParams[BlockParams::VYCRD] + (j+0.5)*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (k+0.5)*blockParams[BlockParams::DVZ];
   return Parameters::q_per_m * (VY*BZ - VZ*BY);
}

template<typename REAL,typename UINT> REAL acceleration_y(const UINT& i,const UINT& k,const REAL* const cellParams,const REAL* const blockParams) {
   const REAL BX = cellParams[CellParams::BX];
   const REAL BZ = cellParams[CellParams::BZ];
   const REAL VX = blockParams[BlockParams::VXCRD] + (i+0.5)*blockParams[BlockParams::DVX];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (k+0.5)*blockParams[BlockParams::DVZ];
   return Parameters::q_per_m * (VZ*BX - VX*BZ);
}

template<typename REAL,typename UINT> REAL acceleration_z(const UINT& i,const UINT& j,const REAL* const cellParams,const REAL* const blockParams) {
   const REAL BX = cellParams[CellParams::BX];
   const REAL BY = cellParams[CellParams::BY];
   const REAL VX = blockParams[BlockParams::VXCRD] + (i+0.5)*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (j+0.5)*blockParams[BlockParams::DVY];
   return Parameters::q_per_m * (VX*BY - VY*BX);
}

template<typename T> T ppm_integral_L(const T& X,const T& AVG,const T& AVG_L,const T& AVG_R) {
   const T ONE = 1.0;
   const T TWO_THIRDS = convert<T>(2.0)/convert<T>(3.0);
   const T SIX = 6.0;
   const T HALF = 0.5;
   
   const T delta_a = AVG_R-AVG_L;
   const T a6      = SIX*(AVG - HALF*(AVG_L + AVG_R));
   return AVG_R - HALF*X*(delta_a - (ONE-TWO_THIRDS*X)*a6);
}

template<typename T> T ppm_integral_R(const T& X,const T& AVG,const T& AVG_L,const T& AVG_R) {
   const T ONE = 1.0;
   const T TWO_THIRDS = convert<T>(2.0)/convert<T>(3.0);
   const T SIX = 6.0;
   const T HALF = 0.5;
   
   const T delta_a = AVG_R-AVG_L;
   const T a6      = SIX*(AVG - HALF*(AVG_L + AVG_R));
   return AVG_L + HALF*X*(delta_a + (ONE-TWO_THIRDS*X)*a6);
}

template<typename T> T calcDelta_a(const T& xl1,const T& xcc,const T& xr1,const T& dx) {
   return (xr1-convert<T>(2.0)*xcc-xl1)/6.0/(dx*dx);
}

//template<typename T> T velDerivs1(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2) {
template<typename T> void velDerivs(T& x_L,T& x_R,const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2,const T& dx) {
   // Calculate first approximations of left/right face values:
   x_L = convert<T>(0.5)*(xl1 + xcc) - (limiter_ppm(xl1,xcc,xr1) - limiter_ppm(xl2,xl1,xcc))/convert<T>(6.0);
   x_R = convert<T>(0.5)*(xcc + xr1) - (limiter_ppm(xcc,xr1,xr2) - limiter_ppm(xl1,xcc,xr1))/convert<T>(6.0);
   /*
   // Discontinuity correction:
   const T EPSILON = 0.01;
   const T ETA1    = 20.0;
   const T ETA2    = 0.05;
   
   const T DELTA_A_NEG = calcDelta_a(xl2,xl1,xcc,dx);
   const T DELTA_A_POS = calcDelta_a(xcc,xr1,xr2,dx);
   const T CONST = fabs(xr1-xl1) - EPSILON * std::min(fabs(xr1),fabs(xl1));
   T eta_worm = 0.0;
   if (DELTA_A_POS*DELTA_A_NEG < 0.0 && CONST > 0.0) {
      eta_worm = dx*dx*(DELTA_A_NEG-DELTA_A_POS)/(xr1-xl1);
   }
   const T ZERO = 0.0;
   const T HALF = 0.5;
   const T ONE  = 1.0;
   const T ETA_J = std::max(ZERO,std::min(ETA1*(eta_worm-ETA2),ONE));
   
   const T x_LD = xl1 + HALF*limiter_ppm(xl2,xl1,xcc);
   const T x_RD = xr1 - HALF*limiter_ppm(xcc,xr1,xr2);
   x_L = x_L*(ONE-ETA_J) + x_LD*ETA_J;
   x_R = x_R*(ONE-ETA_J) + x_RD*ETA_J;
   */
   // Apply monotonicity correction:
   if ((xr1 - xcc)*(xcc - xl1) < convert<T>(0.0)) {
   //if ((x_R - xcc)*(xcc - x_L) < convert<T>(0.0)) {
      x_L = xcc;
      x_R = xcc;
   } else {
      const T value1 = (x_R - x_L)*(xcc - convert<T>(0.5)*(x_L+x_R));
      const T value2 = (x_R - x_L)*(x_R - x_L)/convert<T>(6.0);
      if (value1 > value2) x_L = convert<T>(3.0)*xcc - convert<T>(2.0)*x_R;
      else if (value1 < -value2) x_R = convert<T>(3.0)*xcc - convert<T>(2.0)*x_L;
   }
   //return x_L;
}
/*
template<typename T> T velDerivs2(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2) {
   // Calculate first approximations of left/right face values:
   T x_L = convert<T>(0.5)*(xl1 + xcc) - (limiter_ppm(xl1,xcc,xr1) - limiter_ppm(xl2,xl1,xcc))/convert<T>(6.0);
   T x_R = convert<T>(0.5)*(xcc + xr1) - (limiter_ppm(xcc,xr1,xr2) - limiter_ppm(xl1,xcc,xr1))/convert<T>(6.0);
   
   // Discontinuity correction should be done:
   
   
   // Apply monotonicity correction:
   if ((x_R - xcc)*(xcc - x_L) < convert<T>(0.0)) {
      x_L = xcc;
      x_R = xcc;
   } else {
      const T value1 = (x_R - x_L)*(xcc - convert<T>(0.5)*(x_L+x_R));
      const T value2 = (x_R - x_L)*(x_R - x_L)/convert<T>(6.0);
      if (value1 > value2) x_L = convert<T>(3.0)*xcc - convert<T>(2.0)*x_R;
      else if (value1 < -value2) x_R = convert<T>(3.0)*xcc - convert<T>(2.0)*x_L;
   }
   return x_R;
}
*/
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

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivsX(CELL& cell,const UINT& BLOCK) {
   // Copy volume averages and ghost cells values to a temporary array,
   // which is easier to use to calculate derivatives:
   REAL* const d1x = cell.cpu_d1x + BLOCK*SIZE_DERIV;
   REAL* const d2x = cell.cpu_d2x + BLOCK*SIZE_DERIV;

   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      xl2 = avgs[fullInd(i  ,j+2,k+2)];
      xl1 = avgs[fullInd(i+1,j+2,k+2)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+3,j+2,k+2)];
      xr2 = avgs[fullInd(i+4,j+2,k+2)];
      //d1x[accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      //d2x[accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
      velDerivs(d1x[accIndex(i,j,k)],d2x[accIndex(i,j,k)],xl2,xl1,xcc,xr1,xr2,cell.cpu_blockParams[BLOCK*SIZE_BLOCKPARAMS+BlockParams::DVX]);
   }
}
 
template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivsY(CELL& cell,const UINT& BLOCK) {
   REAL* const d1y = cell.cpu_d1y + BLOCK*SIZE_DERIV;
   REAL* const d2y = cell.cpu_d2y + BLOCK*SIZE_DERIV;
   
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
      
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      xl2 = avgs[fullInd(i+2,j  ,k+2)];
      xl1 = avgs[fullInd(i+2,j+1,k+2)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+2,j+3,k+2)];
      xr2 = avgs[fullInd(i+2,j+4,k+2)];
      //d1y[accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      //d2y[accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
      velDerivs(d1y[accIndex(i,j,k)],d2y[accIndex(i,j,k)],xl2,xl1,xcc,xr1,xr2,cell.cpu_blockParams[BLOCK*SIZE_BLOCKPARAMS+BlockParams::DVY]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivsZ(CELL& cell,const UINT& BLOCK) {
   REAL* const d1z = cell.cpu_d1z + BLOCK*SIZE_DERIV;
   REAL* const d2z = cell.cpu_d2z + BLOCK*SIZE_DERIV;
   
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      xl2 = avgs[fullInd(i+2,j+2,k  )];
      xl1 = avgs[fullInd(i+2,j+2,k+1)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+2,j+2,k+3)];
      xr2 = avgs[fullInd(i+2,j+2,k+4)];
      velDerivs(d1z[accIndex(i,j,k)],d2z[accIndex(i,j,k)],xl2,xl1,xcc,xr1,xr2,cell.cpu_blockParams[BLOCK*SIZE_BLOCKPARAMS+BlockParams::DVZ]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesX(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   // The size of array is (5,4,4)
   const UINT XS = 5;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use to calculate fluxes:
   REAL avgs[SIZE_VELBLOCK + WID2];
   REAL d1x[SIZE_DERIV + SIZE_BDERI];
   REAL d2x[SIZE_DERIV + SIZE_BDERI];
   fetchAveragesX(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   fetchDerivsX(BLOCK,d1x,cell.cpu_d1x,cell.cpu_nbrsVel);
   fetchDerivsX(BLOCK,d2x,cell.cpu_d2x,cell.cpu_nbrsVel);

   REAL* const fx          = cell.cpu_fx          + BLOCK*SIZE_FLUXS;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS; // Does not work for refined blocks
   const REAL DT_PER_DVX = DT/blockParams[BlockParams::DVX];
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      const REAL AX = acceleration_x(j,k,cell.cpu_cellParams,blockParams);
      if (AX > 0.0) {
	 const REAL avg   = avgs[k*YSXS+j*XS+i  ]; // Values are from iv-1
	 const REAL avg_L =  d1x[k*YSXS+j*XS+i  ];
	 const REAL avg_R =  d2x[k*YSXS+j*XS+i  ];
	 fx[accIndex(i,j,k)] = AX * ppm_integral_L( AX*DT_PER_DVX,avg,avg_L,avg_R);
      } else {
	 const REAL avg   = avgs[k*YSXS+j*XS+i+1]; // Values are from iv 
	 const REAL avg_L =  d1x[k*YSXS+j*XS+i+1];
	 const REAL avg_R =  d2x[k*YSXS+j*XS+i+1];
	 fx[accIndex(i,j,k)] = AX * ppm_integral_R(-AX*DT_PER_DVX,avg,avg_L,avg_R);
      }
   }
}
   
template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesY(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   // The size of the array avgs is (4,5,4)
   const UINT XS = 4;
   const UINT YS = 5;
   const UINT YSXS = YS*XS;
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use to calculate fluxes:
   REAL avgs[SIZE_VELBLOCK + WID2];
   REAL d1y[SIZE_DERIV + SIZE_BDERI];
   REAL d2y[SIZE_DERIV + SIZE_BDERI];
   fetchAveragesY(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   fetchDerivsY(BLOCK,d1y,cell.cpu_d1y,cell.cpu_nbrsVel);
   fetchDerivsY(BLOCK,d2y,cell.cpu_d2y,cell.cpu_nbrsVel);
       
   REAL* const fy           = cell.cpu_fy          + BLOCK*SIZE_FLUXS;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   const REAL DT_PER_DVY = DT/blockParams[BlockParams::DVY];
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      const REAL AY = acceleration_y(i,k,cell.cpu_cellParams,blockParams);
      if (AY > 0.0) {
	 const REAL avg   = avgs[k*YSXS+(j  )*XS+i]; // Values are from jv-1
	 const REAL avg_L =  d1y[k*YSXS+(j  )*XS+i];
	 const REAL avg_R =  d2y[k*YSXS+(j  )*XS+i];
	 fy[accIndex(i,j,k)] = AY * ppm_integral_L( AY*DT_PER_DVY,avg,avg_L,avg_R);
      } else {
	 const REAL avg   = avgs[k*YSXS+(j+1)*XS+i]; // Values are from jv
	 const REAL avg_L =  d1y[k*YSXS+(j+1)*XS+i];
	 const REAL avg_R =  d2y[k*YSXS+(j+1)*XS+i];
	 fy[accIndex(i,j,k)] = AY * ppm_integral_R(-AY*DT_PER_DVY,avg,avg_L,avg_R);
      }
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesZ(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   // The size of array avgs is (4,4,5)
   const UINT XS = 4;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;   
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use to calculate fluxes:
   REAL avgs[SIZE_VELBLOCK + WID2];
   REAL d1z[SIZE_DERIV + SIZE_BDERI];
   REAL d2z[SIZE_DERIV + SIZE_BDERI];
   fetchAveragesZ(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   fetchDerivsZ(BLOCK,d1z,cell.cpu_d1z,cell.cpu_nbrsVel);
   fetchDerivsZ(BLOCK,d2z,cell.cpu_d2z,cell.cpu_nbrsVel);
   
   REAL* const fz = cell.cpu_fz + BLOCK*SIZE_FLUXS;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   const REAL DT_PER_DVZ = DT/blockParams[BlockParams::DVZ];
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      const REAL AZ = acceleration_z(i,j,cell.cpu_cellParams,blockParams);
      if (AZ > 0.0) {
	 const REAL avg   = avgs[(k  )*YSXS+j*XS+i]; // Values are from kv
	 const REAL avg_L =  d1z[(k  )*YSXS+j*XS+i];
	 const REAL avg_R =  d2z[(k  )*YSXS+j*XS+i];
	 fz[accIndex(i,j,k)] = AZ * ppm_integral_L(-AZ*DT_PER_DVZ,avg,avg_L,avg_R);
      } else {
	 const REAL avg   = avgs[(k+1)*YSXS+j*XS+i]; // Values are from kv
	 const REAL avg_L =  d1z[(k+1)*YSXS+j*XS+i];
	 const REAL avg_R =  d2z[(k+1)*YSXS+j*XS+i];
	 fz[accIndex(i,j,k)] = AZ * ppm_integral_R(-AZ*DT_PER_DVZ,avg,avg_L,avg_R);
      }
   }      
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVelX(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   REAL* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   
   const UINT XS = 5;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   const REAL cnst = DT / blockParams[BlockParams::DVX];
   fetchFluxesX(BLOCK,flux,cell.cpu_fx,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) 
     avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[k*YSXS+j*XS+(i+1)])*cnst;
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVelY(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   REAL* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;

   const UINT XS = 4;
   const UINT YS = 5;
   const UINT YSXS = YS*XS;
   const REAL cnst = DT / blockParams[BlockParams::DVY];
   fetchFluxesY(BLOCK,flux,cell.cpu_fy,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) 
     avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[k*YSXS+(j+1)*XS+i])*cnst;
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVelZ(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   REAL* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   
   const UINT XS = 4;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   const REAL cnst = DT / blockParams[BlockParams::DVZ];
   fetchFluxesZ(BLOCK,flux,cell.cpu_fz,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i)
     avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[(k+1)*YSXS+j*XS+i])*cnst;
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVel(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   REAL avgs[SIZE_VELBLOCK];
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   
   // Calculate the contribution to d(avg)/dt from vx-fluxes:
   UINT XS = 5;
   UINT YS = 4;
   UINT YSXS = YS*XS;
   REAL cnst = DT / blockParams[BlockParams::DVX];
   fetchFluxesX(BLOCK,flux,cell.cpu_fx,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] = (flux[k*YSXS+j*XS+i] - flux[k*YSXS+j*XS+(i+1)])*cnst;
   }
   // Calculate the contribution to d(avg)/dt from vy-fluxes:
   XS = 4;
   YS = 5;
   YSXS = YS*XS;
   cnst = DT / blockParams[BlockParams::DVY];
   fetchFluxesY(BLOCK,flux,cell.cpu_fy,cell.cpu_nbrsVel);
   //for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
   //   avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[k*YSXS+(j+1)*XS+i])*cnst;
   //}
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) avgs[accIndex(i,j,k)] += flux[k*YSXS+j*XS+i]*cnst;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) avgs[accIndex(i,j,k)] -= flux[k*YSXS+(j+1)*XS+i]*cnst;
   // Calculate the contribution to d(avg)/dt from vz-fluxes:
   XS = 4;
   YS = 4;
   YSXS = YS*XS;
   cnst = DT / blockParams[BlockParams::DVZ];
   fetchFluxesZ(BLOCK,flux,cell.cpu_fz,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[(k+1)*YSXS+j*XS+i])*cnst;
   }
   // Store results:
   REAL* const cpu_avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      cpu_avgs[accIndex(i,j,k)] += avgs[accIndex(i,j,k)];
   }
}
#endif
