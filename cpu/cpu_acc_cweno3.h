/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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
#include "cell_spatial.h"
#include "cpu_common.h"
#include "project.h"

const float C1o2   = 0.5;
const float C1o8   = 1.0/8.0;
const float C2     = 2.0;
const float C3     = 3.0;
const float C4     = 4.0;
const float C3o10  = 3.0/10.0;
const float C3o5   = 3.0/5.0;
const float C1o10  = 1.0/10.0;
const float C13o12 = 13.0/12.0;
const float C1o3   = 1.0/3.0;
const float C1o4   = 1.0/4.0;
const float C5o6   = 5.0/6.0;
const float C1o6   = 1.0/6.0;
const float C7o6   = 7.0/6.0;
const float C11o6  = 11.0/6.0;
const float C13o3  = 13.0/3.0;
const float EPS    = 1.0e-6;

template<typename T> T accIndex(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}
template<typename T> T fullInd(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> T isBoundary(const T& STATE,const T& BND) {return STATE & BND;}

template<typename CELL> void copyAverages(CELL& cell) {
   for (unsigned int i=0; i<SIZE_VELBLOCK; ++i) cell.cpu_d1y[i] = cell.cpu_avgs[i];
}

template<typename CELL> void sumAverages(CELL& cell) {
   for (unsigned int i=0; i<SIZE_VELBLOCK; ++i) cell.cpu_avgs[i] = 0.5*(cell.cpu_avgs[i]+cell.cpu_d1y[i]);
}

template<typename UINT,typename REAL> REAL accelerationX(const UINT& j,const UINT& k,const REAL* const cellParams,const REAL* const blockParams) {
   const REAL EX = cellParams[CellParams::EX];
   const REAL BY = cellParams[CellParams::BY];
   const REAL BZ = cellParams[CellParams::BZ];
   const REAL VY = blockParams[BlockParams::VYCRD] + (j+C1o2)*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (k+C1o2)*blockParams[BlockParams::DVZ];
   return Parameters::q_per_m*(EX + VY*BZ - VZ*BY);
}

template<typename UINT,typename REAL> REAL accelerationY(const UINT& i,const UINT& k,const REAL* const cellParams,const REAL* const blockParams) {
   const REAL EY = cellParams[CellParams::EY];
   const REAL BX = cellParams[CellParams::BX];
   const REAL BZ = cellParams[CellParams::BZ];
   const REAL VX = blockParams[BlockParams::VXCRD] + (i+C1o2)*blockParams[BlockParams::DVX];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (k+C1o2)*blockParams[BlockParams::DVZ];
   return Parameters::q_per_m*(EY + VZ*BX - VX*BZ);
}

template<typename UINT,typename REAL> REAL accelerationZ(const UINT& i,const UINT& j,const REAL* const cellParams,const REAL* const blockParams) {
   const REAL EZ = cellParams[CellParams::EZ];
   const REAL BX = cellParams[CellParams::BX];
   const REAL BY = cellParams[CellParams::BY];
   const REAL VX = blockParams[BlockParams::VXCRD] + (i+C1o2)*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (j+C1o2)*blockParams[BlockParams::DVY];
   return Parameters::q_per_m*(EZ + VX*BY - VY*BX);
}

template<typename T> T cwenoReconstruct3(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2,T& x_left,T& x_right) {
   // Calculate smoothness factors:
   //const T beta0 = C13o12*(xcc-C2*xr1+xr2)*(xcc-C2*xr1+xr2) + C1o4*(C3*xcc-C4*xr1+xr2)*(C3*xcc-C4*xr1+xr2);
   //const T beta1 = C13o12*(xl1-C2*xcc+xr1)*(xl1-C2*xcc+xr1) + C1o4*(xl1+xr1)*(xl1+xr1);
   //const T beta2 = C13o12*(xl2-C2*xl1+xcc)*(xl2-C2*xl1+xcc) + C1o4*(xl2-C4*xl1+C3*xcc)*(xl2-C4*xl1+C3*xcc);

   const T beta_L = (xcc-xl1)*(xcc-xl1);
   const T beta_C = C13o3*(xr1-2*xcc+xl1)*(xr1-2*xcc+xl1) + C1o4*(xr1-xl1)*(xr1-xl1);
   const T beta_R = (xr1-xcc)*(xr1-xcc);
   
   // Calculate weight factors:
   //const T alpha0 = C3o10 / ((EPS+beta0)*(EPS+beta0));
   //const T alpha1 = C3o5  / ((EPS+beta1)*(EPS+beta1));
   //const T alpha2 = C1o10 / ((EPS+beta2)*(EPS+beta2));
   
   const T alpha_L = C1o4 / ((EPS + beta_L)*(EPS + beta_L));
   const T alpha_C = C1o2 / ((EPS + beta_C)*(EPS + beta_C));
   const T alpha_R = C1o4 / ((EPS + beta_R)*(EPS + beta_R));
   
   //const T omega0 = alpha0 / (alpha0+alpha1+alpha2);
   //const T omega1 = alpha1 / (alpha0+alpha1+alpha2);
   //const T omega2 = alpha2 / (alpha0+alpha1+alpha2);
   
   const T omega_L = alpha_L / (alpha_L + alpha_C + alpha_R);
   const T omega_C = alpha_C / (alpha_L + alpha_C + alpha_R);
   const T omega_R = alpha_R / (alpha_L + alpha_C + alpha_R);
   
   // Calculate left/right face 3rd order CWENO reconstructed values:
   /*
   x_right = omega0*(C1o3*xcc  + C5o6*xr1 - C1o6*xr2)
           + omega1*(-C1o6*xl1 + C5o6*xcc + C1o3*xr1)
           + omega2*(C1o3*xl2  - C7o6*xl1 + C11o6*xcc);
   
   x_left  = omega2*(C11o6*xcc - C7o6*xr1 + C1o3*xr2)
           + omega1*(C1o3*xl1  + C5o6*xcc - C1o6*xr1)
           + omega0*(-C1o6*xl2 + C5o6*xl1 + C1o3*xcc);
   */
   x_right = omega_L*(-0.5*xl1 + 1.5*xcc)
           + omega_R*(0.5*xcc + 0.5*xr1)
           + omega_C*(-1.0*xl1 + 4.0*xcc + 5.0*xr1)/12.0;
   
   x_left = omega_L*(0.5*xl1 + 0.5*xcc)
          + omega_R*(1.5*xcc - 0.5*xr1)
          + omega_C*(5.0*xl1 + 4.0*xcc - 1.0*xr1)/12.0;
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

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivsX(CELL& cell,const UINT& BLOCK) {
   // Copy volume averages and ghost cells values to a temporary array,
   // which is easier to use to calculate derivatives:
   REAL* const d1x = cell.cpu_d1x + BLOCK*SIZE_DERIV; // Left face value at cell i,j,k
   REAL* const d2x = cell.cpu_d2x + BLOCK*SIZE_DERIV; // Right face value at cell i,j,k
   
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      // Reconstruct left/right vx face values using CWENO:
      xl2 = avgs[fullInd(i  ,j+2,k+2)];
      xl1 = avgs[fullInd(i+1,j+2,k+2)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+3,j+2,k+2)];
      xr2 = avgs[fullInd(i+4,j+2,k+2)];
      cwenoReconstruct3(xl2,xl1,xcc,xr1,xr2,d1x[accIndex(i,j,k)],d2x[accIndex(i,j,k)]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivsY(CELL& cell,const UINT& BLOCK) {
   // Copy volume averages and ghost cells values to a temporary array,
   // which is easier to use to calculate derivatives:
   REAL* const d1x = cell.cpu_d1x + BLOCK*SIZE_DERIV; // Left face value at cell i,j,k
   REAL* const d2x = cell.cpu_d2x + BLOCK*SIZE_DERIV; // Right face value at cell i,j,k
   
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      // Reconstruct left/right vy face values using CWENO:
      xl2 = avgs[fullInd(i+2,j  ,k+2)];
      xl1 = avgs[fullInd(i+2,j+1,k+2)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+2,j+3,k+2)];
      xr2 = avgs[fullInd(i+2,j+4,k+2)];      
      cwenoReconstruct3(xl2,xl1,xcc,xr1,xr2,d1x[accIndex(i,j,k)],d2x[accIndex(i,j,k)]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivsZ(CELL& cell,const UINT& BLOCK) {
   // Copy volume averages and ghost cells values to a temporary array,
   // which is easier to use to calculate derivatives:
   REAL* const d1x = cell.cpu_d1x + BLOCK*SIZE_DERIV; // Left face value at cell i,j,k
   REAL* const d2x = cell.cpu_d2x + BLOCK*SIZE_DERIV; // Right face value at cell i,j,k
   
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      // Reconstruct left/right vz face values using CWENO:
      xl2 = avgs[fullInd(i+2,j+2,k  )];
      xl1 = avgs[fullInd(i+2,j+2,k+1)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+2,j+2,k+3)];
      xr2 = avgs[fullInd(i+2,j+2,k+4)];	 
      cwenoReconstruct3(xl2,xl1,xcc,xr1,xr2,d1x[accIndex(i,j,k)],d2x[accIndex(i,j,k)]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesX(CELL& cell,const UINT& BLOCK) {
   // The size of array is (5,4,4)
   const UINT XS = 5;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   
   REAL* const Fx = cell.cpu_fx + BLOCK*SIZE_FLUXS;
   const REAL* const cellParams = cell.cpu_cellParams;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use to calculate fluxes:
   REAL left[WID3 + WID2];
   REAL rght[WID3 + WID2];
   fetchDerivsX(BLOCK,left,cell.cpu_d1x,cell.cpu_nbrsVel);
   fetchDerivsX(BLOCK,rght,cell.cpu_d2x,cell.cpu_nbrsVel);
   
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      //const REAL AX = accelerationX(j,k,cellParams,blockParams);
      const REAL AX = 1.0;
      
      Fx[accIndex(i,j,k)] = C1o2*AX*(left[k*YSXS+j*XS+i+1] + rght[k*YSXS+j*XS+i  ])
	- C1o2*fabs(AX)*(left[k*YSXS+j*XS+i+1] - rght[k*YSXS+j*XS+i  ]);
   }
}
   
template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesY(CELL& cell,const UINT& BLOCK) {
   // The size of the array avgs is (4,5,4)
   const UINT XS = 4;
   const UINT YS = 5;
   const UINT YSXS = YS*XS;

   REAL* const Fx = cell.cpu_fy + BLOCK*SIZE_FLUXS;
   const REAL* const cellParams = cell.cpu_cellParams;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   
   REAL left[SIZE_DERIV + SIZE_BDERI];
   REAL rght[SIZE_DERIV + SIZE_BDERI];
   fetchDerivsY(BLOCK,left,cell.cpu_d1x,cell.cpu_nbrsVel);
   fetchDerivsY(BLOCK,rght,cell.cpu_d2x,cell.cpu_nbrsVel);
       
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      //const REAL AY = accelerationY(i,k,cellParams,blockParams);
      const REAL AY = 0.0;
      
      Fx[accIndex(i,j,k)] = C1o2*AY*(left[k*YSXS+(j+1)*XS+i] + rght[k*YSXS+j*XS+i  ])
	- C1o2*fabs(AY)*(left[k*YSXS+(j+1)*XS+i] - rght[k*YSXS+j*XS+i  ]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesZ(CELL& cell,const UINT& BLOCK) {
   // The size of array avgs is (4,4,5)
   const UINT XS = 4;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;

   REAL* const Fx = cell.cpu_fz + BLOCK*SIZE_FLUXS;
   const REAL* const cellParams = cell.cpu_cellParams;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   
   REAL left[SIZE_DERIV + SIZE_BDERI];
   REAL rght[SIZE_DERIV + SIZE_BDERI];
   fetchDerivsZ(BLOCK,left,cell.cpu_d1x,cell.cpu_nbrsVel);
   fetchDerivsZ(BLOCK,rght,cell.cpu_d2x,cell.cpu_nbrsVel);
   
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      //const REAL AZ = accelerationZ(i,k,cellParams,blockParams);
      const REAL AZ = 0.0;
      
      Fx[accIndex(i,j,k)] = C1o2*AZ*(left[(k+1)*YSXS+j*XS+i] + rght[k*YSXS+j*XS+i  ])
	- C1o2*fabs(AZ)*(left[(k+1)*YSXS+j*XS+i] - rght[k*YSXS+j*XS+i  ]);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVelX(CELL& cell,const UINT& BLOCK,REAL* arr,const REAL& DT) {
   REAL* avgs = arr + BLOCK*SIZE_VELBLOCK;
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   const UINT XS = 5;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   const REAL cnst = DT / blockParams[BlockParams::DVX];
   
   fetchFluxesX(BLOCK,flux,cell.cpu_fx,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[k*YSXS+j*XS+(i+1)])*cnst;
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVelY(CELL& cell,const UINT& BLOCK,REAL* arr,const REAL& DT) {
   REAL* avgs = arr + BLOCK*SIZE_VELBLOCK;
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   const UINT XS = 4;
   const UINT YS = 5;
   const UINT YSXS = YS*XS;
   const REAL cnst = DT / blockParams[BlockParams::DVY];
   
   fetchFluxesY(BLOCK,flux,cell.cpu_fy,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[k*YSXS+(j+1)*XS+i])*cnst;
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVelZ(CELL& cell,const UINT& BLOCK,REAL* arr,const REAL& DT) {
   REAL* avgs = arr + BLOCK*SIZE_VELBLOCK;
   REAL flux[SIZE_FLUXS + SIZE_BFLUX];
   
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   const UINT XS = 4;
   const UINT YS = 4;
   const UINT YSXS = YS*XS;
   const REAL cnst = DT / blockParams[BlockParams::DVZ];
   
   fetchFluxesZ(BLOCK,flux,cell.cpu_fz,cell.cpu_nbrsVel);
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[(k+1)*YSXS+j*XS+i])*cnst;
   }
}


#endif
