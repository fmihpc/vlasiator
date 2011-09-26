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

#include "../definitions.h"
#include "../common.h"
#include "../cell_spatial.h"
#include "cpu_common.h"
#include "../project.h"

template<typename T> T accIndex(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}
template<typename T> T fullInd(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> T isBoundary(const T& STATE,const T& BND) {return STATE & BND;}

template<typename T> T velDerivs1(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2) {
   return superbee(xl1,xcc,xr1);
}

template<typename T> T velDerivs2(const T& xl2,const T& xl1,const T& xcc,const T& xr1,const T& xr2) {
   return convert<T>(0.0);
}

template<typename REAL,typename UINT> void fetchAllAverages(const UINT& BLOCK,REAL* const avgs,const REAL* const cpu_avgs,const UINT* const nbrsVel) {
   UINT nbrBlock;
   for (UINT i=0; i<8*WID3; ++i) avgs[i] = 0.0;
   const UINT STATE = nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS];
   
   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<2; ++i) {
	 avgs[fullInd(i  ,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK; 
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<2; ++i) {
	    avgs[fullInd(i  ,j+2,k+2)] = tmp[accIndex(i+2,j,k)];
	 }
      }
   }
   // Copy averages from +x neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<2; ++i) {
	 avgs[fullInd(i+6,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<2; ++i) {
	    avgs[fullInd(i+6,j+2,k+2)] = tmp[accIndex(i  ,j,k)];
	 }
      }
   }
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j  ,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j  ,k+2)] = tmp[accIndex(i,j+2,k)];
	 }
      }
   }
   // Copy averages from +y neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+6,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j+6,k+2)] = tmp[accIndex(i,j  ,k)];
	 }
      }
   }
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k  )] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j+2,k  )] = tmp[accIndex(i,j,k+2)];
	 }
      }
   }
   // Copy averages from +z neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k+6)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    avgs[fullInd(i+2,j+2,k+6)] = tmp[accIndex(i,j,k)];
	 }
      }
   }
   // Copy volume averages of this block:
   creal* const tmp = cpu_avgs + BLOCK*SIZE_VELBLOCK;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
      #pragma ivdep
      for (UINT i=0; i<WID; ++i) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VX_NEG_BND) == 0) {
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
   const uint XS=4;
   const uint YS=5;
   const uint YSXS = YS*XS;
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VY_NEG_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VZ_NEG_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VX_NEG_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VY_NEG_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VZ_NEG_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VX_POS_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VY_POS_BND) == 0) {
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
   if (isBoundary(nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::NBRFLAGS],NbrsVel::VZ_POS_BND) == 0) {
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

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelDerivs(CELL& cell,const UINT& BLOCK) {
   // Copy volume averages and ghost cells values to a temporary array,
   // which is easier to use to calculate derivatives:
   REAL* const d1x = cell.cpu_d1x + BLOCK*SIZE_DERIV;
   REAL* const d1y = cell.cpu_d1y + BLOCK*SIZE_DERIV;
   REAL* const d1z = cell.cpu_d1z + BLOCK*SIZE_DERIV;
   REAL* const d2x = cell.cpu_d2x + BLOCK*SIZE_DERIV;
   REAL* const d2y = cell.cpu_d2y + BLOCK*SIZE_DERIV;
   REAL* const d2z = cell.cpu_d2z + BLOCK*SIZE_DERIV;
   
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);
   
   REAL xl2,xl1,xcc,xr1,xr2;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      // Calculate 1st & 2nd derivatives to vx-direction:
      xl2 = avgs[fullInd(i  ,j+2,k+2)];
      xl1 = avgs[fullInd(i+1,j+2,k+2)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+3,j+2,k+2)];
      xr2 = avgs[fullInd(i+4,j+2,k+2)];
      d1x[accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2x[accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
      // Calculate 1st & 2nd derivatives to vy-direction:
      xl2 = avgs[fullInd(i+2,j  ,k+2)];
      xl1 = avgs[fullInd(i+2,j+1,k+2)];
      xr1 = avgs[fullInd(i+2,j+3,k+2)];
      xr2 = avgs[fullInd(i+2,j+4,k+2)];
      d1y[accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2y[accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
      // Calculate 1st & 2nd derivatives to vz-direction:
      xl2 = avgs[fullInd(i+2,j+2,k  )];
      xl1 = avgs[fullInd(i+2,j+2,k+1)];
      xr1 = avgs[fullInd(i+2,j+2,k+3)];
      xr2 = avgs[fullInd(i+2,j+2,k+4)];
      d1z[accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2z[accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesX(CELL& cell,const UINT& BLOCK) {
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
   
   // Reconstruct volume averages at negative and positive side of the face,
   // and calculate vx-flux for each cell in the block:
   REAL avg_neg,avg_pos;
   REAL* const fx           = cell.cpu_fx          + BLOCK*SIZE_FLUXS;
   creal* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   for (UINT k=0; k<WID; ++k) {
      REAL K = convert<REAL>(1.0)*k;
      for (UINT j=0; j<WID; ++j) {
	 REAL J = convert<REAL>(1.0)*j;
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    avg_neg = reconstruct_neg(avgs[k*YSXS+j*XS+(i  )],d1x[k*YSXS+j*XS+(i  )],d2x[k*YSXS+j*XS+(i  )]);
	    avg_pos = reconstruct_pos(avgs[k*YSXS+j*XS+(i+1)],d1x[k*YSXS+j*XS+(i+1)],d2x[k*YSXS+j*XS+(i+1)]);
	    fx[accIndex(i,j,k)] = velocityFluxX(J,K,avg_neg,avg_pos,cell.cpu_cellParams,blockParams);
	 }
      }
   }
}
   
template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesY(CELL& cell,const UINT& BLOCK) {
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
       
   // Reconstruct volume averages at negative and positive side of the face,
   // and calculate vy-flux for each cell in the block:
   REAL avg_neg,avg_pos;
   REAL* const fy           = cell.cpu_fy          + BLOCK*SIZE_FLUXS;
   creal* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   for (UINT k=0; k<WID; ++k) {
      const REAL K = convert<REAL>(1.0)*k;
      for (UINT j=0; j<WID; ++j) {
	 for (UINT i=0; i<WID; ++i) {
	    const REAL I = convert<REAL>(1.0)*i;
	    avg_neg = reconstruct_neg(avgs[k*YSXS+(j  )*XS+i],d1y[k*YSXS+(j  )*XS+i],d2y[k*YSXS+(j  )*XS+i]);
	    avg_pos = reconstruct_pos(avgs[k*YSXS+(j+1)*XS+i],d1y[k*YSXS+(j+1)*XS+i],d2y[k*YSXS+(j+1)*XS+i]);
	    fy[accIndex(i,j,k)] = velocityFluxY(I,K,avg_neg,avg_pos,cell.cpu_cellParams,blockParams);
	 }
      }
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxesZ(CELL& cell,const UINT& BLOCK) {
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
   
   // Reconstruct volume averages at negative and positive side of the face,
   // and calculate vz-flux:
   REAL avg_neg,avg_pos;
   REAL* const fz           = cell.cpu_fz          + BLOCK*SIZE_FLUXS;
   creal* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS;
   for (UINT k=0; k<WID; ++k) {
      for (UINT j=0; j<WID; ++j) {
	 const REAL J = convert<REAL>(1.0)*j;
	 for (UINT i=0; i<WID; ++i) {
	    const REAL I = convert<REAL>(1.0)*i;
	    avg_neg = reconstruct_neg(avgs[(k  )*YSXS+j*XS+i],d1z[(k  )*YSXS+j*XS+i],d2z[(k  )*YSXS+j*XS+i]);
	    avg_pos = reconstruct_pos(avgs[(k+1)*YSXS+j*XS+i],d1z[(k+1)*YSXS+j*XS+i],d2z[(k+1)*YSXS+j*XS+i]);
	    fz[accIndex(i,j,k)] = velocityFluxZ(I,J,avg_neg,avg_pos,cell.cpu_cellParams,blockParams);
	 }
      }
   }      
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
