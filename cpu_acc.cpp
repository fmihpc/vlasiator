#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "definitions.h"
#include "common.h"
#include "cell_spatial.h"
#include "grid.h"
#include "project.h"
#include "cpu_common.h"

using namespace std;

uint accIndex(cuint& i,cuint& j,cuint& k) {return k*WID2+j*WID+i;}
inline uint fullInd(cuint& i,cuint& j,cuint& k) {return k*64+j*8+i;}

real velDerivs1(creal& xl2,creal& xl1,creal& xcc,creal& xr1,creal& xr2) {
   return superbee(xl1,xcc,xr1);
   //return 0.5*(xr1-xl1);
}

real velDerivs2(creal& xl2,creal& xl1,creal& xcc,creal& xr1,creal& xr2) {
   return convert<real>(0.0);
   //return xr1 - 2.0*xcc - xl1;
}

void fetchAllAverages(cuint& BLOCK,real* avgs,real* cpu_avgs,uint* nbrsVel) {
   uint nbrBlock;

   for (uint i=0; i<8*WID3; ++i) avgs[i] = 0.0;
   
   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VX_NEG_BND == NbrsVel::VX_NEG_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<2; ++i) {
	 avgs[fullInd(i  ,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<2; ++i) {
	 avgs[fullInd(i  ,j+2,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i+2,j,k)];
      }
   }

   // Copy averages from +x neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VX_POS_BND == NbrsVel::VX_POS_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<2; ++i) {
	 avgs[fullInd(i+6,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXPOS];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<2; ++i) {
	 avgs[fullInd(i+6,j+2,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i  ,j,k)];
      }
   }
   
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VY_NEG_BND == NbrsVel::VY_NEG_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<2; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j  ,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<2; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j  ,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j+2,k)];
      }
   }
   // Copy averages from +y neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VY_POS_BND == NbrsVel::VY_POS_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<2; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+6,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYPOS];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<2; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+6,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j  ,k)];
      }
   }
   
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VZ_NEG_BND == NbrsVel::VZ_NEG_BND) {
      for (uint k=0; k<2; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k  )] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      for (uint k=0; k<2; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k  )] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,k+2)];
      }
   }
   // Copy averages from +z neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VZ_NEG_BND == NbrsVel::VZ_NEG_BND) {
      for (uint k=0; k<2; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k+6)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZPOS];
      for (uint k=0; k<2; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k+6)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,k)];
      }
   }

   // Copy volume averages of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[fullInd(i+2,j+2,k+2)] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
      //avgs[fullInd(i+2,j+2,k+2)] = accIndex(i,j,k);
   }
}

void fetchAveragesX(cuint& BLOCK,real* avgs,real* cpu_avgs,uint* nbrsVel) {
   // The size of array avgs is (5,4,4).
   cuint XS=5;
   cuint YS=4;
   cuint YSXS = YS*XS;
   uint nbrBlock;
   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VX_NEG_BND == NbrsVel::VX_NEG_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) {
	 avgs[k*YSXS+j*XS] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) {
	 avgs[k*YSXS+j*XS] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(3,j,k)];
      }
   }
   // Copy volume averages of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[k*YSXS+j*XS+(i+1)] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
   }
}

void fetchAveragesY(cuint& BLOCK,real* avgs,real* cpu_avgs,uint* nbrsVel) {
   // The size of array avgs (4,5,4).
   cuint XS=4;
   cuint YS=5;
   cuint YSXS = YS*XS;
   uint nbrBlock;
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VY_NEG_BND == NbrsVel::VY_NEG_BND) {
      for (uint k=0; k<WID; ++k) for (uint i=0; i<WID; ++i) {
	 avgs[k*YSXS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      for (uint k=0; k<WID; ++k) for (uint i=0; i<WID; ++i) {
	 avgs[k*YSXS+i] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,3,k)];
      }
   }
   // Copy volume averages of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[k*YSXS+(j+1)*XS+i] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
   }
}

void fetchAveragesZ(cuint& BLOCK,real* avgs,real* cpu_avgs,uint* nbrsVel) {
   // The size of array avgs (4,4,5).
   cuint XS=4;
   cuint YS=4;
   cuint YSXS=YS*XS;
   uint nbrBlock;
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VZ_NEG_BND == NbrsVel::VZ_NEG_BND) {
      for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[j*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 avgs[j*XS+i] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,3)];
      }
   }
   // Copy volume averages of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[(k+1)*YSXS+j*XS+i] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
   }
}

void fetchDerivsX(cuint& BLOCK,real* d1x,real* cpu_d1x,uint* nbrsVel) {
   // The size of array avgs (5,4,4).
   cuint XS=5;
   cuint YS=4;
   cuint YSXS=YS*XS;
   uint nbrBlock;
   // Copy derivatives from -x neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VX_NEG_BND == NbrsVel::VX_NEG_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) {
	 d1x[k*YSXS+j*XS] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) {
	 d1x[k*YSXS+j*XS] = cpu_d1x[nbrBlock*SIZE_DERIV + accIndex(3,j,k)];
      }
   }
   // Copy derivatives for this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      d1x[k*YSXS+j*XS+(i+1)] = cpu_d1x[BLOCK*SIZE_DERIV + accIndex(i,j,k)];
   }
}

void fetchDerivsY(cuint& BLOCK,real* d1y,real* cpu_d1y,uint* nbrsVel) {
   // The size of array avgs (4,5,4).
   cuint XS=4;
   cuint YS=5;
   cuint YSXS=YS*XS;
   uint nbrBlock;
   // Copy derivatives from -y neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VY_NEG_BND == NbrsVel::VY_NEG_BND) {
      for (uint k=0; k<WID; ++k) for (uint i=0; i<WID; ++i) {
	 d1y[k*YSXS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      for (uint k=0; k<WID; ++k) for (uint i=0; i<WID; ++i) {
	 d1y[k*YSXS+i] = cpu_d1y[nbrBlock*SIZE_DERIV + accIndex(i,3,k)];
      }
   }
   // Copy derivatives for this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      d1y[k*YSXS+(j+1)*XS+i] = cpu_d1y[BLOCK*SIZE_DERIV + accIndex(i,j,k)];
   }
}

void fetchDerivsZ(cuint& BLOCK,real* d1z,real* cpu_d1z,uint* nbrsVel) {
   // The size of array avgs (4,4,5)
   cuint XS=4;
   cuint YS=4;
   cuint YSXS=YS*XS;
   uint nbrBlock;
   // Copy derivatives from -z neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VZ_NEG_BND == NbrsVel::VZ_NEG_BND) {
      for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 d1z[j*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 d1z[j*XS+i] = cpu_d1z[nbrBlock*SIZE_DERIV + accIndex(i,j,3)];
      }
   }
   // Copy derivatives for this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      d1z[(k+1)*YSXS+j*XS+i] = cpu_d1z[BLOCK*SIZE_DERIV + accIndex(i,j,k)];
   }
}

void fetchFluxesX(cuint& BLOCK,real* flux,real* cpu_fx,uint* nbrsVel) {
   cuint XS = 5;
   cuint YS = 4;
   cuint YSXS = YS*XS;
   uint nbrBlock;
   // Fetch fluxes from +x neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VX_POS_BND == NbrsVel::VX_POS_BND) {
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) {
	 flux[k*YSXS+j*XS+4] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXPOS];
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) {
	 flux[k*YSXS+j*XS+4] = cpu_fx[nbrBlock*SIZE_FLUXS + accIndex(0,j,k)];
      }
   }
   // Fetch the fluxes of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      flux[k*YSXS+j*XS+i] = cpu_fx[BLOCK*SIZE_FLUXS + accIndex(i,j,k)];
   }
}

void fetchFluxesY(cuint& BLOCK,real* flux,real* cpu_fy,uint* nbrsVel) {
   cuint XS = 4;
   cuint YS = 5;
   cuint YSXS = YS*XS;
   uint nbrBlock;
   // Fetch fluxes from +y neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VY_POS_BND == NbrsVel::VY_POS_BND) {
      for (uint k=0; k<WID; ++k) for (uint i=0; i<WID; ++i) {
	 flux[k*YSXS+4*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYPOS];
      for (uint k=0; k<WID; ++k) for (uint i=0; i<WID; ++i) {
	 flux[k*YSXS+4*XS+i] = cpu_fy[nbrBlock*SIZE_FLUXS + accIndex(i,0,k)];
      }
   }   
   // Fetch the fluxes of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      flux[k*YSXS+j*XS+i] = cpu_fy[BLOCK*SIZE_FLUXS + accIndex(i,j,k)];
   }
}

void fetchFluxesZ(cuint& BLOCK,real* flux,real* cpu_fz,uint* nbrsVel) {
   cuint XS = 4;
   cuint YS = 4;
   cuint YSXS = YS*XS;
   uint nbrBlock;
   // Fetch fluxes from +z neighbour, or calculate using a boundary function:
   if (nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE] & NbrsVel::VZ_POS_BND == NbrsVel::VZ_POS_BND) {
      for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 flux[4*YSXS+j*XS+i] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZPOS];
      for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
	 flux[4*YSXS+j*XS+i] = cpu_fz[nbrBlock*SIZE_FLUXS + accIndex(i,j,0)];
      }
   }
   // Fetch the fluxes of this block:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      flux[k*YSXS+j*XS+i] = cpu_fz[BLOCK*SIZE_FLUXS + accIndex(i,j,k)];
   }
}

void cpu_calcVelDerivs(cuint& BLOCK,Grid& grid) {

   // Copy volume averages and ghost cells values to a temporary array, 
   // which is easier to use when calculating derivatives:
   real* d1x = grid.getD1x();
   real* d1y = grid.getD1y();
   real* d1z = grid.getD1z();
   real* d2x = grid.getD2x();
   real* d2y = grid.getD2y();
   real* d2z = grid.getD2z();
   
   real avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,grid.getBlockArray(),grid.getNbrsVel());
   
   // Calculate derivatives for each cell in the block:
   real xl2,xl1,xcc,xr1,xr2;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Calculate 1st & 2nd derivatives to vx-direction:
      xl2 = avgs[fullInd(i  ,j+2,k+2)];
      xl1 = avgs[fullInd(i+1,j+2,k+2)];
      xcc = avgs[fullInd(i+2,j+2,k+2)];
      xr1 = avgs[fullInd(i+3,j+2,k+2)];
      xr2 = avgs[fullInd(i+4,j+2,k+2)];
      d1x[BLOCK*SIZE_DERIV + accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2x[BLOCK*SIZE_DERIV + accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
      // Calculate 1st & 2nd derivatives to vy-direction:
      xl2 = avgs[fullInd(i+2,j  ,k+2)];
      xl1 = avgs[fullInd(i+2,j+1,k+2)];
      xr1 = avgs[fullInd(i+2,j+3,k+2)];
      xr2 = avgs[fullInd(i+2,j+4,k+2)];
      d1y[BLOCK*SIZE_DERIV + accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2y[BLOCK*SIZE_DERIV + accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
      // Calculate 1st & 2nd derivatives to vz-direction:
      xl2 = avgs[fullInd(i+2,j+2,k  )];
      xl1 = avgs[fullInd(i+2,j+2,k+1)];
      xr1 = avgs[fullInd(i+2,j+2,k+3)];
      xr2 = avgs[fullInd(i+2,j+2,k+4)];
      d1z[BLOCK*SIZE_DERIV + accIndex(i,j,k)] = velDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2z[BLOCK*SIZE_DERIV + accIndex(i,j,k)] = velDerivs2(xl2,xl1,xcc,xr1,xr2);
   }
}

void calcVelFluxesX(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   // The size of array avgs is (5,4,4)
   real* fx = grid.getFx();
   cuint XS = 5;
   cuint YS = 4;
   cuint YSXS = YS*XS;   
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use when calculating fluxes:
   real avgs[SIZE_VELBLOCK + WID2];
   real d1x[SIZE_DERIV + SIZE_BDERI];
   real d2x[SIZE_DERIV + SIZE_BDERI];
   fetchAveragesX(BLOCK,avgs,grid.getBlockArray(),grid.getNbrsVel());
   fetchDerivsX(BLOCK,d1x,grid.getD1x(),grid.getNbrsVel());
   fetchDerivsX(BLOCK,d2x,grid.getD2x(),grid.getNbrsVel());
   
   // Reconstruct volume averages at negative and positive side of the face,
   // and calculate vx-flux:
   real avg_neg,avg_pos;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avg_neg = reconstruct_neg(avgs[k*YSXS+j*XS+(i  )],d1x[k*YSXS+j*XS+(i  )],d2x[k*YSXS+j*XS+(i  )]);
      avg_pos = reconstruct_pos(avgs[k*YSXS+j*XS+(i+1)],d1x[k*YSXS+j*XS+(i+1)],d2x[k*YSXS+j*XS+(i+1)]);
      fx[BLOCK*SIZE_FLUXS + accIndex(i,j,k)] = velocityFluxX(i,j,k,avg_neg,avg_pos,grid.getCellParams()+SPATCELL*SIZE_CELLPARAMS,
							  grid.getBlockParams()+BLOCK*SIZE_BLOCKPARAMS);
   }
}

void calcVelFluxesY(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   // The size of array avgs is (4,5,4)
   real* fy = grid.getFy();
   cuint XS = 4;
   cuint YS = 5;
   cuint YSXS = YS*XS;   
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use when calculating fluxes:
   real avgs[SIZE_VELBLOCK + WID2];
   real d1y[SIZE_DERIV + SIZE_BDERI];
   real d2y[SIZE_DERIV + SIZE_BDERI];
   fetchAveragesY(BLOCK,avgs,grid.getBlockArray(),grid.getNbrsVel());
   fetchDerivsY(BLOCK,d1y,grid.getD1y(),grid.getNbrsVel());
   fetchDerivsY(BLOCK,d2y,grid.getD2y(),grid.getNbrsVel());
   
   // Reconstruct volume averages at negative and positive side of the face,
   // and calculate vy-flux:
   real avg_neg,avg_pos;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avg_neg = reconstruct_neg(avgs[k*YSXS+(j  )*XS+i],d1y[k*YSXS+(j  )*XS+i],d2y[k*YSXS+(j  )*XS+i]);
      avg_pos = reconstruct_pos(avgs[k*YSXS+(j+1)*XS+i],d1y[k*YSXS+(j+1)*XS+i],d2y[k*YSXS+(j+1)*XS+i]);
      fy[BLOCK*SIZE_FLUXS + accIndex(i,j,k)] = velocityFluxY(i,j,k,avg_neg,avg_pos,grid.getCellParams()+SPATCELL*SIZE_CELLPARAMS,grid.getBlockParams()+BLOCK*SIZE_BLOCKPARAMS);
   }
}

void calcVelFluxesZ(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   // The size of array avgs is (4,4,5)
   real* fz = grid.getFz();
   cuint XS = 4;
   cuint YS = 4;
   cuint YSXS = YS*XS;   
   // Copy volume averages and fluxes to temporary arrays, which are easier
   // to use when calculating fluxes:
   real avgs[SIZE_VELBLOCK + WID2];
   real d1z[SIZE_DERIV + SIZE_BDERI];
   real d2z[SIZE_DERIV + SIZE_BDERI];
   fetchAveragesZ(BLOCK,avgs,grid.getBlockArray(),grid.getNbrsVel());
   fetchDerivsZ(BLOCK,d1z,grid.getD1z(),grid.getNbrsVel());
   fetchDerivsZ(BLOCK,d2z,grid.getD2z(),grid.getNbrsVel());
   
   // Reconstruct volume averages at negative and positive side of the face,
   // and calculate vz-flux:
   real avg_neg,avg_pos;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avg_neg = reconstruct_neg(avgs[(k  )*YSXS+j*XS+i],d1z[(k  )*YSXS+j*XS+i],d2z[(k  )*YSXS+j*XS+i]);
      avg_pos = reconstruct_pos(avgs[(k+1)*YSXS+j*XS+i],d1z[(k+1)*YSXS+j*XS+i],d2z[(k+1)*YSXS+j*XS+i]);
      fz[BLOCK*SIZE_FLUXS + accIndex(i,j,k)] = velocityFluxZ(i,j,k,avg_neg,avg_pos,grid.getCellParams()+SPATCELL*SIZE_CELLPARAMS,grid.getBlockParams()+BLOCK*SIZE_BLOCKPARAMS);
   }
}

void propagateVel(cuint& SPATCELL,cuint& BLOCK,Grid& grid,creal& DT) {
   real avgs[SIZE_VELBLOCK];
   real flux[SIZE_FLUXS + SIZE_BFLUX];
   
   // Calculate the contribution from vx-fluxes:
   uint XS = 5;
   uint YS = 4;
   uint YSXS = YS*XS;
   real cnst = DT / grid.getBlockParams()[BLOCK*SIZE_BLOCKPARAMS + BlockParams::DVX];
   fetchFluxesX(BLOCK,flux,grid.getFx(),grid.getNbrsVel());
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] = (flux[k*YSXS+j*XS+i] - flux[k*YSXS+j*XS+(i+1)])*cnst;
   }
   
   // Calculate the contribution from vy-fluxes:
   XS = 4;
   YS = 5;
   YSXS = YS*XS;
   cnst = DT / grid.getBlockParams()[BLOCK*SIZE_BLOCKPARAMS + BlockParams::DVY];
   fetchFluxesY(BLOCK,flux,grid.getFy(),grid.getNbrsVel());
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[k*YSXS+(j+1)*XS+i])*cnst;
   }

   // Calculate the contribution from vz-fluxes:
   XS = 4;
   YS = 4;
   YSXS = YS*XS;
   cnst = DT / grid.getBlockParams()[BLOCK*SIZE_BLOCKPARAMS + BlockParams::DVZ];
   fetchFluxesZ(BLOCK,flux,grid.getFz(),grid.getNbrsVel());
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      avgs[accIndex(i,j,k)] += (flux[k*YSXS+j*XS+i] - flux[(k+1)*YSXS+j*XS+i])*cnst;
   }
   
   // Store results:
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      grid.getBlockArray()[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)] += avgs[accIndex(i,j,k)];
   }
}

bool cpu_acceleration(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT) {
   bool success = true;

   // Calculate derivatives:
   #pragma omp parallel for
   for (uint block=cell.velBlockIndex; block<cell.velBlockIndex+cell.N_blocks; ++block) {
      cpu_calcVelDerivs(block,grid);
   }

   // Calculate fluxes:
   #pragma omp parallel for
   for (uint block=cell.velBlockIndex; block<cell.velBlockIndex+cell.N_blocks; ++block) {
      calcVelFluxesX(cellIndex,block,grid);
      calcVelFluxesY(cellIndex,block,grid);
      calcVelFluxesZ(cellIndex,block,grid);
   }
      
   // Propagate volume averages:
   #pragma omp parallel for
   for (uint block=cell.velBlockIndex; block<cell.velBlockIndex+cell.N_blocks; ++block) {
      propagateVel(cellIndex,block,grid,DT);
   }

   return success;
}
