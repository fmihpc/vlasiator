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

inline uint trIndex(cuint& i,cuint& j,cuint& k) {return k*WID2+j*WID+i;}
inline uint isBoundary(cuint& STATE,cuint& BND) {return STATE & BND;}

real spatDerivs1(creal& xl2,creal& xl1,creal& xcc,creal& xr1,creal& xr2) {
   return superbee<real>(xl1,xcc,xr1);
   //return vanLeer(xl1,xcc,xr1);
}

real spatDerivs2(creal& xl2,creal& xl1,creal& xcc,creal& xr1,creal& xr2) {
   return convert<real>(0.0);
}

void cpu_calcSpatDerivs(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   creal* const avgs = grid.getBlockArray();
   cuint* const nbrs = grid.getNbrsVel() + BLOCK*SIZE_NBRS_VEL;
   
   real* const d1x = grid.getD1x();
   real* const d2x = grid.getD2x();
   real* const d1y = grid.getD1y();
   real* const d2y = grid.getD2y();
   real* const d1z = grid.getD1z();
   real* const d2z = grid.getD2z();
   
   // Calculate derivatives for each cell in the block:
   real xl2,xl1,xcc,xr1,xr2;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      xl2 = 0.0;
      xcc = avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      xr2 = 0.0;
      // *************** X DERIVATIVES ******************
      // Get vol.avg. from -x neighbour, or calculate it using a boundary function:
      if (nbrs[NbrsVel::STATE] & NbrsVel::X_NEG_BND == NbrsVel::X_NEG_BND) {
	 xl1 = 0.0;
      } else {
	 cuint nbrBlock = grid.getNbrsVel()[BLOCK*SIZE_NBRS_VEL + NbrsVel::XNEG];
	 xl1 = avgs[nbrBlock*SIZE_VELBLOCK + trIndex(i,j,k)];
	 //xl1 = avgs[nbrs[NbrsVel::XNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Get vol.avg. from +x neighbour, or calculate it using a boundary function:
      if (nbrs[NbrsVel::STATE] & NbrsVel::X_POS_BND == NbrsVel::X_POS_BND) {
	 xr1 = 0.0;
      } else {
	 cuint nbrBlock = grid.getNbrsVel()[BLOCK*SIZE_NBRS_VEL + NbrsVel::XPOS];
	 xr1 = avgs[nbrBlock*SIZE_VELBLOCK + trIndex(i,j,k)];
	 //xr1 = avgs[nbrs[NbrsVel::XPOS]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Calculate derivatives in x-direction:
      d1x[BLOCK*SIZE_DERIV + trIndex(i,j,k)] = spatDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2x[BLOCK*SIZE_DERIV + trIndex(i,j,k)] = spatDerivs2(xl2,xl1,xcc,xr1,xr2);
            
      // *************** Y DERIVATIVES ******************
      // Get vol.avg. from -y neighbour or calculate it using a boundary function:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Y_NEG_BND == NbrsVel::Y_NEG_BND) {
	 xl1 = 0.0;
      } else {
	 xl1 = avgs[nbrs[NbrsVel::YNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Get vol.avg. from +y neighbour or calculate it using a boundary function:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Y_POS_BND == NbrsVel::Y_POS_BND) {
	 xr1 = 0.0;
      } else {
	 xr1 = avgs[nbrs[NbrsVel::YPOS]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Calculate derivatives in y-direction:
      d1y[BLOCK*SIZE_DERIV + trIndex(i,j,k)] = spatDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2y[BLOCK*SIZE_DERIV + trIndex(i,j,k)] = spatDerivs2(xl2,xl1,xcc,xr1,xr2);
      
      // *************** Z DERIVATIVES ******************
      // Get vol.avg. from -z neighbour or calculate it using a boundary function:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Z_NEG_BND == NbrsVel::Z_NEG_BND) {
	 xl1 = 0.0;
      } else {
	 xl1 = avgs[nbrs[NbrsVel::ZNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Get vol.avg. from +z neighbour or calculate it using a boundary function:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Z_POS_BND == NbrsVel::Z_POS_BND) {
	 xr1 = 0.0;
      } else {
	 xr1 = avgs[nbrs[NbrsVel::ZPOS]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Calculate derivatives in z-direction:
      d1z[BLOCK*SIZE_DERIV + trIndex(i,j,k)] = spatDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2z[BLOCK*SIZE_DERIV + trIndex(i,j,k)] = spatDerivs2(xl2,xl1,xcc,xr1,xr2);
   }
}

void cpu_calcSpatDerivs(SpatialCell& cell,cuint& BLOCK,std::vector<SpatialCell*>& nbrPtrs) {
   //cuint* const nbrs = cell.cpu_nbrsVel + BLOCK*SIZE_NBRS_VEL;
   //real* const avgs  = cell.cpu_blockArray + BLOCK*SIZE_VELBLOCK;
   real* const d1x   = cell.cpu_d1x + BLOCK*SIZE_DERIV;
   real* const d2x   = cell.cpu_d2x + BLOCK*SIZE_DERIV;
   real* const d1y   = cell.cpu_d1y + BLOCK*SIZE_DERIV;
   real* const d2y   = cell.cpu_d2y + BLOCK*SIZE_DERIV;
   real* const d1z   = cell.cpu_d1z + BLOCK*SIZE_DERIV;
   real* const d2z   = cell.cpu_d2z + BLOCK*SIZE_DERIV;
   cuint STATE = cell.cpu_nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::STATE];

   // Calculate derivatives for each cell in the block:
   real xl2,xl1,xcc,xr1,xr2;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      xcc = cell.cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      // *************** X DERIVATIVES ******************
      xl2 = 0.0;
      xr2 = 0.0;
      // Get vol.avg. from -x neighbour, or calculate it using a boundary function:
      if (nbrPtrs[0] == NULL) {
      //if (isBoundary(STATE,NbrsVel::X_NEG_BND) > 0) {
	 xl1 = 0.0; // Boundary value!
      } else {
	 xl1 = nbrPtrs[0]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Get vol.avg. from +x neighbour, or calculate it using a boundary function:
      if (nbrPtrs[1] == NULL) {
      //if (isBoundary(STATE,NbrsVel::X_POS_BND) > 0) {
	 xr1 = 0.0; // Boundary value!
      } else {
	 xr1 = nbrPtrs[1]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Calculate derivatives in x-direction:
      d1x[trIndex(i,j,k)] = spatDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2x[trIndex(i,j,k)] = spatDerivs2(xl2,xl1,xcc,xr1,xr2);
      
      // *************** Y DERIVATIVES ******************
      // Get vol.avg. from -y neighbour or calculate it using a boundary function:
      if (nbrPtrs[2] == NULL) {
      //if (isBoundary(STATE,NbrsVel::Y_NEG_BND) > 0) {
	 xl1 = 0.0; // Boundary value!
      } else {
	 xl1 = nbrPtrs[2]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Get vol.avg. from +y neighbour or calculate it using a boundary function:
      if (nbrPtrs[3] == NULL) {
      //if (isBoundary(STATE,NbrsVel::Y_POS_BND) > 0) {
	 xr1 = 0.0; // Boundary value!
      } else {
	 xr1 = nbrPtrs[3]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Calculate derivatives in y-direction:
      d1y[trIndex(i,j,k)] = spatDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2y[trIndex(i,j,k)] = spatDerivs2(xl2,xl1,xcc,xr1,xr2);
      
      // *************** Z DERIVATIVES ******************
      // Get vol.avg. from -z neighbour or calculate it using a boundary function:
      if (nbrPtrs[4] == NULL) {
      //if (isBoundary(STATE,NbrsVel::Z_NEG_BND) > 0) {
	 xl1 = 0.0; // Boundary value!
      } else {
	 xl1 = nbrPtrs[4]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Get vol.avg. from +z neighbour or calculate it using a boundary function:
      if (nbrPtrs[5] == NULL) {
      //if (isBoundary(STATE,NbrsVel::Z_POS_BND) > 0) {
	 xr1 = 0.0; // Boundary value!
      } else {
	 xr1 = nbrPtrs[5]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      // Calculate derivatives in z-direction:
      d1z[trIndex(i,j,k)] = spatDerivs1(xl2,xl1,xcc,xr1,xr2);
      d2z[trIndex(i,j,k)] = spatDerivs2(xl2,xl1,xcc,xr1,xr2);      
   }
}

void cpu_calcSpatFluxesX(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   const real* const avgs = grid.getBlockArray();
   const real* const d1x = grid.getD1x();
   const real* const d2x = grid.getD2x();
   const uint* const nbrs = grid.getNbrsVel() + BLOCK*SIZE_NBRS_VEL;
   real* const fx = grid.getFx();
   
   real avg_pos,d1_pos,d2_pos;
   real avg_neg,d1_neg,d2_neg;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Reconstruct volume average at positive side of x-face:
      avg_pos = avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      d1_pos  = d1x[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      d2_pos  = d2x[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      avg_pos = reconstruct_pos(avg_pos,d1_pos,d2_pos);

      // Reconstruct volume average at negative side of x-face.
      // Fetch vol.avg and derivatives from -x neighbour, or calculate them using a boundary func.:
      if (nbrs[NbrsVel::STATE] & NbrsVel::X_NEG_BND == NbrsVel::X_NEG_BND) {
	 avg_neg = 0.0;
	 d1_neg = 0.0;
	 d2_neg = 0.0;
      } else {
	 avg_neg = avgs[nbrs[NbrsVel::XNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d1_neg  = d1x[nbrs[NbrsVel::XNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d2_neg  = d2x[nbrs[NbrsVel::XNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      avg_neg = reconstruct_neg(avg_neg,d1_neg,d2_neg);
      
      // Calculate x-flux:
      fx[BLOCK*SIZE_FLUXS + trIndex(i,j,k)] = spatialFluxX(i,avg_neg,avg_pos,grid.getBlockParams()+BLOCK*SIZE_BLOCKPARAMS);
   }
}

void cpu_calcSpatFluxesX(SpatialCell& cell,cuint& BLOCK,std::vector<SpatialCell*>& nbrPtrs) {
   creal* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   creal* const d1x  = cell.cpu_d1x  + BLOCK*SIZE_DERIV;
   creal* const d2x  = cell.cpu_d2x  + BLOCK*SIZE_DERIV;
   real* const fx    = cell.cpu_fx   + BLOCK*SIZE_FLUXS;
   
   real avg_pos,d1_pos,d2_pos;
   real avg_neg,d1_neg,d2_neg;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Reconstruct volume average at positive side of x-face:
      avg_pos = avgs[trIndex(i,j,k)];
      d1_pos  = d1x[trIndex(i,j,k)];
      d2_pos  = d2x[trIndex(i,j,k)];
      avg_pos = reconstruct_pos(avg_pos,d1_pos,d2_pos);
      
      // Reconstruct volume average at negative side of x-face.
      // Fetch vol.avg and derivatives from -x neighbour, or calculate them using a boundary func.:
      if (nbrPtrs[0] == NULL) {
	 avg_neg = 0.0; // Boundary values !
	 d1_neg = 0.0;
	 d2_neg = 0.0;
      } else {
	 avg_neg = nbrPtrs[0]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d1_neg  = nbrPtrs[0]->cpu_d1x[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
	 d2_neg  = nbrPtrs[0]->cpu_d2x[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      }
      avg_neg = reconstruct_neg(avg_neg,d1_neg,d2_neg);
      
      // Calculate x-flux:
      fx[trIndex(i,j,k)] = spatialFluxX(i,avg_neg,avg_pos,cell.cpu_blockParams+BLOCK*SIZE_BLOCKPARAMS);
      //fx[trIndex(i,j,k)] = avg_neg;
   }
}

void cpu_calcSpatFluxesY(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   creal* const avgs = grid.getBlockArray();
   creal* const d1y = grid.getD1y();
   creal* const d2y = grid.getD2y();
   cuint* const nbrs = grid.getNbrsVel() + BLOCK*SIZE_NBRS_VEL;
   real* const fy = grid.getFy();
   
   real avg_pos,d1_pos,d2_pos;
   real avg_neg,d1_neg,d2_neg;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Reconstruct volume average at positive side of y-face:
      avg_pos = avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      d1_pos  = d1y[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      d2_pos  = d2y[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      avg_pos = reconstruct_pos(avg_pos,d1_pos,d2_pos);
      
      // Reconstruct volume average at negative side of y-face.
      // Fetch vol.avg and derivatives from -y neighbour, or calculate them using a boundary func.:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Y_NEG_BND == NbrsVel::Y_NEG_BND) {
	 avg_neg = 0.0;
	 d1_neg = 0.0;
	 d2_neg = 0.0;
      } else {
	 avg_neg = avgs[nbrs[NbrsVel::YNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d1_neg  = d1y[nbrs[NbrsVel::YNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d2_neg  = d2y[nbrs[NbrsVel::YNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
      }
      avg_neg = reconstruct_neg(avg_neg,d1_neg,d2_neg);
      
      // Calculate y-flux:
      fy[BLOCK*SIZE_FLUXS + trIndex(i,j,k)] = spatialFluxY(j,avg_neg,avg_pos,grid.getBlockParams()+BLOCK*SIZE_BLOCKPARAMS);
   }   
}

void cpu_calcSpatFluxesY(SpatialCell& cell,cuint& BLOCK,std::vector<SpatialCell*>& nbrPtrs) {
   creal* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   creal* const d1y  = cell.cpu_d1y  + BLOCK*SIZE_DERIV;
   creal* const d2y  = cell.cpu_d2y  + BLOCK*SIZE_DERIV;
   real* const fy    = cell.cpu_fy   + BLOCK*SIZE_FLUXS;
   
   real avg_pos,d1_pos,d2_pos;
   real avg_neg,d1_neg,d2_neg;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Reconstruct volume average at positive side of y-face:
      avg_pos = avgs[trIndex(i,j,k)];
      d1_pos  = d1y[trIndex(i,j,k)];
      d2_pos  = d2y[trIndex(i,j,k)];
      avg_pos = reconstruct_pos(avg_pos,d1_pos,d2_pos);

      // Reconstruct volume average at negative side of y-face.
      // Fetch vol.avg and derivatives from -y neighbour, or calculate them using a boundary func.:
      if (nbrPtrs[2] == NULL) {
	 avg_neg = 0.0; // Boundary values !
	 d1_neg = 0.0;
	 d2_neg = 0.0;
      } else {
	 avg_neg = nbrPtrs[2]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d1_neg  = nbrPtrs[2]->cpu_d1y[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
	 d2_neg  = nbrPtrs[2]->cpu_d2y[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      }
      avg_neg = reconstruct_neg(avg_neg,d1_neg,d2_neg);
      
      // Calculate y-flux:
      //fy[trIndex(i,j,k)] = spatialFluxY(j,avg_neg,avg_pos,cell.cpu_blockParams+BLOCK*SIZE_BLOCKPARAMS);
      fy[trIndex(i,j,k)] = spatialFluxY(j,0.0f,0.0f,cell.cpu_blockParams+BLOCK*SIZE_BLOCKPARAMS);
   }
}

void cpu_calcSpatFluxesZ(cuint& SPATCELL,cuint& BLOCK,Grid& grid) {
   creal* const avgs = grid.getBlockArray();
   creal* const d1z = grid.getD1z();
   creal* const d2z = grid.getD2z();
   cuint* const nbrs = grid.getNbrsVel() + BLOCK*SIZE_NBRS_VEL;
   real* const fz = grid.getFz();
   
   real avg_pos,d1_pos,d2_pos;
   real avg_neg,d1_neg,d2_neg;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Reconstruct volume average at positive side of z-face:
      avg_pos = avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
      d1_pos  = d1z[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      d2_pos  = d2z[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      avg_pos = reconstruct_pos(avg_pos,d1_pos,d2_pos);
            
      // Reconstruct volume average at negative side of z-face.
      // Fetch vol.avg and derivatives from -z neighbour, or calculate them using a boundary func.:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Z_NEG_BND == NbrsVel::Z_NEG_BND) {
	 avg_neg = 0.0;
	 d1_neg = 0.0;
	 d2_neg = 0.0;
      } else {
	 avg_neg = avgs[nbrs[NbrsVel::ZNEG]*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d1_neg  = d1z[nbrs[NbrsVel::ZNEG]*SIZE_DERIV + trIndex(i,j,k)];
	 d2_neg  = d2z[nbrs[NbrsVel::ZNEG]*SIZE_DERIV + trIndex(i,j,k)];
      }
      avg_neg = reconstruct_neg(avg_neg,d1_neg,d2_neg);
      
      // Calculate z-flux:
      fz[BLOCK*SIZE_FLUXS + trIndex(i,j,k)] = spatialFluxZ(k,avg_neg,avg_pos,grid.getBlockParams()+BLOCK*SIZE_BLOCKPARAMS);
   }
}

void cpu_calcSpatFluxesZ(SpatialCell& cell,cuint& BLOCK,std::vector<SpatialCell*>& nbrPtrs) {
   creal* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   creal* const d1z  = cell.cpu_d1z  + BLOCK*SIZE_DERIV;
   creal* const d2z  = cell.cpu_d2z  + BLOCK*SIZE_DERIV;
   real* const fz    = cell.cpu_fz   + BLOCK*SIZE_FLUXS;
   
   real avg_pos,d1_pos,d2_pos;
   real avg_neg,d1_neg,d2_neg;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Reconstruct volume average at positive side of z-face:
      avg_pos = avgs[trIndex(i,j,k)];
      d1_pos  = d1z[trIndex(i,j,k)];
      d2_pos  = d2z[trIndex(i,j,k)];
      avg_pos = reconstruct_pos(avg_pos,d1_pos,d2_pos);
      
      // Reconstruct volume average at negative side of z-face.
      // Fetch vol.avg and derivatives from -z neighbour, or calculate them using a boundary func.:
      if (nbrPtrs[4] == NULL) {
	 avg_neg = 0.0; // Boundary values !
	 d1_neg = 0.0;
	 d2_neg = 0.0;
      } else {
	 avg_neg = nbrPtrs[4]->cpu_avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)];
	 d1_neg  = nbrPtrs[4]->cpu_d1z[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
	 d2_neg  = nbrPtrs[4]->cpu_d2z[BLOCK*SIZE_DERIV + trIndex(i,j,k)];
      }
      avg_neg = reconstruct_neg(avg_neg,d1_neg,d2_neg);
      
      // Calculate z-flux:
      fz[trIndex(i,j,k)] = spatialFluxZ(k,avg_neg,avg_pos,cell.cpu_blockParams+BLOCK*SIZE_BLOCKPARAMS);
      //fz[trIndex(i,j,k)] = avg_pos;
   }
}

void cpu_propagateSpat(cuint& SPATCELL,cuint& BLOCK,Grid& grid,creal& DT) {
   real* const avgs = grid.getBlockArray();
   creal* const fx = grid.getFx();
   creal* const fy = grid.getFy();
   creal* const fz = grid.getFz();
   cuint* const nbrs = grid.getNbrsVel() + BLOCK*SIZE_NBRS_VEL;
   
   creal DTDX = DT/( grid.getCellParams()[SPATCELL*SIZE_CELLPARAMS + CellParams::DX] );
   creal DTDY = DT/( grid.getCellParams()[SPATCELL*SIZE_CELLPARAMS + CellParams::DY] );
   creal DTDZ = DT/( grid.getCellParams()[SPATCELL*SIZE_CELLPARAMS + CellParams::DZ] );
   
   real avg,flux_neg,flux_pos;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Calculate the contribution coming from x-fluxes:
      flux_neg = fx[BLOCK*SIZE_FLUXS + trIndex(i,j,k)];
      // Get positive face flux from +x neighbour, or calculate it using a bound.func.:
      if (nbrs[NbrsVel::STATE] & NbrsVel::X_POS_BND == NbrsVel::X_POS_BND) {
	 flux_pos = 0.0;
      } else {
	 flux_pos = fx[nbrs[NbrsVel::XPOS]*SIZE_FLUXS + trIndex(i,j,k)];
      }
      avg = (flux_neg - flux_pos)*DTDX;
     
      // Calculate the contribution coming from y-fluxes:
      flux_neg = fy[BLOCK*SIZE_FLUXS + trIndex(i,j,k)];
      // Get positive face flux from +y neighbour, or calculate it using a bound.func.:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Y_POS_BND == NbrsVel::Y_POS_BND) {
	 flux_pos = 0.0;
      } else {
	 flux_pos = fy[nbrs[NbrsVel::YPOS]*SIZE_FLUXS + trIndex(i,j,k)];
      }
      avg += (flux_neg - flux_pos)*DTDY;
      
      // Calculate the contribution coming from z-fluxes:
      flux_neg = fz[BLOCK*SIZE_FLUXS + trIndex(i,j,k)];
      // Get positive face flux from +z neighbour, or calculate it using a bound.func.:
      if (nbrs[NbrsVel::STATE] & NbrsVel::Z_POS_BND == NbrsVel::Z_POS_BND) {
	 flux_pos = 0.0;
      } else {
	 flux_pos = fz[nbrs[NbrsVel::ZPOS]*SIZE_FLUXS + trIndex(i,j,k)];
      }
      avg += (flux_neg - flux_pos)*DTDZ;

      // Store new volume average:
      avgs[BLOCK*SIZE_VELBLOCK + trIndex(i,j,k)] += avg;
   }
}

void cpu_propagateSpat(SpatialCell& cell,cuint& BLOCK,std::vector<SpatialCell*>& nbrPtrs,creal& DT) {
   real* const avgs = cell.cpu_avgs + BLOCK*SIZE_VELBLOCK;
   creal* const fx = cell.cpu_fx + BLOCK*SIZE_FLUXS;
   creal* const fy = cell.cpu_fy + BLOCK*SIZE_FLUXS;
   creal* const fz = cell.cpu_fz + BLOCK*SIZE_FLUXS;
   
   creal DTDX = DT / cell.cpu_cellParams[CellParams::DX];
   creal DTDY = DT / cell.cpu_cellParams[CellParams::DY];
   creal DTDZ = DT / cell.cpu_cellParams[CellParams::DZ];
   
   real avg,flux_neg,flux_pos;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      // Calculate the contribution coming from x-fluxes:
      flux_neg = fx[trIndex(i,j,k)];
      // Get positive face flux from +x neighbour, or calculate it using a bound.func.:
      if (nbrPtrs[1] == NULL) {
	 flux_pos = 0.0; // Boundary value!
      } else {
	 flux_pos = nbrPtrs[1]->cpu_fx[BLOCK*SIZE_FLUXS + trIndex(i,j,k)];
      }
      avg = (flux_neg - flux_pos)*DTDX;
      
      // Calculate the contribution coming from y-fluxes:
      flux_neg = fy[trIndex(i,j,k)];
      // Get positive face flux from +y neighbour, or calculate it using a bound.func.:
      if (nbrPtrs[3] == NULL) {
	 flux_pos = 0.0; // Boundary value!
      } else {
	 flux_pos = nbrPtrs[3]->cpu_fy[BLOCK*SIZE_FLUXS + trIndex(i,j,k)];
      }
      avg += (flux_neg - flux_pos)*DTDY;
      
      // Calculate the contribution coming from z-fluxes:
      flux_neg = fz[trIndex(i,j,k)];
      // Get positive face flux from +z neighbour, or calculate it using a bound.func.:
      if (nbrPtrs[5] == NULL) {
	 flux_pos = 0.0; // Boundary value!
      } else {
	 flux_pos = nbrPtrs[5]->cpu_fz[BLOCK*SIZE_FLUXS + trIndex(i,j,k)];
      }
      avg += (flux_neg - flux_pos)*DTDZ;

      // Store new volume average:
      avgs[trIndex(i,j,k)] += avg;
   }
}

bool cpu_translation1(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT) {
   #pragma omp parallel for
   for (uint block=cell.velBlockIndex; block<cell.velBlockIndex+cell.N_blocks; ++block) {
      cpu_calcSpatDerivs(cellIndex,block,grid);
   }
   return true;
}

bool cpu_translation2(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT) {
   #pragma omp parallel for
   for (uint block=cell.velBlockIndex; block<cell.velBlockIndex+cell.N_blocks; ++block) {
      cpu_calcSpatFluxesX(cellIndex,block,grid);
      cpu_calcSpatFluxesY(cellIndex,block,grid);
      cpu_calcSpatFluxesZ(cellIndex,block,grid);
   }
   return true;
}

bool cpu_translation3(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT) {
   #pragma omp parallel for
   for (uint block=cell.velBlockIndex; block<cell.velBlockIndex+cell.N_blocks; ++block) {
      cpu_propagateSpat(cellIndex,block,grid,DT);
   }
   return true;
}

bool cpu_translation1(SpatialCell& cell,vector<SpatialCell*>& spatNbrs) {
   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      cpu_calcSpatDerivs(cell,block,spatNbrs);
   }
   return true;
}

bool cpu_translation2(SpatialCell& cell,vector<SpatialCell*>& spatNbrs) {
   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      cpu_calcSpatFluxesX(cell,block,spatNbrs);
      cpu_calcSpatFluxesY(cell,block,spatNbrs);
      cpu_calcSpatFluxesZ(cell,block,spatNbrs);
   }
   return true;
}

bool cpu_translation3(SpatialCell& cell,vector<SpatialCell*>& spatNbrs) {
   creal DT = Parameters::dt;
   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      cpu_propagateSpat(cell,block,spatNbrs,DT);
   }
   return true;
}
      
