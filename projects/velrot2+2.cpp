#include <cstdlib>
#include <iostream>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"
#include "parameters.h"

using namespace std;

bool initializeProject(void) {
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real getDistribValue(creal& OMEGA,creal& x,creal& y,creal& vx,creal& vy,creal& dvx,creal& dvy) {
   
   creal NORM  = 1.0;  // Norm. constant
   creal R0    = 0.6;  // Center point of distribution in r
   creal SIGMA = 0.15; // Sigma of r-distribution
   creal SIGMA2 = SIGMA*SIGMA;
   
   // Calculate r,phi coords. corresponding to x,y. Check that 
   // (r,phi) has non-zero distribution function value:
   creal r = sqrt(x*x+y*y);
   Real phi = acos(x/r);
   if (y < 0.0) phi = 2.0*M_PI - phi;
   
   creal PHI_CENT = 320 / 180.0*M_PI;
   creal PHI_WIDTH = 60.0 / 180.0*M_PI;
   creal PHI0 = PHI_CENT - 0.5*PHI_WIDTH;
   creal PHI1 = PHI_CENT + 0.5*PHI_WIDTH;
   if (phi < PHI0 || phi > PHI1) return 0.0;
   
   // Calculate v,theta coordinates for (r,phi):
   creal v = OMEGA*r;
   Real theta = phi - M_PI/2.0;
   while (theta < 0.0) theta += 2.0*M_PI;
   
   // Calculate vx',vy' corresponding to (v,theta), and 
   // check that vx_dot,vy_dot is inside the given velocity cell:
   creal vx_dot = v*cos(theta);
   creal vy_dot = v*sin(theta);   
   if (vx_dot < vx || vx_dot > vx+dvx) return 0.0;
   if (vy_dot < vy || vy_dot > vy+dvy) return 0.0;
   
   creal arg = (r-R0)*(r-R0)/SIGMA2;
   return NORM*exp(-1.0*arg);
   
   /*
   creal V = sqrt(vx*vx+vy*vy);
   if (V < 0.1) return 0.0;
   
   creal NORM = 1.0;
   creal V0 = 0.6;
   creal SIGMA = 0.15;
   
   creal THETA_CENT = 320 / 180.0*M_PI;
   creal THETA_WIDTH = 60.0 / 180.0*M_PI;
   creal THETA0 = THETA_CENT - 0.5*THETA_WIDTH;
   creal THETA1 = THETA_CENT + 0.5*THETA_WIDTH;
   
   Real theta = acos(vx/V);
   if (vy < 0.0) theta = 2.0*M_PI-theta;
   //if (V < V0-DELTA || V > V0+DELTA) return 0.0;
   if (theta < THETA0 || theta > THETA1) return 0.0;
   
   creal SIGMA2 = SIGMA*SIGMA;
   creal arg = (V-V0)*(V-V0)/SIGMA2;
   
   return NORM * exp(-1.0*arg)/V;
   */
}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   // Calculate (r,phi) corresponding to (x,y), and check that 
   // the coordinates are within the wedge containing nonzero distrib.f. values:
   creal r = sqrt(x*x+y*y);
   Real phi = acos(x/r);
   creal EPSPHI = acos(dx/r);
   if (y < 0.0) phi = 2.0*M_PI  - phi;
   
   creal PHI_CENT = 320 / 180.0*M_PI;
   creal PHI_WIDTH = 60.0 / 180.0*M_PI;
   creal PHI0 = PHI_CENT - 0.5*PHI_WIDTH;
   creal PHI1 = PHI_CENT + 0.5*PHI_WIDTH;
   if (phi < PHI0-EPSPHI || phi > PHI1+EPSPHI) return 0.0;   
   
   // Calculate (v,theta) corresponding to (vx,vy), and check that 
   // the coordinates are within the wedge of nonzero distrib. f.:
   creal v = sqrt(vx*vx+vy*vy);
   Real theta = acos(vx/v);
   if (vy < 0.0) theta = 2.0*M_PI - theta;
   //if (vz < 0.0) theta = 2.0*M_PI - theta;
   creal EPSTHETA = acos(dvx/v);
   
   Real THETA0 = PHI0 - M_PI/2.0;
   Real THETA1 = PHI1 - M_PI/2.0;
   if (THETA0 < 0.0) THETA0 += 2.0*M_PI;
   if (THETA1 < 0.0) THETA1 += 2.0*M_PI;
   if (theta < THETA0-EPSTHETA || theta > THETA1+EPSTHETA) return 0.0;
   
   // Sample the 4D distribution function.
   creal OMEGA = 1.0;
   cuint N_samples = 100;

   creal d_x = dx / (N_samples+1);
   creal d_y = dy / (N_samples+1);
   creal d_z = dz / (N_samples+1);
   Real avg = 0.0;
   for (uint i=0; i<N_samples; ++i) for (uint j=0; j<N_samples; ++j) {
      avg += getDistribValue(OMEGA,x+i*d_x,y+j*d_y,vx,vy,dvx,dvy);
      //avg += getDistribValue(OMEGA,x+i*d_x,z+j*d_z,vx,vz,dvx,dvz);
   }
   return avg / N_samples/N_samples;

}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = 1.0;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
#else
void calcSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
#endif
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->cpu_cellParams, t);
   }
}

