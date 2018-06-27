/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

using namespace std;

bool initializeProject(void) {return true;}
bool addProjectParameters(){return true;}
bool getProjectParameters(){return true;}

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

/** Integrate the distribution function over the given six-dimensional phase-space cell.
 * @param x Starting value of the x-coordinate of the cell.
 * @param y Starting value of the y-coordinate of the cell.
 * @param z Starting value of the z-coordinate of the cell.
 * @param dx The size of the cell in x-direction.
 * @param dy The size of the cell in y-direction.
 * @param dz The size of the cell in z-direction.
 * @param vx Starting value of the vx-coordinate of the cell.
 * @param vy Starting value of the vy-coordinate of the cell.
 * @param vz Starting value of the vz-coordinate of the cell.
 * @param dvx The size of the cell in vx-direction.
 * @param dvy The size of the cell in vy-direction.
 * @param dvz The size of the cell in vz-direction.
 * @return The volume average of the distribution function in the given phase space cell.
 * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
 */
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

/** Calculate parameters for the given spatial cell at the given time.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
 * 
 * The following array indices contain the coordinates of the "lower left corner" of the cell: 
 * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
 * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
 * @param cellParams Array containing cell parameters.
 * @param t The current value of time. This is passed as a convenience. If you need more detailed information 
 * of the state of the simulation, you can read it from Parameters.
 */
void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::PERBX   ] = 0.0;
   cellParams[CellParams::PERBY   ] = 0.0;
   cellParams[CellParams::PERBZ   ] = 0.0;
   cellParams[CellParams::BGBX   ] = 0.0;
   cellParams[CellParams::BGBY   ] = 0.0;
   cellParams[CellParams::BGBZ   ] = 1.0;
}

void setProjectCell(SpatialCell* cell) {
   // Set up cell parameters:
   calcCellParameters(&((*cell).parameters[0]), 0.0);
   
   cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
   creal dvy_block = SpatialCell::block_dvy; //                    vy
   creal dvz_block = SpatialCell::block_dvz; //                    vz
   creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
      for (uint jv=0; jv<P::vyblocks_ini; ++jv)
         for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
            creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
            creal vy_block = P::vymin + jv*dvy_block; // vy-
            creal vz_block = P::vzmin + kv*dvz_block; // vz-
            
            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) 
               for (uint jc=0; jc<WID; ++jc) 
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vx_cell = vx_block + ic*dvx_blockCell;
                     creal vy_cell = vy_block + jc*dvy_blockCell;
                     creal vz_cell = vz_block + kc*dvz_blockCell;
                     Real average = 
                     calcPhaseSpaceDensity(cell->parameters[CellParams::XCRD],
                                           cell->parameters[CellParams::YCRD],
                                           cell->parameters[CellParams::ZCRD],
                                           cell->parameters[CellParams::DX],
                                           cell->parameters[CellParams::DY],
                                           cell->parameters[CellParams::DZ],
                                           vx_cell,vy_cell,vz_cell,
                                           dvx_blockCell,dvy_blockCell,dvz_blockCell);
                     
                     if(average!=0.0){
                        creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                        creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                        creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                        cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                     }
                  }
         }
         calculateCellVelocityMoments(cell);
         
         //let's get rid of blocks not fulfilling the criteria here to save memory.
         cell->adjustSingleCellVelocityBlocks();
}
