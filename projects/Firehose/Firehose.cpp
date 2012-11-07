/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
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

typedef FirehoseParameters FH;
Real FH::rho[2] = {NAN};
Real FH::Tx[2] = {NAN};
Real FH::Ty[2] = {NAN};
Real FH::Tz[2] = {NAN};
Real FH::Vx[2] = {NAN};
Real FH::Vy[2] = {NAN};
Real FH::Vz[2] = {NAN};
Real FH::Bx = NAN;
Real FH::By = NAN;
Real FH::Bz = NAN;
Real FH::lambda = NAN;
Real FH::amp = NAN;
uint FH::nSpaceSamples = 0;
uint FH::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Firehose.rho1", "Number density, first peak (m^-3)", 0.0);
   RP::add("Firehose.rho2", "Number density, second peak (m^-3)", 0.0);
   RP::add("Firehose.Tx1", "Temperature x, first peak (K)", 0.0);
   RP::add("Firehose.Tx2", "Temperature x, second peak (K)", 0.0);
   RP::add("Firehose.Ty1", "Temperature y, first peak (K)", 0.0);
   RP::add("Firehose.Ty2", "Temperature y, second peak (K)", 0.0);
   RP::add("Firehose.Tz1", "Temperature z, first peak (K)", 0.0);
   RP::add("Firehose.Tz2", "Temperature z, second peak (K)", 0.0);
   RP::add("Firehose.Vx1", "Bulk velocity x component, first peak (m/s)", 0.0);
   RP::add("Firehose.Vx2", "Bulk velocity x component, second peak (m/s)", 0.0);
   RP::add("Firehose.Vy1", "Bulk velocity y component, first peak (m/s)", 0.0);
   RP::add("Firehose.Vy2", "Bulk velocity y component, second peak (m/s)", 0.0);
   RP::add("Firehose.Vz1", "Bulk velocity z component, first peak (m/s)", 0.0);
   RP::add("Firehose.Vz2", "Bulk velocity z component, second peak (m/s)", 0.0);
   RP::add("Firehose.Bx", "Magnetic field x component (T)", 0.0);
   RP::add("Firehose.By", "Magnetic field y component (T)", 0.0);
   RP::add("Firehose.Bz", "Magnetic field z component (T)", 0.0);
   RP::add("Firehose.lambda", "Initial perturbation wavelength (m)", 0.0);
   RP::add("Firehose.amp", "Initial perturbation amplitude (m)", 0.0);
   RP::add("Firehose.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Firehose.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Firehose.rho1", FH::rho[1]);
   RP::get("Firehose.rho2", FH::rho[2]);
   RP::get("Firehose.Tx1", FH::Tx[1]);
   RP::get("Firehose.Tx2", FH::Tx[2]);
   RP::get("Firehose.Ty1", FH::Ty[1]);
   RP::get("Firehose.Ty2", FH::Ty[2]);
   RP::get("Firehose.Tz1", FH::Tz[1]);
   RP::get("Firehose.Tz2", FH::Tz[2]);
   RP::get("Firehose.Vx1", FH::Vx[1]);
   RP::get("Firehose.Vx2", FH::Vx[2]);
   RP::get("Firehose.Vy1", FH::Vy[1]);
   RP::get("Firehose.Vy2", FH::Vy[2]);
   RP::get("Firehose.Vz1", FH::Vz[1]);
   RP::get("Firehose.Vz2", FH::Vz[2]);
   RP::get("Firehose.Bx", FH::Bx);
   RP::get("Firehose.By", FH::By);
   RP::get("Firehose.Bz", FH::Bz);
   RP::get("Firehose.lambda", FH::lambda);
   RP::get("Firehose.amp", FH::amp);
   RP::get("Firehose.nSpaceSamples", FH::nSpaceSamples);
   RP::get("Firehose.nVelocitySamples", FH::nVelocitySamples);
   return true;
}

Real profile(creal top, creal bottom, creal x, creal ) {
     return top * (1.0 + FH::amp*cos(2.0*M_PI*x/FH::lambda));
}

Real getDistribValue(creal& x, creal& y, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   //  creal mu0 = 1.25663706144e-6; // mu_0
   //  creal q = 1.60217653e-19; // q_i
   //  creal gamma = 5./3.;
   
   Real Vx = profile(FH::Vx[1],FH::Vx[1], x, y);
   
   return
   FH::rho[1] * pow(mass / (2.0 * M_PI * k * FH::Tx[1]), 1.5) *
   exp(- mass * (pow(vx - Vx, 2.0) / (2.0 * k * FH::Tx[1]) + 
                 pow(vy - FH::Vy[1], 2.0) / (2.0 * k * FH::Ty[1]) + 
             pow(vz - FH::Vz[1], 2.0) / (2.0 * k * FH::Tz[1]))); 
//   FH::rho[2] * pow(mass / (2.0 * M_PI * k * FH::Tx[2]), 1.5) *
//   exp(- mass * (pow(vx - FH::Vx[2], 2.0) / (2.0 * k * FH::Tx[2]) + 
//                 pow(vy - FH::Vy[2], 2.0) / (2.0 * k * FH::Ty[2]) + 
//           pow(vz - FH::Vz[2], 2.0) / (2.0 * k * FH::Tz[2]))); 
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
Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   creal d_x = dx / (FH::nSpaceSamples-1);
   creal d_y = dy / (FH::nSpaceSamples-1);
   creal d_vx = dvx / (FH::nVelocitySamples-1);
   creal d_vy = dvy / (FH::nVelocitySamples-1);
   creal d_vz = dvz / (FH::nVelocitySamples-1);
   Real avg = 0.0;
//#pragma omp parallel for collapse(6) reduction(+:avg)
   for (uint i=0; i<FH::nSpaceSamples; ++i)
     for (uint j=0; j<FH::nSpaceSamples; ++j)
       for (uint vi=0; vi<FH::nVelocitySamples; ++vi)
        for (uint vj=0; vj<FH::nVelocitySamples; ++vj)
         for (uint vk=0; vk<FH::nVelocitySamples; ++vk)
       {
          avg += getDistribValue(x+i*d_x, y+j*d_y, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
       }
   return avg / pow(FH::nSpaceSamples, 2.0) /  pow(FH::nVelocitySamples, 3.0);
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
   cellParams[CellParams::PERBX   ] = FH::Bx;
   cellParams[CellParams::PERBY   ] = FH::By;
   cellParams[CellParams::PERBZ   ] = FH::Bz;
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

