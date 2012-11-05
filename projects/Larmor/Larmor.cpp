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
#include <iomanip>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

typedef larmorParameters LarmP;
Real LarmP::BX0 = NAN;
Real LarmP::BY0 = NAN;
Real LarmP::BZ0 = NAN;
Real LarmP::VX0 = NAN;
Real LarmP::VY0 = NAN;
Real LarmP::VZ0 = NAN;
Real LarmP::DENSITY = NAN;
Real LarmP::TEMPERATURE = NAN;
Real LarmP::maxwCutoff = NAN;
//uint LarmP::sectorSize = 0;
uint LarmP::nSpaceSamples = 0;
uint LarmP::nVelocitySamples = 0;
Real LarmP::SCA_X = NAN;
Real LarmP::SCA_Y = NAN;

bool initializeProject(void) {return true;}

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Larmor.BX0", "Background field value (T)", 1.0e-9);
   RP::add("Larmor.BY0", "Background field value (T)", 2.0e-9);
   RP::add("Larmor.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("Larmor.VX0", "Bulk velocity in x", 0.0);
   RP::add("Larmor.VY0", "Bulk velocity in y", 0.0);
   RP::add("Larmor.VZ0", "Bulk velocuty in z", 0.0);
   RP::add("Larmor.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Larmor.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Larmor.nSpaceSamples", "Number of sampling points per spatial dimension", 1);
   RP::add("Larmor.nVelocitySamples", "Number of sampling points per velocity dimension", 1);
   RP::add("Larmor.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   RP::add("Larmor.Scale_x", "Scale length in x (m)", 2.0e6);
   RP::add("Larmor.Scale_y", "Scale length in y (m)", 2.0e6);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("Larmor.BX0", LarmP::BX0);
   RP::get("Larmor.BY0", LarmP::BY0);
   RP::get("Larmor.BZ0", LarmP::BZ0);
   RP::get("Larmor.VX0", LarmP::VX0);
   RP::get("Larmor.VY0", LarmP::VY0);
   RP::get("Larmor.VZ0", LarmP::VZ0);
   RP::get("Larmor.rho", LarmP::DENSITY);
   RP::get("Larmor.Temperature", LarmP::TEMPERATURE);
   RP::get("Larmor.nSpaceSamples", LarmP::nSpaceSamples);
   RP::get("Larmor.nVelocitySamples", LarmP::nVelocitySamples);
   RP::get("Larmor.maxwCutoff", LarmP::maxwCutoff);
   RP::get("Larmor.Scale_x", LarmP::SCA_X);
   RP::get("Larmor.Scale_y", LarmP::SCA_Y);
   return true;
}

Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
   
   return exp(- mass * ((vx-LarmP::VX0)*(vx-LarmP::VX0) + (vy-LarmP::VY0)*(vy-LarmP::VY0)+ (vz-LarmP::VZ0)*(vz-LarmP::VZ0)) / (2.0 * k * LarmP::TEMPERATURE))*
      exp(-pow(x-Parameters::xmax/2.5, 2.0)/pow(LarmP::SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/2.0, 2.0)/pow(LarmP::SCA_Y, 2.0));
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
   if(vx < Parameters::vxmin + 0.5 * dvx ||
      vy < Parameters::vymin + 0.5 * dvy ||
      vz < Parameters::vzmin + 0.5 * dvz ||
      vx > Parameters::vxmax - 1.5 * dvx ||
      vy > Parameters::vymax - 1.5 * dvy ||
      vz > Parameters::vzmax - 1.5 * dvz
   ) return 0.0;
   
   creal mass = Parameters::m;
   creal q = Parameters::q;
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0

   creal d_x = dx / (LarmP::nSpaceSamples-1);
   creal d_y = dy / (LarmP::nSpaceSamples-1);
   creal d_z = dz / (LarmP::nSpaceSamples-1);
   creal d_vx = dvx / (LarmP::nVelocitySamples-1);
   creal d_vy = dvy / (LarmP::nVelocitySamples-1);
   creal d_vz = dvz / (LarmP::nVelocitySamples-1);
   Real avg = 0.0;
   
   for (uint i=0; i<LarmP::nSpaceSamples; ++i)
      for (uint j=0; j<LarmP::nSpaceSamples; ++j)
         for (uint k=0; k<LarmP::nSpaceSamples; ++k)
            for (uint vi=0; vi<LarmP::nVelocitySamples; ++vi)
               for (uint vj=0; vj<LarmP::nVelocitySamples; ++vj)
                  for (uint vk=0; vk<LarmP::nVelocitySamples; ++vk)
                  {
                     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
                  }
   
   creal result = avg *LarmP::DENSITY * pow(mass / (2.0 * M_PI * k * LarmP::TEMPERATURE), 1.5) /
      (LarmP::nSpaceSamples*LarmP::nSpaceSamples*LarmP::nSpaceSamples) / 
      (LarmP::nVelocitySamples*LarmP::nVelocitySamples*LarmP::nVelocitySamples);
   
				  
   if(result < LarmP::maxwCutoff) {
      return 0.0;
   } else {
      return result;
   }
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
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BGBX   ] = 0.0;
   cellParams[CellParams::BGBY   ] = 0.0;
   cellParams[CellParams::BGBZ   ] = LarmP::BZ0;
   cellParams[CellParams::PERBX   ] = 0.0;
   cellParams[CellParams::PERBY   ] = 0.0;
   cellParams[CellParams::PERBZ   ] = 0.0;
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
