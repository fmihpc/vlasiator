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

typedef alfvenParameters AP;
Real AP::B0 = NAN;
Real AP::Bx_guiding = NAN;
Real AP::By_guiding = NAN;
Real AP::Bz_guiding = NAN;
Real AP::DENSITY = NAN;
Real AP::ALPHA = NAN;
Real AP::WAVELENGTH = NAN;
Real AP::TEMPERATURE = NAN;
Real AP::A_VEL = NAN;
Real AP::A_MAG = NAN;
uint AP::nSpaceSamples = 0;
uint AP::nVelocitySamples = 0;

bool initializeProject(void) {
   Real norm = sqrt(AP::Bx_guiding*AP::Bx_guiding + AP::By_guiding*AP::By_guiding + AP::Bz_guiding*AP::Bz_guiding);
   AP::Bx_guiding /= norm;
   AP::By_guiding /= norm;
   AP::By_guiding /= norm;
   AP::ALPHA = atan(AP::By_guiding/AP::Bx_guiding);
   
   return true;
} 

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Alfven.B0", "Guiding field value (T)", 1.0e-10);
   RP::add("Alfven.Bx_guiding", "Guiding field x component", 1);
   RP::add("Alfven.By_guiding", "Guiding field y component", 0);
   RP::add("Alfven.Bz_guiding", "Guiding field z component", 0);
   RP::add("Alfven.rho", "Number density (m^-3)", 1.0e8);
   RP::add("Alfven.Wavelength", "Wavelength (m)", 100000.0);
   RP::add("Alfven.Temperature", "Temperature (K)", 0.86456498092);
   RP::add("Alfven.A_mag", "Amplitude of the magnetic perturbation", 0.1);
   RP::add("Alfven.A_vel", "Amplitude of the velocity perturbation", 0.1);
   RP::add("Alfven.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Alfven.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Alfven.B0", AP::B0);
   RP::get("Alfven.Bx_guiding", AP::Bx_guiding);
   RP::get("Alfven.By_guiding", AP::By_guiding);
   RP::get("Alfven.Bz_guiding", AP::Bz_guiding);
   RP::get("Alfven.rho", AP::DENSITY);
   RP::get("Alfven.Wavelength", AP::WAVELENGTH);
   RP::get("Alfven.Temperature", AP::TEMPERATURE);
   RP::get("Alfven.A_mag", AP::A_MAG);
   RP::get("Alfven.A_vel", AP::A_VEL);
   RP::get("Alfven.nSpaceSamples", AP::nSpaceSamples);
   RP::get("Alfven.nVelocitySamples", AP::nVelocitySamples);
   
   return true;
}

/*Real calcPhaseSpaceDensity(creal& z,creal& x,creal& y,creal& dz,creal& dx,creal& dy,
			   creal& vz,creal& vx,creal& vy,creal& dvz,creal& dvx,creal& dvy) {*/
Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0
//   creal q = 1.60217653e-19; // q_i
   creal ALFVEN_VEL = AP::B0 / sqrt(mu0 * AP::DENSITY * mass);

   creal ksi = (x * cos(AP::ALPHA) + y * sin(AP::ALPHA)) / AP::WAVELENGTH;
   creal Vx = AP::A_VEL * ALFVEN_VEL * sin(AP::ALPHA) * sin(2.0 * M_PI * ksi);
   creal Vy = - AP::A_VEL * ALFVEN_VEL * cos(AP::ALPHA) * sin(2.0 * M_PI * ksi);
   creal Vz = - AP::A_VEL * ALFVEN_VEL * cos(2.0 * M_PI * ksi);
  
   creal den = AP::DENSITY * pow(mass / (2.0 * M_PI * k * AP::TEMPERATURE), 1.5) *
   exp(- mass * (pow(vx - Vx, 2.0) + pow(vy - Vy, 2.0) + pow(vz - Vz, 2.0)) / (2.0 * k * AP::TEMPERATURE));
  return den;
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
   creal d_x = dx / (AP::nSpaceSamples-1);
   creal d_y = dy / (AP::nSpaceSamples-1);
   creal d_z = dz / (AP::nSpaceSamples-1);
   creal d_vx = dvx / (AP::nVelocitySamples-1);
   creal d_vy = dvy / (AP::nVelocitySamples-1);
   creal d_vz = dvz / (AP::nVelocitySamples-1);
   Real avg = 0.0;
   for (uint i=0; i<AP::nSpaceSamples; ++i)
      for (uint j=0; j<AP::nSpaceSamples; ++j)
	 for (uint k=0; k<AP::nSpaceSamples; ++k)
	    for (uint vi=0; vi<AP::nVelocitySamples; ++vi)
	       for (uint vj=0; vj<AP::nVelocitySamples; ++vj)
		  for (uint vk=0; vk<AP::nVelocitySamples; ++vk)
		  {
		     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
		  }
   return avg / pow(AP::nSpaceSamples, 3.0) / pow(AP::nVelocitySamples, 3.0);
}

/** Calculate parameters for the given spatial cell at the given time.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ,
 * CellParams::PERBX, CellParams::PERBY, and CellParams::PERBZ.
 * CellParams::BGBX, CellParams::BGBY, and CellParams::BGBZ.
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
   
   Real dBxavg, dByavg, dBzavg;
   dBxavg = dByavg = dBzavg = 0.0;
   Real d_x = dx / (AP::nSpaceSamples - 1);
   Real d_y = dy / (AP::nSpaceSamples - 1);

   for (uint i=0; i<AP::nSpaceSamples; ++i)
      for (uint j=0; j<AP::nSpaceSamples; ++j)
	 for (uint k=0; k<AP::nSpaceSamples; ++k) {
	    Real ksi = ((x + i * d_x)  * cos(AP::ALPHA) + (y + j * d_y) * sin(AP::ALPHA)) / AP::WAVELENGTH;
	    dBxavg += sin(2.0 * M_PI * ksi);
	    dByavg += sin(2.0 * M_PI * ksi);
	    dBzavg += cos(2.0 * M_PI * ksi);
	 }
   cuint nPts = pow(AP::nSpaceSamples, 3.0);
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::PERBX   ] = 0.0;
   cellParams[CellParams::PERBY   ] = 0.0;
   cellParams[CellParams::PERBZ   ] = 0.0;
   cellParams[CellParams::BGBX   ] = AP::B0 * cos(AP::ALPHA) - AP::A_MAG * AP::B0 * sin(AP::ALPHA) * dBxavg / nPts;
   cellParams[CellParams::BGBY   ] = AP::B0 * sin(AP::ALPHA) + AP::A_MAG * AP::B0 * cos(AP::ALPHA) * dByavg / nPts;
   cellParams[CellParams::BGBZ   ] = AP::B0 * AP::A_MAG * dBzavg / nPts;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
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
