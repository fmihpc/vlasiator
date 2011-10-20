
/*
This file is part of Vlasiator.

Copyright 2011 Finnish Meteorological Institute

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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"
#include "parameters.h"

using namespace std;

typedef alfvenParameters AP;
Real AP::B0 = NAN;
Real AP::DENSITY = NAN;
Real AP::ALPHA = NAN;
Real AP::COS_ALPHA = NAN;
Real AP::SIN_ALPHA = NAN;
Real AP::WAVELENGTH = NAN;
Real AP::TEMPERATURE = NAN;
Real AP::A_VEL = NAN;
Real AP::A_MAG = NAN;
uint AP::nSpaceSamples = 0;
uint AP::nVelocitySamples = 0;

bool initializeProject(void) {
   Real Bx_guiding, By_guiding, Bz_guiding;
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
   RP::parse();
   RP::get("Alfven.B0", AP::B0);
   RP::get("Alfven.Bx_guiding", Bx_guiding);
   RP::get("Alfven.By_guiding", By_guiding);
   RP::get("Alfven.Bz_guiding", Bz_guiding);
   RP::get("Alfven.rho", AP::DENSITY);
   RP::get("Alfven.Wavelength", AP::WAVELENGTH);
   RP::get("Alfven.Temperature", AP::TEMPERATURE);
   RP::get("Alfven.A_mag", AP::A_MAG);
   RP::get("Alfven.A_vel", AP::A_VEL);
   RP::get("Alfven.nSpaceSamples", AP::nSpaceSamples);
   RP::get("Alfven.nVelocitySamples", AP::nVelocitySamples);
   
   Real norm = sqrt(Bx_guiding*Bx_guiding + By_guiding*By_guiding + Bz_guiding*Bz_guiding);
   Bx_guiding /= norm;
   By_guiding /= norm;
   By_guiding /= norm;
   AP::ALPHA = atan(By_guiding/Bx_guiding);
   AP::COS_ALPHA=cos(AP::ALPHA);
   AP::SIN_ALPHA=sin(AP::ALPHA);
   
   return true;
} 

bool cellParametersChanged(creal& t) {return false;}

/*Real calcPhaseSpaceDensity(creal& z,creal& x,creal& y,creal& dz,creal& dx,creal& dy,
			   creal& vz,creal& vx,creal& vy,creal& dvz,creal& dvx,creal& dvy) {*/
Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   static creal mass = 1.67262171e-27; // m_p in kg
   static creal k = 1.3806505e-23; // Boltzmann
   static creal mu0 = 1.25663706144e-6; // mu_0
//   creal q = 1.60217653e-19; // q_i
   static creal ALFVEN_VEL = AP::B0 / sqrt(mu0 * AP::DENSITY * mass);
   static creal densityFactor=AP::DENSITY * pow(mass / (2.0 * M_PI * k * AP::TEMPERATURE), 1.5);
   
   creal ksi = (x * AP::COS_ALPHA + y * AP::SIN_ALPHA) / AP::WAVELENGTH;
   creal Vx = AP::A_VEL * ALFVEN_VEL * AP::SIN_ALPHA * sin(2.0 * M_PI * ksi);
   creal Vy = - AP::A_VEL * ALFVEN_VEL * AP::COS_ALPHA * sin(2.0 * M_PI * ksi);
   creal Vz = - AP::A_VEL * ALFVEN_VEL * cos(2.0 * M_PI * ksi);

   creal den = densityFactor *
   exp(- mass * ((vx - Vx)*(vx-Vx) + (vy - Vy)*(vy-Vy) + (vz - Vz)*(vz-Vz)) / (2.0 * k * AP::TEMPERATURE));
  return den;
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal d_x = dx / (AP::nSpaceSamples-1);
   creal d_y = dy / (AP::nSpaceSamples-1);
   creal d_z = dz / (AP::nSpaceSamples-1);
   creal d_vx = dvx / (AP::nVelocitySamples-1);
   creal d_vy = dvy / (AP::nVelocitySamples-1);
   creal d_vz = dvz / (AP::nVelocitySamples-1);
   Real avg = 0.0;
#pragma omp parallel for collapse(6) reduction(+:avg)
   for (uint i=0; i<AP::nSpaceSamples; ++i)
       for (uint j=0; j<AP::nSpaceSamples; ++j)
	 for (uint k=0; k<AP::nSpaceSamples; ++k)
	    for (uint vi=0; vi<AP::nVelocitySamples; ++vi)
	       for (uint vj=0; vj<AP::nVelocitySamples; ++vj)
		  for (uint vk=0; vk<AP::nVelocitySamples; ++vk)
		  {
		     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
		  }
   return avg / AP::nSpaceSamples*AP::nSpaceSamples*AP::nSpaceSamples / AP::nVelocitySamples*AP::nVelocitySamples*AP::nVelocitySamples;

}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   
   Real dBxavg, dByavg, dBzavg, ksi;
   dBxavg = dByavg = dBzavg = 0.0;
   Real d_x = dx / (AP::nSpaceSamples - 1);
   Real d_y = dy / (AP::nSpaceSamples - 1);
   
   for (uint i=0; i<AP::nSpaceSamples; ++i)
      for (uint j=0; j<AP::nSpaceSamples; ++j)
	 for (uint k=0; k<AP::nSpaceSamples; ++k) {
	    ksi = ((x + i * d_x)  * AP::COS_ALPHA + (y + j * d_y) * AP::SIN_ALPHA) / AP::WAVELENGTH;
	    dBxavg += sin(2.0 * M_PI * ksi);
	    dByavg += sin(2.0 * M_PI * ksi);
	    dBzavg += cos(2.0 * M_PI * ksi);
	 }
   cuint nPts = AP::nSpaceSamples*AP::nSpaceSamples*AP::nSpaceSamples;
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = AP::B0 * AP::COS_ALPHA - AP::A_MAG * AP::B0 * AP::SIN_ALPHA * dBxavg / nPts;
   cellParams[CellParams::BY   ] = AP::B0 * AP::SIN_ALPHA + AP::A_MAG * AP::B0 * AP::COS_ALPHA * dByavg / nPts;
   cellParams[CellParams::BZ   ] = AP::B0 * AP::A_MAG * dBzavg / nPts;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
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

