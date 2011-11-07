/*
 This file is part of Vl*asiator.
 
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

typedef diffusionParameters DiffP;
Real DiffP::B0 = NAN;
Real DiffP::DENSITY = NAN;
Real DiffP::SCA_X = NAN;
Real DiffP::SCA_Y = NAN;
Real DiffP::TEMPERATURE = NAN;
uint DiffP::nSpaceSamples = 0;
uint DiffP::nVelocitySamples = 0;

bool initializeProject(void) {
   typedef Readparameters RP;
   RP::add("Diffusion.B0", "Background field value (T)", 1.0e-9);
   RP::add("Diffusion.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Diffusion.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Diffusion.Scale_x", "Scale length in x (m)", 100000.0);
   RP::add("Diffusion.Scale_y", "Scale length in y (m)", 100000.0);
   RP::add("Diffusion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Diffusion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   RP::parse();
   RP::get("Diffusion.B0", DiffP::B0);
   RP::get("Diffusion.rho", DiffP::DENSITY);
   RP::get("Diffusion.Temperature", DiffP::TEMPERATURE);
   RP::get("Diffusion.Scale_x", DiffP::SCA_X);
   RP::get("Diffusion.Scale_y", DiffP::SCA_Y);
   RP::get("Diffusion.nSpaceSamples", DiffP::nSpaceSamples);
   RP::get("Diffusion.nVelocitySamples", DiffP::nVelocitySamples);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}


Real getDistribValue(creal& x,creal& y,creal& z,
		     creal& vx,creal& vy,creal& vz) {
  
  creal mass = 1.67262171e-27; // m_p in kg
  creal k = 1.3806505e-23; // Boltzmann
  creal mu0 = 1.25663706144e-6; // mu_0
  creal q = 1.60217653e-19; // q_i
  
  return DiffP::DENSITY * pow(mass / (2.0 * M_PI * k * DiffP::TEMPERATURE), 1.5) * (
     5.0 * exp(- (pow(x, 2.0) / pow(DiffP::SCA_X, 2.0) +  pow(y, 2.0) / pow(DiffP::SCA_Y, 2.0))) * 
     exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * k * DiffP::TEMPERATURE))
    +
    exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * k * DiffP::TEMPERATURE)));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal d_x = dx / (DiffP::nSpaceSamples-1);
   creal d_y = dy / (DiffP::nSpaceSamples-1);
   creal d_z = dz / (DiffP::nSpaceSamples-1);
   creal d_vx = dvx / (DiffP::nVelocitySamples-1);
   creal d_vy = dvy / (DiffP::nVelocitySamples-1);
   creal d_vz = dvz / (DiffP::nVelocitySamples-1);
   Real avg = 0.0;
   #pragma omp parallel for collapse(6) reduction(+:avg)
   for (uint i=0; i<DiffP::nSpaceSamples; ++i)
      for (uint j=0; j<DiffP::nSpaceSamples; ++j)
	 for (uint k=0; k<DiffP::nSpaceSamples; ++k)
	    for (uint vi=0; vi<DiffP::nVelocitySamples; ++vi)
	       for (uint vj=0; vj<DiffP::nVelocitySamples; ++vj)
		  for (uint vk=0; vk<DiffP::nVelocitySamples; ++vk)
		  {
		     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
		  }
   return avg / (DiffP::nSpaceSamples*DiffP::nSpaceSamples*DiffP::nSpaceSamples) / (DiffP::nVelocitySamples*DiffP::nVelocitySamples*DiffP::nVelocitySamples);
	       
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
  creal x = cellParams[CellParams::XCRD];
  creal dx = cellParams[CellParams::DX];
  
  cellParams[CellParams::EX   ] = 0.0;
  cellParams[CellParams::EY   ] = 0.0;
  cellParams[CellParams::EZ   ] = 0.0;
  cellParams[CellParams::BX   ] = 0.0;
  cellParams[CellParams::BY   ] = 0.0;
  cellParams[CellParams::BZ   ] = DiffP::B0;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->cpu_cellParams, t);
   }
}

