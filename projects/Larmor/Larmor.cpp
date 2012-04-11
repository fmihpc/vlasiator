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

typedef larmorParameters LarmP;
Real LarmP::BX0 = NAN;
Real LarmP::BY0 = NAN;
Real LarmP::BZ0 = NAN;
Real LarmP::VX0 = NAN;
Real LarmP::VY0 = NAN;
Real LarmP::VZ0 = NAN;
Real LarmP::DENSITY = NAN;
Real LarmP::TEMPERATURE = NAN;
Real LarmP::magPertAmp = NAN;
Real LarmP::densityPertAmp = NAN;
Real LarmP::velocityPertAmp = NAN;
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
   RP::add("Larmor.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
   RP::add("Larmor.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("Larmor.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
   RP::add("Larmor.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Larmor.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
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
   RP::get("Larmor.magPertAmp", LarmP::magPertAmp);
   RP::get("Larmor.densityPertAmp", LarmP::densityPertAmp);
   RP::get("Larmor.velocityPertAmp", LarmP::velocityPertAmp);
   RP::get("Larmor.nSpaceSamples", LarmP::nSpaceSamples);
   RP::get("Larmor.nVelocitySamples", LarmP::nVelocitySamples);
   RP::get("Larmor.maxwCutoff", LarmP::maxwCutoff);
   RP::get("Larmor.Scale_x", LarmP::SCA_X);
   RP::get("Larmor.Scale_y", LarmP::SCA_Y);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}


Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
   
   return exp(- mass * ((vx-LarmP::VX0)*(vx-LarmP::VX0) + (vy-LarmP::VY0)*(vy-LarmP::VY0)+ (vz-LarmP::VZ0)*(vz-LarmP::VZ0)) / (2.0 * k * LarmP::TEMPERATURE))*
      exp(-pow(x-Parameters::xmax/2.5, 2.0)/pow(LarmP::SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/2.0, 2.0)/pow(LarmP::SCA_Y, 2.0));
}

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

   static int rndRho = 0;
   static int rndVel[3] = {0};
   int cellID = (int) (x / dx) +
                (int) (y / dy) * Parameters::xcells_ini +
                (int) (z / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   srand(cellID);
   
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
      
void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   int cellID = (int) (x / dx) +
   (int) (y / dy) * Parameters::xcells_ini +
   (int) (z / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   srand(cellID);
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = LarmP::BZ0;
   cellParams[CellParams::BXVOL   ] = 0.0;
   cellParams[CellParams::BYVOL   ] = 0.0;
   cellParams[CellParams::BZVOL   ] = LarmP::BZ0;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

