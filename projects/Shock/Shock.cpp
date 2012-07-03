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

typedef gradbParameters GradBP;
Real GradBP::BX0 = NAN;
Real GradBP::BY0 = NAN;
Real GradBP::BZ0 = NAN;
Real GradBP::EX0 = NAN;
Real GradBP::VX0 = NAN;
Real GradBP::VY0 = NAN;
Real GradBP::VZ0 = NAN;
Real GradBP::DENSITY = NAN;
Real GradBP::TEMPERATURE = NAN;
Real GradBP::magPertAmp = NAN;
Real GradBP::densityPertAmp = NAN;
Real GradBP::velocityPertAmp = NAN;
Real GradBP::maxwCutoff = NAN;
//uint GradBP::sectorSize = 0;
uint GradBP::nSpaceSamples = 0;
uint GradBP::nVelocitySamples = 0;
Real GradBP::SCA_X = NAN;
Real GradBP::SCA_Y = NAN;
Real GradBP::Sharp_Y = NAN;

bool initializeProject(void) {return true;}

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("GradB.BX0", "Background field value (T)", 1.0e-9);
   RP::add("GradB.BY0", "Background field value (T)", 2.0e-9);
   RP::add("GradB.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("GradB.EX0", "Background electric field", 0.0);
   RP::add("GradB.VX0", "Bulk velocity in x", 0.0);
   RP::add("GradB.VY0", "Bulk velocity in y", 0.0);
   RP::add("GradB.VZ0", "Bulk velocuty in z", 0.0);
   RP::add("GradB.rho", "Number density (m^-3)", 1.0e7);
   RP::add("GradB.Temperature", "Temperature (K)", 2.0e6);
   RP::add("GradB.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
   RP::add("GradB.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("GradB.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
   RP::add("GradB.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("GradB.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   RP::add("GradB.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   RP::add("GradB.Scale_x", "Scale length in x (m)", 2.0e6);
   RP::add("GradB.Scale_y", "Scale length in y (m)", 2.0e6);
   RP::add("GradB.Sharp_Y", "Sharpness of tannh", 0.1);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("GradB.BX0", GradBP::BX0);
   RP::get("GradB.BY0", GradBP::BY0);
   RP::get("GradB.BZ0", GradBP::BZ0);
   RP::get("GradB.EX0", GradBP::EX0);
   RP::get("GradB.VX0", GradBP::VX0);
   RP::get("GradB.VY0", GradBP::VY0);
   RP::get("GradB.VZ0", GradBP::VZ0);
   RP::get("GradB.rho", GradBP::DENSITY);
   RP::get("GradB.Temperature", GradBP::TEMPERATURE);
   RP::get("GradB.magPertAmp", GradBP::magPertAmp);
   RP::get("GradB.densityPertAmp", GradBP::densityPertAmp);
   RP::get("GradB.velocityPertAmp", GradBP::velocityPertAmp);
   RP::get("GradB.nSpaceSamples", GradBP::nSpaceSamples);
   RP::get("GradB.nVelocitySamples", GradBP::nVelocitySamples);
   RP::get("GradB.maxwCutoff", GradBP::maxwCutoff);
   RP::get("GradB.Scale_x", GradBP::SCA_X);
   RP::get("GradB.Scale_y", GradBP::SCA_Y);
   RP::get("GradB.Sharp_Y", GradBP::Sharp_Y);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}


Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
	return exp(- mass * ((vx-GradBP::VX0)*(vx-GradBP::VX0) + (vy-GradBP::VY0)*(vy-GradBP::VY0)+ (vz-GradBP::VZ0)*(vz-GradBP::VZ0)) / (2.0 * k * GradBP::TEMPERATURE));
	//*exp(-pow(x-Parameters::xmax/2.0, 2.0)/pow(GradBP::SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/4.0, 2.0)/pow(GradBP::SCA_Y, 2.0));
}

//Real getGradValue(creal& x, creal& y, creal& z) {
//	return x/Parameters::xmax;
//}

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

   creal d_x = dx / (GradBP::nSpaceSamples-1);
   creal d_y = dy / (GradBP::nSpaceSamples-1);
   creal d_z = dz / (GradBP::nSpaceSamples-1);
   creal d_vx = dvx / (GradBP::nVelocitySamples-1);
   creal d_vy = dvy / (GradBP::nVelocitySamples-1);
   creal d_vz = dvz / (GradBP::nVelocitySamples-1);
   Real avg = 0.0;
   
   for (uint i=0; i<GradBP::nSpaceSamples; ++i)
    for (uint j=0; j<GradBP::nSpaceSamples; ++j)
	 for (uint k=0; k<GradBP::nSpaceSamples; ++k)      
     for (uint vi=0; vi<GradBP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<GradBP::nVelocitySamples; ++vj)
	  for (uint vk=0; vk<GradBP::nVelocitySamples; ++vk)
         {
	     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
  	     }
   
   creal result = avg *GradBP::DENSITY * pow(mass / (2.0 * M_PI * k * GradBP::TEMPERATURE), 1.5) /
                    (GradBP::nSpaceSamples*GradBP::nSpaceSamples*GradBP::nSpaceSamples) / 
//   	            (Parameters::vzmax - Parameters::vzmin) / 
                  (GradBP::nVelocitySamples*GradBP::nVelocitySamples*GradBP::nVelocitySamples);
				  
				  
   if(result < GradBP::maxwCutoff) {
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
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = GradBP::BZ0*(3.0 + 2.0*tanh((y - Parameters::ymax/2.0)/(GradBP::Sharp_Y*Parameters::ymax)));
   // cellParams[CellParams::BZ   ] = GradBP::BZ0*(1.0 + GradBP::DelB*y/Parameters::ymax);
   
   cellParams[CellParams::EXVOL   ] = GradBP::EX0;
   //cellParams[CellParams::EXVOL   ] = -1.0e-5;
   
   
   cellParams[CellParams::BXVOL   ] = 0.0;
   cellParams[CellParams::BYVOL   ] = 0.0;
   cellParams[CellParams::BZVOL   ] = GradBP::BZ0*(3.0 + 2.0*tanh((y - 3.0*Parameters::ymax/4.0)/(GradBP::Sharp_Y*Parameters::ymax)));
   //cellParams[CellParams::BZVOL   ] = GradBP::BZ0*(3.0 + 2.0*tanh((y - (x-0.5*Parameters::xmax)/10.0 - 3.0*Parameters::ymax/4.0) / 
   // (GradBP::Sharp_Y*Parameters::ymax)));
   //cellParams[CellParams::BZVOL   ] = GradBP::BZ0*(1.0 + GradBP::DelB*y/Parameters::ymax);
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

