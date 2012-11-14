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

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"

#include "Larmor.h"

namespace projects {
    Larmor::Larmor(): Project() { }
    Larmor::~Larmor() { }


    bool Larmor::initialize(void) {return true;}

    void Larmor::addParameters() {
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
    }

    void Larmor::getParameters() {
      typedef Readparameters RP;
      RP::get("Larmor.BX0", this->BX0);
      RP::get("Larmor.BY0", this->BY0);
      RP::get("Larmor.BZ0", this->BZ0);
      RP::get("Larmor.VX0", this->VX0);
      RP::get("Larmor.VY0", this->VY0);
      RP::get("Larmor.VZ0", this->VZ0);
      RP::get("Larmor.rho", this->DENSITY);
      RP::get("Larmor.Temperature", this->TEMPERATURE);
      RP::get("Larmor.nSpaceSamples", this->nSpaceSamples);
      RP::get("Larmor.nVelocitySamples", this->nVelocitySamples);
      RP::get("Larmor.maxwCutoff", this->maxwCutoff);
      RP::get("Larmor.Scale_x", this->SCA_X);
      RP::get("Larmor.Scale_y", this->SCA_Y);
    }

    Real Larmor::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz) {
      creal k = 1.3806505e-23; // Boltzmann
      creal mass = 1.67262171e-27; // m_p in kg
      
      return exp(- mass * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0)+ (vz-this->VZ0)*(vz-this->VZ0)) / (2.0 * k * this->TEMPERATURE))*
	  exp(-pow(x-Parameters::xmax/2.5, 2.0)/pow(this->SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/2.0, 2.0)/pow(this->SCA_Y, 2.0));
    }


    Real Larmor::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
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

      creal d_x = dx / (this->nSpaceSamples-1);
      creal d_y = dy / (this->nSpaceSamples-1);
      creal d_z = dz / (this->nSpaceSamples-1);
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      
      for (uint i=0; i<this->nSpaceSamples; ++i)
	  for (uint j=0; j<this->nSpaceSamples; ++j)
	    for (uint k=0; k<this->nSpaceSamples; ++k)
		for (uint vi=0; vi<this->nVelocitySamples; ++vi)
		  for (uint vj=0; vj<this->nVelocitySamples; ++vj)
		      for (uint vk=0; vk<this->nVelocitySamples; ++vk)
		      {
			avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
		      }
      
      creal result = avg *this->DENSITY * pow(mass / (2.0 * M_PI * k * this->TEMPERATURE), 1.5) /
	  (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) / 
	  (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
      
				      
      if(result < this->maxwCutoff) {
	  return 0.0;
      } else {
	  return result;
      }
    }


    void Larmor::calcCellParameters(Real* cellParams,creal& t) {
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
      cellParams[CellParams::BGBZ   ] = this->BZ0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
    }
} //namespace projects 
  
