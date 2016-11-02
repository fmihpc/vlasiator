/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "Diffusion.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Diffusion::Diffusion(): Project() { }
   Diffusion::~Diffusion() { }
   
   bool Diffusion::initialize(void) {
      return Project::initialize();
   }
   
   void Diffusion::addParameters() {
      typedef Readparameters RP;
      RP::add("Diffusion.B0", "Background field value (T)", 1.0e-9);
      RP::add("Diffusion.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Diffusion.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Diffusion.Scale_x", "Scale length in x (m)", 100000.0);
      RP::add("Diffusion.Scale_y", "Scale length in y (m)", 100000.0);
      RP::add("Diffusion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Diffusion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }

   void Diffusion::getParameters() {
      Project::getParameters();

      typedef Readparameters RP;
      RP::get("Diffusion.B0", this->B0);
      RP::get("Diffusion.rho", this->DENSITY);
      RP::get("Diffusion.Temperature", this->TEMPERATURE);
      RP::get("Diffusion.Scale_x", this->SCA_X);
      RP::get("Diffusion.Scale_y", this->SCA_Y);
      RP::get("Diffusion.nSpaceSamples", this->nSpaceSamples);
      RP::get("Diffusion.nVelocitySamples", this->nVelocitySamples);
   }

   Real Diffusion::getDistribValue(
      creal& x,creal& y,creal& z,
      creal& vx,creal& vy,creal& vz
   ) {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      return this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) * (
         5.0 * exp(- (pow(x, 2.0) / pow(this->SCA_X, 2.0) +  pow(y, 2.0) / pow(this->SCA_Y, 2.0))) * 
         exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * kb * this->TEMPERATURE))
         +
         exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * kb * this->TEMPERATURE)));
   }
   
   Real Diffusion::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) {
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
      return avg / (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) / (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
   }
   
   void Diffusion::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;

   }

   void Diffusion::setCellBackgroundField(SpatialCell* cell) {
      ConstantField bgField;
      bgField.initialize(0,0,this->B0); //bg bx, by,bz
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
} // namespace projects
