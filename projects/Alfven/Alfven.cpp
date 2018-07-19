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

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "Alfven.h"

using namespace std;

namespace projects {
   Alfven::Alfven(): Project() { }
   Alfven::~Alfven() { }
   
   bool Alfven::initialize(void) {
      bool success = Project::initialize();

      Real norm = sqrt(this->Bx_guiding*this->Bx_guiding + this->By_guiding*this->By_guiding + this->Bz_guiding*this->Bz_guiding);
      this->Bx_guiding /= norm;
      this->By_guiding /= norm;
      this->By_guiding /= norm;
      this->ALPHA = atan(this->By_guiding/this->Bx_guiding);
      
      return success;
   } 

   void Alfven::addParameters() {
      typedef Readparameters RP;
      RP::add("Alfven.B0", "Guiding field value (T)", 1.0e-10);
      RP::add("Alfven.Bx_guiding", "Guiding field x component", 1);
      RP::add("Alfven.By_guiding", "Guiding field y component", 0);
      RP::add("Alfven.Bz_guiding", "Guiding field z component", 0);
      RP::add("Alfven.Wavelength", "Wavelength (m)", 100000.0);
      RP::add("Alfven.A_mag", "Amplitude of the magnetic perturbation", 0.1);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Alfven.rho", "Number density (m^-3)", 1.0e8);
         RP::add(pop + "_Alfven.Temperature", "Temperature (K)", 0.86456498092);
         RP::add(pop + "_Alfven.A_vel", "Amplitude of the velocity perturbation", 0.1);
         RP::add(pop + "_Alfven.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_Alfven.nVelocitySamples", "Number of sampling points per velocity dimension", 5);

      }
   }

   void Alfven::getParameters(){
      Project::getParameters();
      
      typedef Readparameters RP;
      RP::get("Alfven.B0", this->B0);
      RP::get("Alfven.Bx_guiding", this->Bx_guiding);
      RP::get("Alfven.By_guiding", this->By_guiding);
      RP::get("Alfven.Bz_guiding", this->Bz_guiding);
      RP::get("Alfven.Wavelength", this->WAVELENGTH);
      RP::get("Alfven.A_mag", this->A_MAG);
      RP::get("Alfven.nSpaceSamples", this->nSpaceSamples);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         AlfvenSpeciesParameters sP;

         RP::get(pop + "_Alfven.rho", sP.rho);
         RP::get(pop + "_Alfven.Temperature",sP.T);
         RP::get(pop + "_Alfven.A_vel", sP.A_VEL);
         RP::get(pop + "_Alfven.nVelocitySamples", sP.nVelocitySamples);

         speciesParams.push_back(sP);
      }
   }

   Realf Alfven::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, const uint popID) const {
      const AlfvenSpeciesParameters& sP = speciesParams[popID];
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      creal mu0 = physicalconstants::MU_0;
      creal ALFVEN_VEL = this->B0 / sqrt(mu0 * sP.rho * mass);

      creal ksi = (x * cos(this->ALPHA) + y * sin(this->ALPHA)) / this->WAVELENGTH;
      creal Vx = sP.A_VEL * ALFVEN_VEL * sin(this->ALPHA) * sin(2.0 * M_PI * ksi);
      creal Vy = - sP.A_VEL * ALFVEN_VEL * cos(this->ALPHA) * sin(2.0 * M_PI * ksi);
      creal Vz = - sP.A_VEL * ALFVEN_VEL * cos(2.0 * M_PI * ksi);
   
      crealf den = sP.rho * pow(mass / (2.0 * M_PI * kb * sP.T), 1.5) *
      exp(- mass * (pow(vx - Vx, 2.0) + pow(vy - Vy, 2.0) + pow(vz - Vz, 2.0)) / (2.0 * kb * sP.T));
      return den;
   }
   
   Realf Alfven::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      return sampleVelocitySpace(x, y, z, dx, dy, dz, vx, vy, vz, dvx, dvy, dvz, popID, this->nSpaceSamples, speciesParams[popID].nVelocitySamples);
   }
   
   void Alfven::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      
      Real dBxavg, dByavg, dBzavg;
      dBxavg = dByavg = dBzavg = 0.0;
      Real d_x = dx / (this->nSpaceSamples - 1);
      Real d_y = dy / (this->nSpaceSamples - 1);

      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint j=0; j<this->nSpaceSamples; ++j)
      for (uint k=0; k<this->nSpaceSamples; ++k) {
         Real ksi = ((x + i * d_x)  * cos(this->ALPHA) + (y + j * d_y) * sin(this->ALPHA)) / this->WAVELENGTH;
         dBxavg += sin(2.0 * M_PI * ksi);
         dByavg += sin(2.0 * M_PI * ksi);
         dBzavg += cos(2.0 * M_PI * ksi);
      }
      cuint nPts = pow(this->nSpaceSamples, 3.0);
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      //Field below could laso be set as background field
      cellParams[CellParams::PERBX   ] = this->B0 * cos(this->ALPHA) - this->A_MAG * this->B0 * sin(this->ALPHA) * dBxavg / nPts;
      cellParams[CellParams::PERBY   ] = this->B0 * sin(this->ALPHA) + this->A_MAG * this->B0 * cos(this->ALPHA) * dByavg / nPts;
      cellParams[CellParams::PERBZ   ] = this->B0 * this->A_MAG * dBzavg / nPts;
   }
} // namespace projects
