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
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../object_wrapper.h"
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

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Diffusion.rho", "Number density (m^-3)", 1.0e7);
         RP::add(pop + "_Diffusion.Temperature", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Diffusion.Scale_x", "Scale length in x (m)", 100000.0);
         RP::add(pop + "_Diffusion.Scale_y", "Scale length in y (m)", 100000.0);
      }
   }

   void Diffusion::getParameters() {
      Project::getParameters();

      typedef Readparameters RP;
      RP::get("Diffusion.B0", this->B0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        DiffusionSpeciesParameters sP;

        RP::get(pop + "_Diffusion.rho", sP.DENSITY);
        RP::get(pop + "_Diffusion.Temperature", sP.TEMPERATURE);
        RP::get(pop + "_Diffusion.Scale_x", sP.SCA_X);
        RP::get(pop + "_Diffusion.Scale_y", sP.SCA_Y);

        speciesParams.push_back(sP);
      }
   }

   Real Diffusion::getDistribValue(
      creal& x,creal& y,creal& z,
      creal& vx,creal& vy,creal& vz,
      const uint popID
   ) const {
      const DiffusionSpeciesParameters& sP = speciesParams[popID];
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      
      return sP.DENSITY * pow(mass / (2.0 * M_PI * kb * sP.TEMPERATURE), 1.5) * (
         5.0 * exp(- (pow(x, 2.0) / pow(sP.SCA_X, 2.0) +  pow(y, 2.0) / pow(sP.SCA_Y, 2.0))) * 
         exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * kb * sP.TEMPERATURE))
         +
         exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * kb * sP.TEMPERATURE)));
   }
   
   Real Diffusion::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      return getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, popID);
   }
   
   void Diffusion::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void Diffusion::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(0,0,this->B0); //bg bx, by,bz
      setBackgroundField(bgField, BgBGrid);
   }
} // namespace projects
