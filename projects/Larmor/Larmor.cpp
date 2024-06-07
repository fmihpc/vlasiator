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
#include "../../object_wrapper.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "Larmor.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
    Larmor::Larmor(): Project() { }
    Larmor::~Larmor() { }


   bool Larmor::initialize(void) {return Project::initialize();}

    void Larmor::addParameters() {
      typedef Readparameters RP;
      RP::add("Larmor.BX0", "Background field value (T)", 0.0);
      RP::add("Larmor.BY0", "Background field value (T)", 0.0);
      RP::add("Larmor.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Larmor.VX0", "Bulk velocity in x", 0.0);
      RP::add("Larmor.VY0", "Bulk velocity in y", 0.0);
      RP::add("Larmor.VZ0", "Bulk velocuty in z", 0.0);
      RP::add("Larmor.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Larmor.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Larmor.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("Larmor.Scale_x", "Scale length in x (m)", 2.0e6);
      RP::add("Larmor.Scale_y", "Scale length in y (m)", 2.0e6);
    }

    void Larmor::getParameters() {
       Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }

      RP::get("Larmor.BX0", this->BX0);
      RP::get("Larmor.BY0", this->BY0);
      RP::get("Larmor.BZ0", this->BZ0);
      RP::get("Larmor.VX0", this->VX0);
      RP::get("Larmor.VY0", this->VY0);
      RP::get("Larmor.VZ0", this->VZ0);
      RP::get("Larmor.rho", this->DENSITY);
      RP::get("Larmor.Temperature", this->TEMPERATURE);
      RP::get("Larmor.maxwCutoff", this->maxwCutoff);
      RP::get("Larmor.Scale_x", this->SCA_X);
      RP::get("Larmor.Scale_y", this->SCA_Y);
    }

    Real Larmor::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, const uint popID) const {
      creal kb = physicalconstants::K_B;
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      
      return exp(- mass * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0)+ (vz-this->VZ0)*(vz-this->VZ0)) / (2.0 * kb * this->TEMPERATURE))*
      exp(-pow(x-Parameters::xmax/2.5, 2.0)/pow(this->SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/2.0, 2.0)/pow(this->SCA_Y, 2.0));
    }

    Real Larmor::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, 
            creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
       const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
      vmesh::MeshParameters& meshParams = getObjectWrapper().velocityMeshes[meshID];
      if (vx < meshParams.meshMinLimits[0] + 0.5*dvx ||
          vy < meshParams.meshMinLimits[1] + 0.5*dvy ||
          vz < meshParams.meshMinLimits[2] + 0.5*dvz ||
          vx > meshParams.meshMaxLimits[0] - 1.5*dvx ||
          vy > meshParams.meshMaxLimits[1] - 1.5*dvy ||
          vz > meshParams.meshMaxLimits[2] - 1.5*dvz) {
         return 0.0;
      }

      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;

      Real avg = getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, popID);
      creal result = avg *this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5);
      
      if(result < this->maxwCutoff) {
         return 0.0;
      } else {
         return result;
      }
   }


   void Larmor::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

    void Larmor::setProjectBField(
       FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
       FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
       FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
    ) {
      ConstantField bgField;
      bgField.initialize(this->BX0,
                         this->BY0,
                         this->BZ0);
      
      setBackgroundField(bgField, BgBGrid);
   }
} //namespace projects 
  
