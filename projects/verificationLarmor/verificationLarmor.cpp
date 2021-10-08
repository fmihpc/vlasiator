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
#include "verificationLarmor.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   verificationLarmor::verificationLarmor(): Project() { }
   verificationLarmor::~verificationLarmor() { }
   bool verificationLarmor::initialize(void) {return Project::initialize();}

   void verificationLarmor::addParameters() {
      typedef Readparameters RP;
      RP::add("VerificationLarmor.BX0", "Background field value (T)", 0.0);
      RP::add("VerificationLarmor.BY0", "Background field value (T)", 0.0);
      RP::add("VerificationLarmor.BZ0", "Background field value (T)", 0.0);
      RP::add("VerificationLarmor.VX0", "Bulk velocity in x", 0.0);
      RP::add("VerificationLarmor.VY0", "Bulk velocity in y", 0.0);
      RP::add("VerificationLarmor.VZ0", "Bulk velocity in z", 0.0);
      RP::add("VerificationLarmor.X0", "Initial Position", 0.0);
      RP::add("VerificationLarmor.Y0", "Initial Position", 0.0);
      RP::add("VerificationLarmor.Z0", "Initial Position", 0.0);
      RP::add("VerificationLarmor.rho", "Number density (m^-3)", 1.0e7);
   }

   void verificationLarmor::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("VerificationLarmor.BX0", this->BX0);
      RP::get("VerificationLarmor.BY0", this->BY0);
      RP::get("VerificationLarmor.BZ0", this->BZ0);
      RP::get("VerificationLarmor.VX0", this->VX0);
      RP::get("VerificationLarmor.VY0", this->VY0);
      RP::get("VerificationLarmor.VZ0", this->VZ0);
      RP::get("VerificationLarmor.X0", this->X0);
      RP::get("VerificationLarmor.Y0", this->Y0);
      RP::get("VerificationLarmor.Z0", this->Z0);
      RP::get("VerificationLarmor.rho", this->DENSITY);
   }

   Real verificationLarmor::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz,
           creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {

      static bool isSet=false;
      //static variables should be threadprivate
   #pragma omp threadprivate(isSet)

      const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
      const vmesh::MeshParameters& meshParams = getObjectWrapper().velocityMeshes[meshID];
      if (vx < meshParams.meshMinLimits[0] + 0.5*dvx ||
          vy < meshParams.meshMinLimits[1] + 0.5*dvy ||
          vz < meshParams.meshMinLimits[2] + 0.5*dvz ||
          vx > meshParams.meshMaxLimits[0] - 1.5*dvx ||
          vy > meshParams.meshMaxLimits[1] - 1.5*dvy ||
          vz > meshParams.meshMaxLimits[2] - 1.5*dvz) {
         return 0.0;
      }

      if(isSet)
         return 0.0; //exactly one value to be set

      if( fabs(vx-this->VX0)<dvx &&
         fabs(vy-this->VY0)<dvy &&
         fabs(vz-this->VZ0)<dvz &&
         fabs(x-this->X0)<dx &&
         fabs(y-this->Y0)<dy &&
         fabs(z-this->Z0)<dz){
         isSet=true;
         return this->DENSITY/(dvx*dvy*dvz);
      }

      return 0.0;
      
   }


   void verificationLarmor::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void verificationLarmor::setProjectBField(
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
