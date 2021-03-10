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

#include "VelocityBox.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   VelocityBox::VelocityBox(): Project() { }
   VelocityBox::~VelocityBox() { }


   bool VelocityBox::initialize(void) {return Project::initialize();}

   void VelocityBox::addParameters(){
      typedef Readparameters RP;
      RP::add("VelocityBox.rho", "Number density in full 6 dimensions (m^-6 s^3)", 0.0);
      RP::add("VelocityBox.Vx1", "Box min x (m/s)", 0.0);
      RP::add("VelocityBox.Vx2", "Box max x (m/s)", 0.0);
      RP::add("VelocityBox.Vy1", "Box min y (m/s)", 0.0);
      RP::add("VelocityBox.Vy2", "Box max y (m/s)", 0.0);
      RP::add("VelocityBox.Vz1", "Box min z (m/s)", 0.0);
      RP::add("VelocityBox.Vz2", "Box max z (m/s)", 0.0);
      RP::add("VelocityBox.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("VelocityBox.By", "Magnetic field y component (T)", 0.0);
      RP::add("VelocityBox.Bz", "Magnetic field z component (T)", 0.0);
   }

   void VelocityBox::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("VelocityBox.rho", this->rho);
      RP::get("VelocityBox.Vx1", this->Vx[0]);
      RP::get("VelocityBox.Vx2", this->Vx[1]);
      RP::get("VelocityBox.Vy1", this->Vy[0]);
      RP::get("VelocityBox.Vy2", this->Vy[1]);
      RP::get("VelocityBox.Vz1", this->Vz[0]);
      RP::get("VelocityBox.Vz2", this->Vz[1]);
      RP::get("VelocityBox.Bx", this->Bx);
      RP::get("VelocityBox.By", this->By);
      RP::get("VelocityBox.Bz", this->Bz);
   }

  Real VelocityBox::getDistribValue(creal& vx, creal& vy, creal& vz, const uint popID) const {
     if (vx >= this->Vx[0] && vx <= this->Vx[1] &&
         vy >= this->Vy[0] && vy <= this->Vy[1] &&
         vz >= this->Vz[0] && vz <= this->Vz[1])
       return this->rho;
     else
       return 0.0;
   }



  Real VelocityBox::calcPhaseSpaceDensity(
     creal& x, creal& y, creal& z,
     creal& dx, creal& dy, creal& dz,
     creal& vx, creal& vy, creal& vz,
     creal& dvx, creal& dvy, creal& dvz,const uint popID
  ) const {
    return getDistribValue(vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, popID);
  }


  
   void VelocityBox::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void VelocityBox::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField, BgBGrid);
   }
   
}// namespace projects
