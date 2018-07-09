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

#include "verificationLarmor.h"

namespace projects {
   verificationLarmor::verificationLarmor(): Project() { }
   verificationLarmor::~verificationLarmor() { }
   bool verificationLarmor::initialize(void) {return true;}

   void verificationLarmor::addParameters() {
      typedef Readparameters RP;
      RP::add("Larmor.BX0", "Background field value (T)", 0.0);
      RP::add("Larmor.BY0", "Background field value (T)", 0.0);
      RP::add("Larmor.BZ0", "Background field value (T)", 0.0);
      RP::add("Larmor.VX0", "Bulk velocity in x", 0.0);
      RP::add("Larmor.VY0", "Bulk velocity in y", 0.0);
      RP::add("Larmor.VZ0", "Bulk velocity in z", 0.0);
      RP::add("Larmor.X0", "Initial Position", 0.0);
      RP::add("Larmor.Y0", "Initial Position", 0.0);
      RP::add("Larmor.Z0", "Initial Position", 0.0);
      RP::add("Larmor.rho", "Number density (m^-3)", 1.0e7);
   }

   void verificationLarmor::getParameters() {
      typedef Readparameters RP;
      RP::get("Larmor.BX0", this->BX0);
      RP::get("Larmor.BY0", this->BY0);
      RP::get("Larmor.BZ0", this->BZ0);
      RP::get("Larmor.VX0", this->VX0);
      RP::get("Larmor.VY0", this->VY0);
      RP::get("Larmor.VZ0", this->VZ0);
      RP::get("Larmor.X0", this->X0);
      RP::get("Larmor.Y0", this->Y0);
      RP::get("Larmor.Z0", this->Z0);
      RP::get("Larmor.rho", this->DENSITY);
   }

   Real verificationLarmor::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {

      static bool isSet=false;
      //static variables should be threadprivate
   #pragma omp threadprivate(isSet)

      if(vx < Parameters::vxmin + 0.5 * dvx ||
         vy < Parameters::vymin + 0.5 * dvy ||
         vz < Parameters::vzmin + 0.5 * dvz ||
         vx > Parameters::vxmax - 1.5 * dvx ||
         vy > Parameters::vymax - 1.5 * dvy ||
         vz > Parameters::vzmax - 1.5 * dvz
      ) return 0.0;

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


   void verificationLarmor::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
      cellParams[CellParams::BGBX   ] = this->BX0;
      cellParams[CellParams::BGBY   ] = this->BY0;
      cellParams[CellParams::BGBZ   ] = this->BZ0;

   }

} //namespace projects
