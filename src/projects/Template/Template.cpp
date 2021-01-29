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
#include "../../backgroundfield/dipole.hpp"

#include "Template.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Template::Template(): TriAxisSearch() { }
   Template::~Template() { }
   
   void Template::addParameters() {
      typedef Readparameters RP;
      RP::add("Template.param", "This is my project's parameter. Default is 0.0", 0.0);
   }
   
   void Template::getParameters(){
      Parameters::getParameters();

      typedef Readparameters RP;
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      RP::get("Template.param", this->param);
   }
   
   bool Template::initialize() {
      this->param += 1.0;
      return Project::initialize();
   }

   Real Template::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      creal rho = 1.0;
      creal T = 1.0;
      const std::array<Real, 3> V0 = this->getV0(x, y, z, popID)[0];
      creal Vx0 = V0[0];
      creal Vy0 = V0[1];
      creal Vz0 = V0[2];
      return rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
      exp(- physicalconstants::MASS_PROTON * ((vx-Vx0)*(vx-Vx0) + (vy-Vy0)*(vy-Vy0) + (vz-Vz0)*(vz-Vz0)) / (2.0 * physicalconstants::K_B * T));
   }
   
   void Template::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      Dipole bgField;
      bgField.initialize(8e15, 0.0, 0.0, 0.0, 0.0); //set dipole moment and location
      setBackgroundField(bgField, BgBGrid);
   }
   
   vector<std::array<Real, 3>> Template::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> centerPoints;
      std::array<Real, 3> point {{0.0, 0.0, 0.0}};
      if(x < 0.0) point[1] = 1.0;
      centerPoints.push_back(point);
      return centerPoints;
   }
   
} // namespace projects

