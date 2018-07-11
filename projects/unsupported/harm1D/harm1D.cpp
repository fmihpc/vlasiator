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

#include "harm1D.h"

namespace projects {
   harm1D::harm1D(): Project() { }
   harm1D::~harm1D() { }
   
   using namespace std;
   
   bool harm1D::initialize(void) {return true;}
   void harm1D::addParameters() { }
   void harm1D::getParameters() { }

   Real harm1D::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz
   ) {
      /*
      creal VX0 = 0.5;
      creal VY0 = 0.0;
      creal VZ0 = 0.0;
      creal SIGMA = 0.4714;
      creal INVSIG2 = 1.0/(SIGMA*SIGMA);
      return 1.5*exp(-INVSIG2*(vx-VX0)*(vx-VX0))*exp(-INVSIG2*(vy-VY0)*(vy-VY0))*exp(-INVSIG2*(vz-VZ0)*(vz-VZ0));
      */
      creal X0 = 1.0/14.0;
      creal Y0 = 1.0/14.0;
      
      creal VX0 = -0.4;
      creal VY0 = -0.4;
      creal VZ0 = 0.0;
      creal DVX = 0.1;
      creal DVY = 0.1;
      creal DVZ = 0.1;
      creal VSIGMA = 0.2;
      creal INVVSIG2 = 1.0/(VSIGMA*VSIGMA);

      if (fabs(x + 0.6) > dx) return 1e-10;
      if (fabs(vx) > 0.051) return 1e-10;
      if (fabs(vy) > 0.8) return 1e-10;
      if (fabs(vz) > 0.8) return 1e-10;
      //if (fabs(x) > X0 || fabs(y) > Y0) return 0.0;
      //if (fabs(vy-VY0) > DVY) return 0.0;
      //if (fabs(vz-VZ0) > DVZ) return 0.0;
      //return 5.0*exp(-INVVSIG2*(vx-VX0)*(vx-VX0))*exp(-INVVSIG2*(vy-VY0)*(vy-VY0))*exp(-INVVSIG2*(vz-VZ0)*(vz-VZ0));
      return 1.0;
   }

   void harm1D::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];

      // Setting these is not needed for correct propagation, 
      // but may be a good idea for visualization:
      cellParams[CellParams::EX   ] = -1.0*(x+0.5*dx);
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
      cellParams[CellParams::BGBX   ] = 0.0;
      cellParams[CellParams::BGBY   ] = 0.0;
      cellParams[CellParams::BGBZ   ] = 0.0;
   }
} // namespace projects
