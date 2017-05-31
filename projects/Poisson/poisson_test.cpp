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
 * 
 * File:   poisson_test.cpp
 * Author: sandroos
 *
 * Created on January 14, 2015, 3:17 PM
 */

#include <cstdlib>
#include <iostream>

#include "../../readparameters.h"

#include "../../poisson_solver/poisson_solver.h"
#include "poisson_test.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   
   static Real radius = 0;
   
   PoissonTest::PoissonTest(): Project() { }
   
   PoissonTest::~PoissonTest() { }
   
   void PoissonTest::addParameters() {
      typedef Readparameters RP;
      RP::add("Poisson.solver","Name of the Poisson solver",string("SOR"));
      RP::add("Poisson.radius","Radius where charge density is non-zero",(Real)15e3);
      RP::add("Poisson.max_iterations","Maximum number of iterations",(uint)1000);
      RP::add("Poisson.min_relative_change","Potential is iterated until it the relative change is less than this value",(Real)1e-5);
      RP::add("Poisson.is_2D","If true then system is two-dimensional in xy-plane",true);
   }

   void PoissonTest::getParameters() {
      typedef Readparameters RP;
      RP::get("Poisson.solver",poisson::Poisson::solverName);
      RP::get("Poisson.radius",radius);
      RP::get("Poisson.max_iterations",poisson::Poisson::maxIterations);
      RP::get("Poisson.min_relative_change",poisson::Poisson::minRelativePotentialChange);
      RP::get("Poisson.is_2D",poisson::Poisson::is2D);
   }

   bool PoissonTest::initialize() {
      return true;
   }
   
   void PoissonTest::setCellBackgroundField(SpatialCell* cell) {

   }
   
   void PoissonTest::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      Real dx = cellParams[CellParams::DX];
      Real dy = cellParams[CellParams::DY];
      Real dz = cellParams[CellParams::DZ];
      Real x = cellParams[CellParams::XCRD] + 0.5*dx;
      Real y = cellParams[CellParams::YCRD] + 0.5*dy;
      Real z = cellParams[CellParams::ZCRD] + 0.5*dz;
      Real R = sqrt(x*x + y*y + z*z);
      
      cellParams[CellParams::RHOQ_TOT] = 0;
      
      if (R > radius) return;
      
      const Real volume = dx*dy*dz;
      cellParams[CellParams::RHOQ_TOT] = physicalconstants::CHARGE
              / volume / physicalconstants::EPS_0;
      
      if (Parameters::isRestart == true) return;
      cellParams[CellParams::PHI] = 0;
   }

   Real PoissonTest::calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,
            const uint popID) const {
      return 0.0;
   }

   std::vector<uint> PoissonTest::findBlocksToInitialize(SpatialCell* cell) {
      vector<uint> blockList;
      return blockList;
   }

} // namespace projects
