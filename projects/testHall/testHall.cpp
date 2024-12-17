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
#include "../../object_wrapper.h"

#include "testHall.h"

using namespace std;

namespace projects {
   TestHall::TestHall(): Project() { }
   TestHall::~TestHall() { }
   
   bool TestHall::initialize(void) {
      bool success = Project::initialize();
      this->constBgB[0] = 0.0;
      this->constBgB[1] = 0.0;
      this->constBgB[2] = 0.0;
      this->dipoleScalingFactor = 1.0;
      this->dipoleTilt = 0.0;
      this->noDipoleInSW = 0;
      return success;
   }
   
   void TestHall::addParameters(){
      typedef Readparameters RP;
      RP::add("TestHall.BX0", "Magnetic field x (T)", 1.0e-9);
      RP::add("TestHall.BY0", "Magnetic field y (T)", 1.0e-9);
      RP::add("TestHall.BZ0", "Magnetic field z (T)", 1.0e-9);
      RP::add("TestHall.VX0", "velocity x (m/s)", -1.0e3);
      RP::add("TestHall.VY0", "velocity y (m/s)", 1.0e3);
      RP::add("TestHall.VZ0", "velocity z (m/s)", 1.0e3);
      RP::add("TestHall.Temperature", "Temperature (K)", 1.0e6);
      RP::add("TestHall.rho", "Number density (m^-3)", 1.0e6);
   }
   
   void TestHall::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("TestHall.BX0", this->BX0);
      RP::get("TestHall.BY0", this->BY0);
      RP::get("TestHall.BZ0", this->BZ0);
      RP::get("TestHall.VX0", this->VX0);
      RP::get("TestHall.VY0", this->VY0);
      RP::get("TestHall.VZ0", this->VZ0);
      RP::get("TestHall.Temperature", this->TEMPERATURE);
      RP::get("TestHall.rho", this->DENSITY);
   }
   
   Real TestHall::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz,const uint popID
   ) const {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      return this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) * (
         exp(- mass * (pow(vx + 0.5 * dvx - this->VX0, 2.0) + pow(vy + 0.5 * dvy - this->VY0, 2.0) + pow(vz + 0.5 * dvz - this->VZ0, 2.0)) / (2.0 * kb * this->TEMPERATURE)));
   }
   
   void TestHall::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void TestHall::setProjectBField(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                   std::span<std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                   fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      setBackgroundFieldToZero(bgb);

      if(!P::isRestart) {
         const auto& localSize = technicalGrid.getLocalSize();

#pragma omp parallel for collapse(3)
         for (auto x = 0; x < localSize[0]; ++x) {
            for (auto y = 0; y < localSize[1]; ++y) {
               for (auto z = 0; z < localSize[2]; ++z) {
                  const auto xyz = technicalGrid.getPhysicalCoords(x, y, z);
                  const auto stencil = technicalGrid.makeStencil(x, y, z);
                  auto& cell = perb[stencil.center()];

                  cell[fsgrids::bfield::PERBX] = this->BX0 * cos(2.0 * M_PI * 1.0 * xyz[0] / (P::xmax - P::xmin)) *
                                                 cos(2.0 * M_PI * 1.0 * xyz[1] / (P::ymax - P::ymin)) *
                                                 cos(2.0 * M_PI * 1.0 * xyz[2] / (P::zmax - P::zmin));
                  cell[fsgrids::bfield::PERBY] = this->BY0 * cos(2.0 * M_PI * 1.0 * xyz[0] / (P::xmax - P::xmin)) *
                                                 cos(2.0 * M_PI * 1.0 * xyz[1] / (P::ymax - P::ymin)) *
                                                 cos(2.0 * M_PI * 1.0 * xyz[2] / (P::zmax - P::zmin));
                  cell[fsgrids::bfield::PERBZ] = this->BZ0 * cos(2.0 * M_PI * 1.0 * xyz[0] / (P::xmax - P::xmin)) *
                                                 cos(2.0 * M_PI * 1.0 * xyz[1] / (P::ymax - P::ymin)) *
                                                 cos(2.0 * M_PI * 1.0 * xyz[2] / (P::zmax - P::zmin));
               }
            }
         }
      }
   }

} // namespace projects
