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

   Realf TestHall::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      //const speciesParameters& sP = this->speciesParams[popID];
      // Fetch spatial cell center coordinates
      // const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->DENSITY;
      Real initT = this->TEMPERATURE;
      const Real initV0X = this->VX0;
      const Real initV0Y = this->VY0;
      const Real initV0Z = this->VZ0;

      #ifdef USE_GPU
      vmesh::VelocityMesh *vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>(
         {WID, WID, WID, nRequested},
         ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) {
            vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
            Realf* bufferData = VBC->getData();
            const vmesh::GlobalID blockGID = GIDlist[initIndex];
            // Calculate parameters for new block
            Real blockCoords[6];
            vmesh->getBlockInfo(blockGID,&blockCoords[0]);
            creal vxBlock = blockCoords[0];
            creal vyBlock = blockCoords[1];
            creal vzBlock = blockCoords[2];
            creal dvxCell = blockCoords[3];
            creal dvyCell = blockCoords[4];
            creal dvzCell = blockCoords[5];
            ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
               creal vx = vxBlock + (i+0.5)*dvxCell - initV0X;
               creal vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
               creal vz = vzBlock + (k+0.5)*dvzCell - initV0Z;
               const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };
         }, rhosum);
      return rhosum;
   }

   void TestHall::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void TestHall::setProjectBField(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                   std::span<std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                   std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid) {
      setBackgroundFieldToZero(bgb);

      if(!P::isRestart) {
         // local copies for lambda capture
         const auto BX0_l = this->BX0;
         const auto BY0_l = this->BY0;
         const auto BZ0_l = this->BZ0;
         fsgrid.parallel_for_coords([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                                    phiprof::initializeTimer("setProjectBField-loop"), technical,
                                    [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer, const std::array<Real, 3> xyz) {
            auto& cell = perb[stencil.ooo()];

            cell[fsgrids::bfield::PERBX] = BX0_l * cos(2.0 * M_PI * 1.0 * xyz[0] / (P::xmax - P::xmin)) *
                                           cos(2.0 * M_PI * 1.0 * xyz[1] / (P::ymax - P::ymin)) *
                                           cos(2.0 * M_PI * 1.0 * xyz[2] / (P::zmax - P::zmin));
            cell[fsgrids::bfield::PERBY] = BY0_l * cos(2.0 * M_PI * 1.0 * xyz[0] / (P::xmax - P::xmin)) *
                                           cos(2.0 * M_PI * 1.0 * xyz[1] / (P::ymax - P::ymin)) *
                                           cos(2.0 * M_PI * 1.0 * xyz[2] / (P::zmax - P::zmin));
            cell[fsgrids::bfield::PERBZ] = BZ0_l * cos(2.0 * M_PI * 1.0 * xyz[0] / (P::xmax - P::xmin)) *
                                           cos(2.0 * M_PI * 1.0 * xyz[1] / (P::ymax - P::ymin)) *
                                           cos(2.0 * M_PI * 1.0 * xyz[2] / (P::zmax - P::zmin));
         });
      }
   }

} // namespace projects
