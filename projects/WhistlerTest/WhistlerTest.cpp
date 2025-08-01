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

#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "WhistlerTest.h"

using namespace spatial_cell;

namespace projects {
WhistlerTest::WhistlerTest() : Project() {}
WhistlerTest::~WhistlerTest() {}

bool WhistlerTest::initialize(void) {
   bool success = Project::initialize();

   creal m = physicalconstants::MASS_PROTON;
   creal e = physicalconstants::CHARGE;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal mu0 = physicalconstants::MU_0;

   kx = 2 * M_PI / (Parameters::xmax-Parameters::xmin);
   ky = 2 * M_PI / (Parameters::ymax-Parameters::ymin);
   kz = 2 * M_PI / (Parameters::zmax-Parameters::zmin);

   Rinv[3][3] = {
      {1.0/3, 2.0/3, 2.0/3},
      {-0.89442719, 0.4472136, 0},
      {-0.298142397, -0.59628479, 0.74535599}
   };

   rho0 = m * n0; // Mass density

   // Calculate Alfvén speed
   VA = B0 / sqrt(mu0 * rho0);

   angle_rad = angle * M_PI / 180.0;

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cout << "Initialized multi-wave turbulence simulation\n";
         std::cout << "Background field strength: " << B0 << " T\n";
         std::cout << "Alfvén speed: " << VA << " m/s\n";
         std::cout << "Angle: " << angle << " degrees\n";
         std::cout << "Amplitude: " << amplitude << " T\n";
      }
   }

   return success;
}

void WhistlerTest::addParameters() {
   typedef Readparameters RP;

   RP::add("WhistlerTest.amplitude", "Velocity amplitude (m/s)",0.0);
   
   
   RP::add("WhistlerTest.n0", "Background density (1/m^3)", 1e6);
   RP::add("WhistlerTest.B0", "Background magnetic field strength (T)", 1e-8);
   RP::add("WhistlerTest.T0", "Temperature (K)", 1e6);
   RP::add("WhistlerTest.verbose", "Verbose output", 1);
   RP::add("WhistlerTest.angle", "Wave-B0 angle (degrees)",0.0);
}

void WhistlerTest::getParameters() {
   typedef Readparameters RP;
   Project::getParameters();

   RP::get("WhistlerTest.amplitude", amplitude);

   // Get scalar parameters
   RP::get("WhistlerTest.n0", n0);
   RP::get("WhistlerTest.B0", B0);
   RP::get("WhistlerTest.T0", T0);
   RP::get("WhistlerTest.verbose", verbose);
   RP::get("WhistlerTest.angle", angle);
}

void WhistlerTest::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

Realf WhistlerTest::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {

      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];

      creal mass = physicalconstants::MASS_PROTON;
      creal mu0 = physicalconstants::MU_0;
      Real ux = 0.0, uy = 0.0, uz = 0.0;

      // for (int idx = 0; idx < nWaves; idx++) {
      //    Real cosalpha = cos(angle);
      //    Real sinalpha = sin(angle);
      //    Real kwave = 2 * M_PI / wavelength.at(idx);
      //    Real xpar = x * cosalpha + y * sinalpha;
         
      //    Real uperp = amplitude.at(idx) * sin(kwave * xpar + phase.at(idx) * M_PI / 180);
      //    Real upara = amplitude.at(idx) * cos(kwave * xpar + phase.at(idx) * M_PI / 180);
         
      //    ux += -uperp * sinalpha;
      //    uy += uperp * cosalpha;
      //    uz += upara;
      // }
      creal initV0X = ux;
      creal initV0Y = uy;
      creal initV0Z = uz;

      Real initRho = n0;
      Real initT = T0;

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

void WhistlerTest::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                    FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Set background field
   ConstantField bgField;

   Real Bx0 = 0.0; By0 = 0.0, Bz0 = 0.0;
   Real Bksi0 = B0 * cos(angle_rad);
   Real Beta0 = B0 * sin(angle_rad);
   Bx0 = Bksi0 * Rinv[0][0] + Beta0 * Rinv[0][1];
   By0 = Bksi0 * Rinv[1][0] + Beta0 * Rinv[1][1];
   Bz0 = Bksi0 * Rinv[2][0] + Beta0 * Rinv[2][1];

   bgField.initialize(Bx0, By0, Bz0); // Set background field according to angle
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

   creal mu0 = physicalconstants::MU_0;

#pragma omp parallel for collapse(3)
      for (int i = 0; i < localSize[0]; ++i) {
         for (int j = 0; j < localSize[1]; ++j) {
            for (int k = 0; k < localSize[2]; ++k) {
               const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

               Real Bx = 0.0, By = 0.0, Bz = 0.0;

               // Sum contributions from all waves
               Real Bksi = 0.0;
               Real Beta = amplitude * sin(kx * x[0] + ky * x[1] + kz * x[2]);
               Real Bzeta = amplitude * cos(kx * x[0] + ky * x[1] + kz * x[2]);
               
               Bx += Beta * Rinv[0][1] + Bzeta * Rinv[0][2];
               By += Beta * Rinv[1][1] + Bzeta * Rinv[1][2];
               Bz += Beta * Rinv[2][1] + Bzeta * Rinv[2][2];

               cell->at(fsgrids::bfield::PERBX) = Bx;
               cell->at(fsgrids::bfield::PERBY) = By;
               cell->at(fsgrids::bfield::PERBZ) = Bz;
            }
         }
      }
   }
}

} // namespace projects