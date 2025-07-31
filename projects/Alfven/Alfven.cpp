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
#include "../../object_wrapper.h"

#include "Alfven.h"

using namespace std;

namespace projects {
   Alfven::Alfven(): Project() { }
   Alfven::~Alfven() { }

   bool Alfven::initialize(void) {
      bool success = Project::initialize();

      Real norm = sqrt(this->Bx_guiding*this->Bx_guiding + this->By_guiding*this->By_guiding + this->Bz_guiding*this->Bz_guiding);
      this->Bx_guiding /= norm;
      this->By_guiding /= norm;
      this->Bz_guiding /= norm;
      this->ALPHA = atan(this->By_guiding/this->Bx_guiding);

      return success;
   }

   void Alfven::addParameters() {
      typedef Readparameters RP;
      RP::add("Alfven.B0", "Guiding field value (T)", 1.0e-10);
      RP::add("Alfven.Bx_guiding", "Guiding field x component", 1);
      RP::add("Alfven.By_guiding", "Guiding field y component", 0);
      RP::add("Alfven.Bz_guiding", "Guiding field z component", 0);
      RP::add("Alfven.Wavelength", "Wavelength (m)", 100000.0);
      RP::add("Alfven.A_mag", "Amplitude of the magnetic perturbation", 0.1);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Alfven.rho", "Number density (m^-3)", 1.0e8);
         RP::add(pop + "_Alfven.Temperature", "Temperature (K)", 0.86456498092);
         RP::add(pop + "_Alfven.A_vel", "Amplitude of the velocity perturbation", 0.1);

      }
   }

   void Alfven::getParameters(){
      Project::getParameters();

      typedef Readparameters RP;
      RP::get("Alfven.B0", this->B0);
      RP::get("Alfven.Bx_guiding", this->Bx_guiding);
      RP::get("Alfven.By_guiding", this->By_guiding);
      RP::get("Alfven.Bz_guiding", this->Bz_guiding);
      RP::get("Alfven.Wavelength", this->WAVELENGTH);
      RP::get("Alfven.A_mag", this->A_MAG);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         AlfvenSpeciesParameters sP;

         RP::get(pop + "_Alfven.rho", sP.rho);
         RP::get(pop + "_Alfven.Temperature",sP.T);
         RP::get(pop + "_Alfven.A_vel", sP.A_VEL);

         speciesParams.push_back(sP);
      }
   }

   Realf Alfven::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const AlfvenSpeciesParameters& sP = this->speciesParams[popID];

      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal mu0 = physicalconstants::MU_0;
      creal ALFVEN_VEL = this->B0 / sqrt(mu0 * sP.rho * mass);

      creal ksi = (x * cos(this->ALPHA) + y * sin(this->ALPHA)) / this->WAVELENGTH;
      creal initV0X = sP.A_VEL * ALFVEN_VEL * sin(this->ALPHA) * sin(2.0 * M_PI * ksi);
      creal initV0Y = - sP.A_VEL * ALFVEN_VEL * cos(this->ALPHA) * sin(2.0 * M_PI * ksi);
      creal initV0Z = - sP.A_VEL * ALFVEN_VEL * cos(2.0 * M_PI * ksi);

      Real initRho = sP.rho;
      Real initT = sP.T;

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

   void Alfven::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      //Real* cellParams = cell->get_cell_parameters();
      //creal x = cellParams[CellParams::XCRD];
      //creal dx = cellParams[CellParams::DX];
      //creal y = cellParams[CellParams::YCRD];
      //creal dy = cellParams[CellParams::DY];
      //
      //Real ksi = ((x + 0.5 * dx)  * cos(this->ALPHA) + (y + 0.5 * dy) * sin(this->ALPHA)) / this->WAVELENGTH;
      //Real dBxavg = sin(2.0 * M_PI * ksi);
      //Real dByavg = sin(2.0 * M_PI * ksi);
      //Real dBzavg = cos(2.0 * M_PI * ksi);
   }

   void Alfven::setProjectBField(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                 std::span<std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                 std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid) {
      const auto& gridSpacing = fsgrid.getGridSpacing();
      setBackgroundFieldToZero(fsgrid, technical, bgb);

      if (!P::isRestart) {
         const auto* localSize = &fsgrid.getLocalSize()[0];

         // local copies for lambda capture
         const auto ALPHA_l = this->ALPHA;
         const auto WAVELENGTH_l = this->WAVELENGTH;
         const auto B0_l = this->B0;
         const auto A_MAG_l = this->A_MAG;

         fsgrid.parallel_for_coords([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                                    phiprof::initializeTimer("setProjectBField"), technical,
                                    [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer, const std::array<Real, 3> xyz) {
            auto& cell = perb[stencil.ooo()];

            const Real dx = gridSpacing[0];
            const Real dy = gridSpacing[1];
            const Real ksi = ((xyz[0] + 0.5 * dx) * cos(ALPHA_l) + (xyz[1] + 0.5 * dy) * sin(ALPHA_l)) / WAVELENGTH_l;
            const Real dBxavg = sin(2.0 * M_PI * ksi);
            const Real dByavg = sin(2.0 * M_PI * ksi);
            const Real dBzavg = cos(2.0 * M_PI * ksi);

            cell[fsgrids::bfield::PERBX] = B0_l * cos(ALPHA_l) - A_MAG_l * B0_l * sin(ALPHA_l) * dBxavg;
            cell[fsgrids::bfield::PERBY] = B0_l * sin(ALPHA_l) + A_MAG_l * B0_l * cos(ALPHA_l) * dByavg;
            cell[fsgrids::bfield::PERBZ] = B0_l * A_MAG_l * dBzavg;
         });
      }
   }

} // namespace projects
