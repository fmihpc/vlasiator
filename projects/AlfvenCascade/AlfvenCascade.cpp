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

#include <climits>
#include <cstddef>
#include <string>
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

#include "AlfvenCascade.h"

using namespace spatial_cell;

namespace projects {
AlfvenCascade::AlfvenCascade() : Project() {}
AlfvenCascade::~AlfvenCascade() {}

   bool AlfvenCascade::initialize(void) {
      bool success = Project::initialize();

      creal m = physicalconstants::MASS_PROTON;
      creal e = physicalconstants::CHARGE;
      creal kB = physicalconstants::K_B;
      creal gamma = 5.0 / 3.0;
      creal mu0 = physicalconstants::MU_0;

      rho0 = m * n0; // Mass density
      p0 = n0 * kB * T; // pressure

      // Calculate Alfvén speed
      VA = B / sqrt(mu0 * rho0);

      if (verbose) {
         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
         if (myRank == MASTER_RANK) {
            std::cout << "Initialized multi-wave turbulence simulation\n";
            std::cout << "Number of waves: " << this->nWaves << "\n";
            std::cout << "Background field strength: " << this->B << " T\n";
            std::cout << "Alfvén speed: " << VA << " m/s\n";
            std::cout << "Gaussian mask state: " << this->gaussianMask << "\n";
            std::cout << "Center of gaussian mask: " << this->gaussianMaskLocation << "\n";
            std::cout << "Width of gaussian mask: " << this->gaussianMaskWidth << "\n";
            
            for (int idx = 0; idx < nWaves; idx++) {
               std::cout << "\nWave " << idx + 1 << ":\n";
               std::cout << "Wavelength: " << this->wavesParams.at(idx).wavelength << " m\n";
               std::cout << "Amplitude: " << this->wavesParams.at(idx).amplitude << " m/s\n";
               std::cout << "Phase: " << this->wavesParams.at(idx).phase << " rad\n";
               std::cout << "Angle: " << this->angle * 180/M_PI << " degrees\n";
            }
         }
      }

      return success;
   }

   void AlfvenCascade::addParameters() {
      typedef Readparameters RP;
      
      RP::add<int>("AlfvenCascade.numberOfWaves", "Number of waves in the simulation", this->nWaves, 1000);
      RP::add<Real>("AlfvenCascade.n0", "Background density (1/m^3)", this->n0, 1e6);
      RP::add<Real>("AlfvenCascade.B", "Background magnetic field strength (T)", this->B, 1e-8);
      RP::add<Real>("AlfvenCascade.T", "Temperature (K)", this->T, 1e6);
      RP::add<Real>("AlfvenCascade.spectralIndex", "Power law index for initial spectrum", this->spectralIndex, -5.0/3.0);
      RP::add<int>("AlfvenCascade.randomSeed", "Seed for random phase generation", this->randomSeed, 12345);
      RP::add<bool>("AlfvenCascade.verbose", "Verbose output", this->verbose, true);
      RP::add<Real>("AlfvenCascade.angle", "Wave angle (rad)", this->angle,0.0);
      RP::add<bool>("AlfvenCascade.gaussianMask", "True if using gaussian mask to initial perturbation", this->gaussianMask, false);
      RP::add<Real>("AlfvenCascade.gaussianMaskLocation", "Location of gaussian mask", this->gaussianMaskLocation, 0.0);
      RP::add<Real>("AlfvenCascade.gaussianMaskWidth", "Width of gaussian mask", this->gaussianMaskWidth, 1.0);

      // per wave parameters
      // so, at this point, this->nWaves isnt actually set, so the default value of 1000 is used
      // so this initializes 1000 waves, however in the next function only nWaves waves are actually parsed for use, so no problemo
      for (int i=0; i<this->nWaves; i++) {
         WaveParameters* wv = new WaveParameters();
         this->wavesParamsRead.push_back(wv);
         RP::add<Real>("wave" + std::to_string(i+1) + ".wavelength", "Wavelength of wave (m)", wv->wavelength, 1.0);
         RP::add<Real>("wave" + std::to_string(i+1) + ".amplitude", "Velocity amplitude (m/s)", wv->amplitude, 1.0);
         RP::add<Real>("wave" + std::to_string(i+1) + ".phase", "Initial phase (rad)", wv->phase, 0.0);
      }
   }

   void AlfvenCascade::getParameters() {
      for (int i=0; i<this->nWaves; i++) {
         this->wavesParams.push_back(*this->wavesParamsRead.at(i));
      }
   }

   void AlfvenCascade::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

   Realf AlfvenCascade::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      // const AlfvenSpeciesParameters& sP = this->speciesParams[popID];

      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      creal mass = physicalconstants::MASS_PROTON;
      creal mu0 = physicalconstants::MU_0;
      Real ux = 0.0, uy = 0.0, uz = 0.0;

      for (int idx = 0; idx < this->nWaves; idx++) {
         WaveParameters wave = this->wavesParams.at(idx);

         Real cosalpha = cos(this->angle);
         Real sinalpha = sin(this->angle);
         Real kwave = 2 * M_PI / wave.wavelength;
         Real xpar = x * cosalpha + y * sinalpha;

         Real uperp = 0.0, upara = 0.0;

         if (!gaussianMask) {
            uperp = wave.amplitude * sin(kwave * xpar + wave.phase);
            upara = wave.amplitude * cos(kwave * xpar + wave.phase);
         } else {
            Real gaussianVal = exp(-((x - this->gaussianMaskLocation) * (x - this->gaussianMaskLocation)) / (2 * this->gaussianMaskWidth * this->gaussianMaskWidth));
            uperp = wave.amplitude * sin(kwave * xpar + wave.phase) * gaussianVal;
            upara = wave.amplitude * cos(kwave * xpar + wave.phase) * gaussianVal;
         }

         ux += -uperp * sinalpha;
         uy += uperp * cosalpha;
         uz += upara;
      }
      creal initV0X = ux;
      creal initV0Y = uy;
      creal initV0Z = uz;

      Real initRho = this->n0;
      Real initT = this->T;

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

   void AlfvenCascade::setProjectBField(fsgrids::perbspan perb,
                                 fsgrids::bgbspan bgb,
                                 fsgrids::technicalspan technical, FieldSolverGrid &fsgrid) {
      // Set background field
      setBackgroundFieldToZero(fsgrid, technical, bgb);

      if (!P::isRestart) {

         creal mu0 = physicalconstants::MU_0;
         // local copy for lambda capture
         const bool gaussianMask_l = this->gaussianMask;
         creal gaussianMaskLocation_l = this->gaussianMaskLocation;
         creal gaussianMaskWidth_l = this->gaussianMaskWidth;
         creal nWaves_l = this->nWaves;
         creal angle_l = this->angle;
         creal rho0_l = this->rho0;
         const auto waveParams_l = this->wavesParams;
         
         fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                              phiprof::initializeTimer("setProjectBField"), technical,
                              [=](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
            const std::array<Real, 3> xyz = coordinates.getPhysicalCoords(stencil.i, stencil.j, stencil.k);
            const std::array<Real, 3> gridSpacing = coordinates.physicalGridSpacing;
            auto& cell = perb[stencil.ooo()];

            Real Bx = 0.0, By = 0.0, Bz = 0.0;

            // Sum contributions from all waves
            for (int idx = 0; idx < nWaves_l; idx++) {
               const auto wave = waveParams_l.at(idx);

               Real cosalpha = cos(angle_l);
               Real sinalpha = sin(angle_l);
               Real kwave = 2 * M_PI / wave.wavelength;
               Real xpar = xyz[0] * cosalpha + xyz[1] * sinalpha;

               // Calculate B1 from v1 using Alfvén wave relation
               Real B1 = std::pow(-1.0,idx) * wave.amplitude * sqrt(mu0 * rho0_l);

               Real Bperp = 0.0, Bpara = 0.0;

               if (!gaussianMask_l) {
                  Bperp = B1 * sin(kwave * xpar + wave.phase);
                  Bpara = B1 * cos(kwave * xpar + wave.phase);
               } else {
                  Real gaussianVal = exp(-((xyz[0] - gaussianMaskLocation_l) * (xyz[0] - gaussianMaskLocation_l)) / (2 * gaussianMaskWidth_l * gaussianMaskWidth_l));
                  Bperp = B1 * sin(kwave * xpar + wave.phase) * gaussianVal;
                  Bpara = B1 * cos(kwave * xpar + wave.phase) * gaussianVal;
               }
                  Bx += -Bperp * sinalpha;
                  By += Bperp * cosalpha;
                  Bz += Bpara;
            }

            cell[fsgrids::bfield::PERBX] = Bx;
            cell[fsgrids::bfield::PERBY] = By;
            cell[fsgrids::bfield::PERBZ] = Bz;
         });
      }
   }
} // namespace projects