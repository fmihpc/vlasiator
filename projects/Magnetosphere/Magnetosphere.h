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

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "MultiPeak.h"

using namespace std;
using namespace spatial_cell;

Real projects::MultiPeak::rhoRnd;

namespace projects {
   MultiPeak::MultiPeak() : TriAxisSearch() {}

   MultiPeak::~MultiPeak() {}

   bool MultiPeak::initialize(void) { return Project::initialize(); }

   void MultiPeak::addParameters() {
      typedef Readparameters RP;

      RP::add("MultiPeak.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("MultiPeak.By", "Magnetic field y component (T)", 0.0);
      RP::add("MultiPeak.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("MultiPeak.dBx", "Magnetic field x component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBy", "Magnetic field y component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBz", "Magnetic field z component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.magXPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along x (T)", 1.0e-9);
      RP::add("MultiPeak.magYPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along y (T)", 1.0e-9);
      RP::add("MultiPeak.magZPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along z (T)", 1.0e-9);
      RP::add("MultiPeak.lambda", "B cosine perturbation wavelength (m)", 1.0);
      RP::add("MultiPeak.densityModel", "Which spatial density model is used?", string("uniform"));

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         RP::add(pop + "_MultiPeak.n", "Number of peaks to create", 0);
         RP::addComposing(pop + "_MultiPeak.rho", "Number density (m^-3)");
         RP::addComposing(pop + "_MultiPeak.Tx", "Temperature (K)");
         RP::addComposing(pop + "_MultiPeak.Ty", "Temperature");
         RP::addComposing(pop + "_MultiPeak.Tz", "Temperature");
         RP::addComposing(pop + "_MultiPeak.Vx", "Bulk velocity x component (m/s)");
         RP::addComposing(pop + "_MultiPeak.Vy", "Bulk velocity y component (m/s)");
         RP::addComposing(pop + "_MultiPeak.Vz", "Bulk velocity z component (m/s)");
         RP::addComposing(pop + "_MultiPeak.rhoPertAbsAmp", "Absolute amplitude of the density perturbation");
      }
   }

   void MultiPeak::getParameters() {

      typedef Readparameters RP;
      Project::getParameters();
      RP::get("MultiPeak.Bx", this->Bx);
      RP::get("MultiPeak.By", this->By);
      RP::get("MultiPeak.Bz", this->Bz);
      RP::get("MultiPeak.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("MultiPeak.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("MultiPeak.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("MultiPeak.dBx", this->dBx);
      RP::get("MultiPeak.dBy", this->dBy);
      RP::get("MultiPeak.dBz", this->dBz);
      RP::get("MultiPeak.lambda", this->lambda);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         MultiPeakSpeciesParameters sP;
         RP::get(pop + "_MultiPeak.n", sP.numberOfPeaks);
         RP::get(pop + "_MultiPeak.rho", sP.rho);
         RP::get(pop + "_MultiPeak.Tx", sP.Tx);
         RP::get(pop + "_MultiPeak.Ty", sP.Ty);
         RP::get(pop + "_MultiPeak.Tz", sP.Tz);
         RP::get(pop + "_MultiPeak.Vx", sP.Vx);
         RP::get(pop + "_MultiPeak.Vy", sP.Vy);
         RP::get(pop + "_MultiPeak.Vz", sP.Vz);

         RP::get(pop + "_MultiPeak.rhoPertAbsAmp", sP.rhoPertAbsAmp);
         if (!sP.isConsistent()) {
            cerr << "You should define all parameters (MultiPeak.rho, MultiPeak.Tx, MultiPeak.Ty, MultiPeak.Tz, MultiPeak.Vx, MultiPeak.Vy, MultiPeak.Vz, MultiPeak.rhoPertAbsAmp) for all "
                 << sP.numberOfPeaks << " peaks of population " << pop << "." << endl;
            abort();
         }

         speciesParams.push_back(sP);
      }

      string densModelString;
      RP::get("MultiPeak.densityModel", densModelString);

      if (densModelString == "uniform")
         densityModel = Uniform;
      else if (densModelString == "testcase")
         densityModel = TestCase;
   }

   Realf MultiPeak::fillPhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, const uint nRequested) const {
      const MultiPeakSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x = cell->parameters[CellParams::XCRD] + 0.5 * cell->parameters[CellParams::DX];
      const Real y = cell->parameters[CellParams::YCRD] + 0.5 * cell->parameters[CellParams::DY];
      const Real z = cell->parameters[CellParams::ZCRD] + 0.5 * cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;

      Real rhoFactor = 1.0;
      switch (densityModel) {
      case Uniform:
         rhoFactor = 1.0;
         break;
      case TestCase:
         rhoFactor = 1.0;
         if ((x >= 3.9e5 && x <= 6.1e5) && (y >= 3.9e5 && y <= 6.1e5)) {
            rhoFactor = 1.5;
         }
         break;
      default:
         rhoFactor = 1.0;
         break;
      }
      // device-accessable variables
      const Real rhoRndDev = rhoRnd;
      uint nPeaks = sP.numberOfPeaks;
      #define MAXPEAKS 10
      if (nPeaks > MAXPEAKS) {
         std::cerr << " ERROR in " << __FILE__ << ":" << __LINE__ << ": max number of supported peaks is " << MAXPEAKS << " (got " << nPeaks << ")" << std::endl;
         std::cerr << " Truncating peaks at " << MAXPEAKS << "!" << std::endl;
         nPeaks = MAXPEAKS;
      }
      Real VxDev[MAXPEAKS], VyDev[MAXPEAKS], VzDev[MAXPEAKS], TxDev[MAXPEAKS], TyDev[MAXPEAKS], TzDev[MAXPEAKS], rhoDev[MAXPEAKS], rhoPertAbsAmpDev[MAXPEAKS];
      for (uint i = 0; i < MAXPEAKS; ++i) {
         if (i >= nPeaks) {
            break;
         }
         VxDev[i] = sP.Vx[i];
         VyDev[i] = sP.Vy[i];
         VzDev[i] = sP.Vz[i];
         TxDev[i] = sP.Tx[i];
         TyDev[i] = sP.Ty[i];
         TzDev[i] = sP.Tz[i];
         rhoDev[i] = sP.rho[i];
         rhoPertAbsAmpDev[i] = sP.rhoPertAbsAmp[i];
      }

      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>({WID, WID, WID, nRequested},
                                        ARCH_LOOP_LAMBDA(const uint i, const uint j, const uint k, const uint initIndex, Realf* lsum) {
                                           vmesh::GlobalID* GIDlist = vmesh->getGrid()->data();
                                           Realf* bufferData = VBC->getData();
                                           const vmesh::GlobalID blockGID = GIDlist[initIndex];
                                           // Calculate parameters for new block
                                           Real blockCoords[6];
                                           vmesh->getBlockInfo(blockGID, &blockCoords[0]);
                                           creal vxBlock = blockCoords[0];
                                           creal vyBlock = blockCoords[1];
                                           creal vzBlock = blockCoords[2];
                                           creal dvxCell = blockCoords[3];
                                           creal dvyCell = blockCoords[4];
                                           creal dvzCell = blockCoords[5];
                                           ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
                                              Realf value = 0;
                                              for (uint ipeak = 0; ipeak < nPeaks; ++ipeak) {
                                                 creal vx = vxBlock + (i + 0.5) * dvxCell - VxDev[ipeak];
                                                 creal vy = vyBlock + (j + 0.5) * dvyCell - VyDev[ipeak];
                                                 creal vz = vzBlock + (k + 0.5) * dvzCell - VzDev[ipeak];
                                                 value += TriMaxwellianPhaseSpaceDensity(
                                                     vx, vy, vz, TxDev[ipeak], TyDev[ipeak], TzDev[ipeak], (rhoDev[ipeak] + rhoPertAbsAmpDev[ipeak] * rhoRndDev) * rhoFactor, mass);
                                              }
                                              bufferData[initIndex * WID3 + k * WID2 + j * WID + i] = value;
                                              // lsum[0] += value;
                                           };
                                        },
                                        rhosum);
      return rhosum;
   }

   /* Evaluates local SpatialCell properties for the project and population,
      then evaluates the phase-space density at the given coordinates.
      Used as a probe for projectTriAxisSearch.
   */
   Realf MultiPeak::probePhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, Real vx_in, Real vy_in, Real vz_in) const {
      const MultiPeakSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x = cell->parameters[CellParams::XCRD] + 0.5 * cell->parameters[CellParams::DX];
      const Real y = cell->parameters[CellParams::YCRD] + 0.5 * cell->parameters[CellParams::DY];
      const Real z = cell->parameters[CellParams::ZCRD] + 0.5 * cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;

      Real rhoFactor = 1.0;
      switch (densityModel) {
      case Uniform:
         rhoFactor = 1.0;
         break;
      case TestCase:
         rhoFactor = 1.0;
         if ((x >= 3.9e5 && x <= 6.1e5) && (y >= 3.9e5 && y <= 6.1e5)) {
            rhoFactor = 1.5;
         }
         break;
      default:
         rhoFactor = 1.0;
         break;
      }

      Realf value = 0;
      for (uint i = 0; i < sP.numberOfPeaks; ++i) {
         creal vx = vx_in - sP.Vx[i];
         creal vy = vy_in - sP.Vy[i];
         creal vz = vz_in - sP.Vz[i];
         value += TriMaxwellianPhaseSpaceDensity(vx, vy, vz, sP.Tx[i], sP.Ty[i], sP.Tz[i], sP.rho[i] + sP.rhoPertAbsAmp[i] * rhoRnd * rhoFactor, mass);
      }
      return value;
   }

   void MultiPeak::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {
      std::default_random_engine rndState;
      setRandomCellSeed(cell, rndState);
      rhoRnd = 0.5 - getRandomNumber(rndState);
   }

   void MultiPeak::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      ConstantField bgField;
      bgField.initialize(this->Bx, this->By, this->Bz);

      setBackgroundField(bgField, BgBGrid);

      if (!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
                  std::default_random_engine rndState;
                  setRandomSeed(cellid, rndState);

                  if (this->lambda != 0.0) {
                     cell->at(fsgrids::bfield::PERBX) = this->dBx * cos(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBY) = this->dBy * sin(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBZ) = this->dBz * cos(2.0 * M_PI * xyz[0] / this->lambda);
                  }

                  cell->at(fsgrids::bfield::PERBX) += this->magXPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBY) += this->magYPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBZ) += this->magZPertAbsAmp * (0.5 - getRandomNumber(rndState));
               }
            }
         }
      }
   }

   std::vector<std::array<Real, 3>> MultiPeak::getV0(creal x, creal y, creal z, const uint popID) const {
      const MultiPeakSpeciesParameters& sP = speciesParams[popID];
      vector<std::array<Real, 3>> centerPoints;
      for (uint i = 0; i < sP.numberOfPeaks; i++) {
         array<Real, 3> point{{sP.Vx[i], sP.Vy[i], sP.Vz[i]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }

} // namespace projects
