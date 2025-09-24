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
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "Flowthrough.h"

using namespace std;

/** Enumerates spatial density models Flowthrough project supports.
 * In most cases you want to use 'Maxwellian'. However, test package
 * uses 'SheetMaxwellian'.*/

enum DensityModel { Maxwellian, SheetMaxwellian, Square, Triangle, Sinewave };

static DensityModel densityModel;

namespace projects {
   Flowthrough::Flowthrough() : TriAxisSearch() {}
   Flowthrough::~Flowthrough() {}

   bool Flowthrough::initialize(void) { return Project::initialize(); }

   void Flowthrough::addParameters() {
      typedef Readparameters RP;
      RP::add("Flowthrough.emptyBox", "Is the simulation domain empty initially?", false);
      RP::add("Flowthrough.densityModel", "Plasma density model, 'Maxwellian' or 'SheetMaxwellian'", string("Maxwellian"));
      RP::add("Flowthrough.densityWidth", "Width of signal around origin", 6.e7);
      RP::add("Flowthrough.rescaleDensity", "Rescale VDF to match spatial ", false);
      RP::add("Flowthrough.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("Flowthrough.By", "Magnetic field y component (T)", 0.0);
      RP::add("Flowthrough.Bz", "Magnetic field z component (T)", 0.0);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         RP::add(pop + "_Flowthrough.rho", "Number density (m^-3)", 0.0);
         RP::add(pop + "_Flowthrough.rhoBase", "Background number density (m^-3)", 0.0);
         RP::add(pop + "_Flowthrough.T", "Temperature (K)", 0.0);
         RP::add(pop + "_Flowthrough.VX0", "Initial bulk velocity in x-direction", 0.0);
         RP::add(pop + "_Flowthrough.VY0", "Initial bulk velocity in y-direction", 0.0);
         RP::add(pop + "_Flowthrough.VZ0", "Initial bulk velocity in z-direction", 0.0);
      }
   }

   void Flowthrough::getParameters() {
      Project::getParameters();
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      typedef Readparameters RP;

      RP::get("Flowthrough.emptyBox", emptyBox);
      RP::get("Flowthrough.Bx", this->Bx);
      RP::get("Flowthrough.By", this->By);
      RP::get("Flowthrough.Bz", this->Bz);
      string densityModelString;
      RP::get("Flowthrough.densityModel", densityModelString);
      if (densityModelString == "Maxwellian")
         densityModel = Maxwellian;
      else if (densityModelString == "SheetMaxwellian")
         densityModel = SheetMaxwellian;
      else if (densityModelString == "Square")
         densityModel = Square;
      else if (densityModelString == "Triangle")
         densityModel = Triangle;
      else if (densityModelString == "Sinewave")
         densityModel = Sinewave;
      else {
         if (myRank == MASTER_RANK)
            cerr << __FILE__ << ":" << __LINE__ << " ERROR: Unknown option value!" << endl;
         exit(1);
      }
      RP::get("Flowthrough.densityWidth", this->densityWidth);
      RP::get("Flowthrough.rescaleDensity", this->rescaleDensityFlag);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         FlowthroughSpeciesParameters sP;

         RP::get(pop + "_Flowthrough.rho", sP.rho);
         RP::get(pop + "_Flowthrough.rhoBase", sP.rhoBase);
         RP::get(pop + "_Flowthrough.T", sP.T);
         RP::get(pop + "_Flowthrough.VX0", sP.V0[0]);
         RP::get(pop + "_Flowthrough.VY0", sP.V0[1]);
         RP::get(pop + "_Flowthrough.VZ0", sP.V0[2]);

         speciesParams.push_back(sP);
      }
   }
   Real Flowthrough::getCorrectNumberDensity(spatial_cell::SpatialCell* cell, const uint popID) const {
      const FlowthroughSpeciesParameters& sP = speciesParams[popID];
      Real rvalue;
      const Real x = cell->parameters[CellParams::XCRD] + 0.5 * cell->parameters[CellParams::DX];
      const Real y = cell->parameters[CellParams::YCRD] + 0.5 * cell->parameters[CellParams::DY];
      const Real z = cell->parameters[CellParams::ZCRD] + 0.5 * cell->parameters[CellParams::DZ];
      switch (densityModel) {
      case Maxwellian:
         rvalue = sP.rho;
         break;
      case SheetMaxwellian:
         rvalue = sqrt(x * x + y * y + z * z);
         if (rvalue <= 0.5 * densityWidth) {
            rvalue = 4 * sP.rho;
         } else {
            rvalue = 0;
         }
         break;
      case Square:
         if (abs(x) < 0.5 * densityWidth) {
            rvalue = 4 * sP.rho;
         } else {
            rvalue = 4 * sP.rhoBase;
            // rvalue = 0;
         }
         break;
      case Triangle:
         if (abs(x) < 0.5 * densityWidth) {
            rvalue = 4;
            rvalue *= (sP.rhoBase + (sP.rho - sP.rhoBase) * (1. - abs(x) / (0.5 * densityWidth)));
         } else {
            rvalue = 4 * sP.rhoBase;
            // rvalue = 0;
         }
         break;
      case Sinewave:
         if (abs(x) < 0.5 * densityWidth) {
            rvalue = 4;
            rvalue *= (sP.rhoBase + (sP.rho - sP.rhoBase) * (0.5 + 0.5 * cos(M_PI * x / (0.5 * densityWidth))));
         } else {
            rvalue = 4 * sP.rhoBase;
            // rvalue = 0;
         }
         break;
      default:
         rvalue = sP.rho;
         break;
      }
      return rvalue;
   }

   Realf Flowthrough::fillPhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, const uint nRequested) const {
      const FlowthroughSpeciesParameters& sP = speciesParams[popID];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->getCorrectNumberDensity(cell, popID);
      Real initT = sP.T;
      const Real initV0X = sP.V0[0];
      const Real initV0Y = sP.V0[1];
      const Real initV0Z = sP.V0[2];

      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif

      if (emptyBox == true) {
         Realf* bufferData = cell->get_velocity_blocks(popID)->getData();
         std::memset(bufferData, 0, nRequested * WID3 * sizeof(Realf));
         return 0;
      }

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
                                              creal vx = vxBlock + (i + 0.5) * dvxCell - initV0X;
                                              creal vy = vyBlock + (j + 0.5) * dvyCell - initV0Y;
                                              creal vz = vzBlock + (k + 0.5) * dvzCell - initV0Z;
                                              const Realf value = MaxwellianPhaseSpaceDensity(vx, vy, vz, initT, initRho, mass);
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
   Realf Flowthrough::probePhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, Real vx_in, Real vy_in, Real vz_in) const {
      const FlowthroughSpeciesParameters& sP = speciesParams[popID];
      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->getCorrectNumberDensity(cell, popID);
      Real initT = sP.T;
      const Real initV0X = sP.V0[0];
      const Real initV0Y = sP.V0[1];
      const Real initV0Z = sP.V0[2];

      if (emptyBox == true) {
         return 0;
      }
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx, vy, vz, initT, initRho, mass);
      return value;
   }

   void Flowthrough::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

   void Flowthrough::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                      FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      ConstantField bgField;
      bgField.initialize(Bx, By, Bz); // bg bx, by,bz
      setBackgroundField(bgField, BgBGrid);
   }

   std::vector<std::array<Real, 3>> Flowthrough::getV0(creal x, creal y, creal z, const uint popID) const {
      const FlowthroughSpeciesParameters& sP = speciesParams[popID];
      vector<std::array<Real, 3>> centerPoints;
      std::array<Real, 3> point{{sP.V0[0], sP.V0[1], sP.V0[2]}};
      centerPoints.push_back(point);
      return centerPoints;
   }

} // namespace projects
