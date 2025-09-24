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

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../object_wrapper.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../common.h"
#include "../../parameters.h"
#include "../../readparameters.h"
#include "../../vlasovsolver/vlasovmover.h"
#include "../project.h"
#include "Shocktest.h"

using namespace std;
using namespace spatial_cell;

namespace projects {

   Shocktest::Shocktest() : TriAxisSearch() {} // Constructor
   Shocktest::~Shocktest() {}                  // Destructor

   bool Shocktest::initialize(void) { return Project::initialize(); }

   void Shocktest::addParameters() {
      typedef Readparameters RP;
      RP::add("Shocktest.rho1", "Number density, left state (m^-3)", 0.0);
      RP::add("Shocktest.rho2", "Number density, right state (m^-3)", 0.0);
      RP::add("Shocktest.T1", "Temperature, left state (K)", 0.0);
      RP::add("Shocktest.T2", "Temperature, right state (K)", 0.0);
      RP::add("Shocktest.Vx1", "Bulk velocity x component, left state (m/s)", 0.0);
      RP::add("Shocktest.Vx2", "Bulk velocity x component, right state (m/s)", 0.0);
      RP::add("Shocktest.Vy1", "Bulk velocity y component, left state (m/s)", 0.0);
      RP::add("Shocktest.Vy2", "Bulk velocity y component, right state (m/s)", 0.0);
      RP::add("Shocktest.Vz1", "Bulk velocity z component, left state (m/s)", 0.0);
      RP::add("Shocktest.Vz2", "Bulk velocity z component, right state (m/s)", 0.0);
      RP::add("Shocktest.Bx1", "Magnetic field x component, left state (T)", 0.0);
      RP::add("Shocktest.Bx2", "Magnetic field x component, right state (T)", 0.0);
      RP::add("Shocktest.By1", "Magnetic field y component, left state (T)", 0.0);
      RP::add("Shocktest.By2", "Magnetic field y component, right state (T)", 0.0);
      RP::add("Shocktest.Bz1", "Magnetic field z component, left state (T)", 0.0);
      RP::add("Shocktest.Bz2", "Magnetic field z component, right state (T)", 0.0);
   }

   void Shocktest::getParameters() {
      Project::getParameters();

      if (getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      this->rho[this->LEFT] = {NAN};
      this->T[this->LEFT] = {NAN};
      this->Vx[this->LEFT] = {NAN};
      this->Vy[this->LEFT] = {NAN};
      this->Vz[this->LEFT] = {NAN};
      this->Bx[this->LEFT] = {NAN};
      this->By[this->LEFT] = {NAN};
      this->Bz[this->LEFT] = {NAN};
      this->rho[this->RIGHT] = {NAN};
      this->T[this->RIGHT] = {NAN};
      this->Vx[this->RIGHT] = {NAN};
      this->Vy[this->RIGHT] = {NAN};
      this->Vz[this->RIGHT] = {NAN};
      this->Bx[this->RIGHT] = {NAN};
      this->By[this->RIGHT] = {NAN};
      this->Bz[this->RIGHT] = {NAN};

      typedef Readparameters RP;
      RP::get("Shocktest.rho1", this->rho[this->LEFT]);
      RP::get("Shocktest.rho2", this->rho[this->RIGHT]);
      RP::get("Shocktest.T1", this->T[this->LEFT]);
      RP::get("Shocktest.T2", this->T[this->RIGHT]);
      RP::get("Shocktest.Vx1", this->Vx[this->LEFT]);
      RP::get("Shocktest.Vx2", this->Vx[this->RIGHT]);
      RP::get("Shocktest.Vy1", this->Vy[this->LEFT]);
      RP::get("Shocktest.Vy2", this->Vy[this->RIGHT]);
      RP::get("Shocktest.Vz1", this->Vz[this->LEFT]);
      RP::get("Shocktest.Vz2", this->Vz[this->RIGHT]);
      RP::get("Shocktest.Bx1", this->Bx[this->LEFT]);
      RP::get("Shocktest.Bx2", this->Bx[this->RIGHT]);
      RP::get("Shocktest.By1", this->By[this->LEFT]);
      RP::get("Shocktest.By2", this->By[this->RIGHT]);
      RP::get("Shocktest.Bz1", this->Bz[this->LEFT]);
      RP::get("Shocktest.Bz2", this->Bz[this->RIGHT]);
   }

   Realf Shocktest::fillPhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, const uint nRequested) const {
      // const speciesParameters& sP = this->speciesParams[popID];
      //  Fetch spatial cell center coordinates
      const Real x = cell->parameters[CellParams::XCRD] + 0.5 * cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;
      Real initRho = this->rho[side];
      Real initT = this->T[side];
      const Real initV0X = this->Vx[side];
      const Real initV0Y = this->Vy[side];
      const Real initV0Z = this->Vz[side];

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
   Realf Shocktest::probePhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, Real vx_in, Real vy_in, Real vz_in) const {
      // const speciesParameters& sP = this->speciesParams[popID];
      //  Fetch spatial cell center coordinates
      const Real x = cell->parameters[CellParams::XCRD] + 0.5 * cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;
      Real initRho = this->rho[side];
      Real initT = this->T[side];
      const Real initV0X = this->Vx[side];
      const Real initV0Y = this->Vy[side];
      const Real initV0Z = this->Vz[side];

      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx, vy, vz, initT, initRho, mass);
      return value;
   }

   /** Returns the center coordinates of the maxwellian distribution
   @ param x The x coordinate of the given spatial cell
   @ param y The x coordinate of the given spatial cell
   @ param z The x coordinate of the given spatial cell
   */
   vector<std::array<Real, 3>> Shocktest::getV0(creal x, creal y, creal z, const uint popID) const {
      vector<std::array<Real, 3>> centerPoints;
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;
      std::array<Real, 3> V0{{this->Vx[side], this->Vy[side], this->Vz[side]}};
      centerPoints.push_back(V0);
      return centerPoints;
   }

   /** Calculate parameters for the given spatial cell at the given time.
    * @param cellParams Array containing cell parameters.
    * @param t The current value of time. This is passed as a convenience. If you need more detailed information
    * of the state of the simulation, you can read it from Parameters.
    */
   void Shocktest::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

   void Shocktest::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      setBackgroundFieldToZero(BgBGrid);

      if (!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

         #pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);

                  cell->at(fsgrids::bfield::PERBX) = (xyz[0] < 0.0) ? this->Bx[this->LEFT] : this->Bx[this->RIGHT];
                  cell->at(fsgrids::bfield::PERBY) = (xyz[0] < 0.0) ? this->By[this->LEFT] : this->By[this->RIGHT];
                  cell->at(fsgrids::bfield::PERBZ) = (xyz[0] < 0.0) ? this->Bz[this->LEFT] : this->Bz[this->RIGHT];
               }
            }
         }
      }
   }

} // Namespace projects
