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
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../object_wrapper.h"
#include "../../velocity_mesh_parameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "verificationLarmor.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   verificationLarmor::verificationLarmor(): Project() { }
   verificationLarmor::~verificationLarmor() { }
   bool verificationLarmor::initialize(void) {return Project::initialize();}

   void verificationLarmor::addParameters() {
      typedef Readparameters RP;
      RP::add("VerificationLarmor.BX0", "Background field value (T)", 0.0);
      RP::add("VerificationLarmor.BY0", "Background field value (T)", 0.0);
      RP::add("VerificationLarmor.BZ0", "Background field value (T)", 0.0);
      RP::add("VerificationLarmor.VX0", "Bulk velocity in x", 0.0);
      RP::add("VerificationLarmor.VY0", "Bulk velocity in y", 0.0);
      RP::add("VerificationLarmor.VZ0", "Bulk velocity in z", 0.0);
      RP::add("VerificationLarmor.X0", "Initial Position", 0.0);
      RP::add("VerificationLarmor.Y0", "Initial Position", 0.0);
      RP::add("VerificationLarmor.Z0", "Initial Position", 0.0);
      RP::add("VerificationLarmor.rho", "Number density (m^-3)", 1.0e7);
   }

   void verificationLarmor::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("VerificationLarmor.BX0", this->BX0);
      RP::get("VerificationLarmor.BY0", this->BY0);
      RP::get("VerificationLarmor.BZ0", this->BZ0);
      RP::get("VerificationLarmor.VX0", this->VX0);
      RP::get("VerificationLarmor.VY0", this->VY0);
      RP::get("VerificationLarmor.VZ0", this->VZ0);
      RP::get("VerificationLarmor.X0", this->X0);
      RP::get("VerificationLarmor.Y0", this->Y0);
      RP::get("VerificationLarmor.Z0", this->Z0);
      RP::get("VerificationLarmor.rho", this->DENSITY);
   }

   Realf verificationLarmor::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      // Fetch spatial cell low corner coordinates
      const Real x  = cell->parameters[CellParams::XCRD];
      const Real dx = cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD];
      const Real dy = cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD];
      const Real dz = cell->parameters[CellParams::DZ];
      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->DENSITY;

      // NOTE: This fill function does not have a GPU-supported version.
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
      Realf* bufferData = VBC->getData();

      // Values are only set in cell at X0,Y0,Z0. Otherwise return empty.
      if (fabs(x-this->X0)>=dx ||
          fabs(y-this->Y0)>=dy ||
          fabs(z-this->Z0)>=dz) {
         std::memset(bufferData, 0, nRequested*WID3*sizeof(Realf));
         return 0;
      }

      static bool isSet=false;
      //static variables should be threadprivate
      #pragma omp threadprivate(isSet)

      // Loop over blocks
      Realf rhosum = 0;
      for (uint blockLID=0; blockLID<nRequested; ++blockLID) {
         vmesh::GlobalID blockGID = GIDlist[blockLID];
         // Calculate parameters for block
         Real blockCoords[6];
         vmesh->getBlockInfo(blockGID,&blockCoords[0]);
         creal vxBlock = blockCoords[0];
         creal vyBlock = blockCoords[1];
         creal vzBlock = blockCoords[2];
         creal dvxCell = blockCoords[3];
         creal dvyCell = blockCoords[4];
         creal dvzCell = blockCoords[5];
         for (uint kc=0; kc<WID; ++kc) {
            for (uint jc=0; jc<WID; ++jc) {
               for (uint ic=0; ic<WID; ++ic) {
                  creal vx = vxBlock + (ic+0.5)*dvxCell - this->VX0;
                  creal vy = vyBlock + (jc+0.5)*dvyCell - this->VY0;
                  creal vz = vzBlock + (kc+0.5)*dvzCell - this->VZ0;
                  Realf value=0;
                  if (isSet) {
                     continue;
                  }
                  if (fabs(vx)<dvxCell &&
                      fabs(vy)<dvyCell &&
                      fabs(vz)<dvzCell) {
                      isSet = true;
                      value = initRho/(dvxCell*dvyCell*dvzCell);
                  }
                  bufferData[blockLID*WID3 + kc*WID2 + jc*WID + ic] = value;
                  rhosum += value;
               }
            }
         }
      } // End loop over blocks
      return rhosum;
   }

   void verificationLarmor::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void verificationLarmor::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->BX0,
                         this->BY0,
                         this->BZ0);

      setBackgroundField(bgField, BgBGrid);
   }

} //namespace projects
