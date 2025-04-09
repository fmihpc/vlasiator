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

#include "Larmor.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
    Larmor::Larmor(): Project() { }
    Larmor::~Larmor() { }


   bool Larmor::initialize(void) {return Project::initialize();}

    void Larmor::addParameters() {
      typedef Readparameters RP;
      RP::add("Larmor.BX0", "Background field value (T)", 0.0);
      RP::add("Larmor.BY0", "Background field value (T)", 0.0);
      RP::add("Larmor.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Larmor.VX0", "Bulk velocity in x", 0.0);
      RP::add("Larmor.VY0", "Bulk velocity in y", 0.0);
      RP::add("Larmor.VZ0", "Bulk velocuty in z", 0.0);
      RP::add("Larmor.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Larmor.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Larmor.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("Larmor.Scale_x", "Scale length in x (m)", 2.0e6);
      RP::add("Larmor.Scale_y", "Scale length in y (m)", 2.0e6);
    }

    void Larmor::getParameters() {
       Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }

      RP::get("Larmor.BX0", this->BX0);
      RP::get("Larmor.BY0", this->BY0);
      RP::get("Larmor.BZ0", this->BZ0);
      RP::get("Larmor.VX0", this->VX0);
      RP::get("Larmor.VY0", this->VY0);
      RP::get("Larmor.VZ0", this->VZ0);
      RP::get("Larmor.rho", this->DENSITY);
      RP::get("Larmor.Temperature", this->TEMPERATURE);
      RP::get("Larmor.maxwCutoff", this->maxwCutoff);
      RP::get("Larmor.Scale_x", this->SCA_X);
      RP::get("Larmor.Scale_y", this->SCA_Y);
    }

   Realf Larmor::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      //const speciesParameters& sP = this->speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->DENSITY;
      Real initT = this->TEMPERATURE;
      const Real initV0X = this->VX0;
      const Real initV0Y = this->VY0;
      const Real initV0Z = this->VZ0;
      initRho = initRho * exp(-pow(x-Parameters::xmax/2.5, 2.0)/pow(this->SCA_X, 2.0)) * exp(-pow(y-Parameters::ymax/2.0, 2.0)/pow(this->SCA_Y, 2.0));

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

   void Larmor::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

    void Larmor::setProjectBField(
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
