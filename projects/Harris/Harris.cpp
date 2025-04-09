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

#include "Harris.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Harris::Harris(): TriAxisSearch() { }
   Harris::~Harris() { }

   bool Harris::initialize(void) {return Project::initialize();}

   void Harris::addParameters(){
      typedef Readparameters RP;
      RP::add("Harris.Scale_size", "Harris sheet scale size (m)", 150000.0);
      RP::add("Harris.BX0", "Magnetic field at infinity (T)", 8.33061003094e-8);
      RP::add("Harris.BY0", "Magnetic field at infinity (T)", 8.33061003094e-8);
      RP::add("Harris.BZ0", "Magnetic field at infinity (T)", 8.33061003094e-8);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Harris.Temperature", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Harris.rho", "Number density at infinity (m^-3)", 1.0e7);
      }
   }

   void Harris::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;
      RP::get("Harris.Scale_size", this->SCA_LAMBDA);
      RP::get("Harris.BX0", this->BX0);
      RP::get("Harris.BY0", this->BY0);
      RP::get("Harris.BZ0", this->BZ0);


      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         HarrisSpeciesParameters sP;

         RP::get(pop + "_Harris.Temperature", sP.TEMPERATURE);
         RP::get(pop + "_Harris.rho", sP.DENSITY);

         speciesParams.push_back(sP);
      }
   }

   Realf Harris::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const HarrisSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY;
      Real initT = sP.TEMPERATURE;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

      initRho *= (1.0 + 5.0 / pow(cosh(x / (this->SCA_LAMBDA)), 2.0));

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

   /* Evaluates local SpatialCell properties for the project and population,
      then evaluates the phase-space density at the given coordinates.
      Used as a probe for projectTriAxisSearch.
   */
   Realf Harris::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      const HarrisSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY;
      Real initT = sP.TEMPERATURE;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

      initRho *= (1.0 + 5.0 / pow(cosh(x / (this->SCA_LAMBDA)), 2.0));
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
      return value;
   }

   void Harris::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   vector<std::array<Real, 3>> Harris::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> V0;
      std::array<Real, 3> v = {{0.0, 0.0, 0.0 }};
      V0.push_back(v);
      return V0;
   }

   void Harris::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);

      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

         #pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);

                  cell->at(fsgrids::bfield::PERBX) = this->BX0 * tanh((xyz[1] + 0.5 * perBGrid.DY) / this->SCA_LAMBDA);
                  cell->at(fsgrids::bfield::PERBY) = this->BY0 * tanh((xyz[2] + 0.5 * perBGrid.DZ) / this->SCA_LAMBDA);
                  cell->at(fsgrids::bfield::PERBZ) = this->BZ0 * tanh((xyz[0] + 0.5 * perBGrid.DX) / this->SCA_LAMBDA);
               }
            }
         }
      }
   }

} // namespace projects
