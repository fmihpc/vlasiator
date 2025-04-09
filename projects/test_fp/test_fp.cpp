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

#include "test_fp.h"

enum cases {BXCASE,BYCASE,BZCASE,BALLCASE};

using namespace std;

namespace projects {
   test_fp::test_fp(): TriAxisSearch() { }
   test_fp::~test_fp() { }


   /*typedef test_fpParameters tfP;
   Real this->B0 = NAN;
   Real this->DENSITY = NAN;
   Real this->TEMPERATURE = NAN;
   Real this->ALPHA = NAN;
   int  this->CASE = 5;
   bool this->shear = false;
   */

   bool test_fp::initialize(void) {
      Project::initialize();
      this->ALPHA *= M_PI / 4.0;
      return true;
   }

   void test_fp::addParameters(void){
      typedef Readparameters RP;
      RP::add("test_fp.V0", "Velocity magnitude (m/s)", 1.0e6);
      RP::add("test_fp.B0", "Magnetic field value in the non-zero patch (T)", 1.0e-9);
      RP::add("test_fp.rho", "Number density (m^-3)", 1.0e7);
      RP::add("test_fp.Temperature", "Temperature (K)", 1.0e-6);
      RP::add("test_fp.angle", "Orientation of the propagation expressed in pi/4", 0.0);
      RP::add("test_fp.Bdirection", "Direction of the magnetic field (0:x, 1:y, 2:z, 3:all)", 0);
      RP::add("test_fp.shear", "Add a shear (if false, V=0.5 everywhere).", true);
   }

   void test_fp::getParameters(void){
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("test_fp.B0", this->B0);
      RP::get("test_fp.V0", this->V0);
      RP::get("test_fp.rho", this->DENSITY);
      RP::get("test_fp.Temperature", this->TEMPERATURE);
      RP::get("test_fp.angle", this->ALPHA);
      RP::get("test_fp.Bdirection", this->CASE);
      RP::get("test_fp.shear", this->shear);
   }

   Real test_fp::sign(creal value) const {
      if (abs(value) < 1e-5) return 0.0;
      else return value / abs(value);
   }

   Realf test_fp::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      //const speciesParameters& sP = this->speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->DENSITY;
      Real initT = this->TEMPERATURE;

      std::array<Real, 3> initV0 = this->getV0(x, y, z, popID)[0];
      const Real initV0X = initV0[0];
      const Real initV0Y = initV0[1];
      const Real initV0Z = initV0[2];

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
   Realf test_fp::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->DENSITY;
      Real initT = this->TEMPERATURE;

      std::array<Real, 3> initV0 = this->getV0(x, y, z, popID)[0];
      const Real initV0X = initV0[0];
      const Real initV0Y = initV0[1];
      const Real initV0Z = initV0[2];

      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
      return value;
   }

   void test_fp::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);

      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

         creal dx = perBGrid.DX * 3.5;
         creal dy = perBGrid.DY * 3.5;
         creal dz = perBGrid.DZ * 3.5;

         Real areaFactor = 1.0;

         #pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t i = 0; i < localSize[0]; ++i) {
            for (FsGridTools::FsIndex_t j = 0; j < localSize[1]; ++j) {
               for (FsGridTools::FsIndex_t k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

                  creal x = xyz[0] + 0.5 * perBGrid.DX;
                  creal y = xyz[1] + 0.5 * perBGrid.DY;
                  creal z = xyz[2] + 0.5 * perBGrid.DZ;

                  switch (this->CASE) {
                     case BXCASE:
                        cell->at(fsgrids::bfield::PERBX) = 0.1 * this->B0 * areaFactor;
                        //areaFactor = (CellParams::DY * CellParams::DZ) / (dy * dz);
                        if (y >= -dy && y <= dy)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBX) = this->B0 * areaFactor;
                        break;
                     case BYCASE:
                        cell->at(fsgrids::bfield::PERBY) = 0.1 * this->B0 * areaFactor;
                        //areaFactor = (CellParams::DX * CellParams::DZ) / (dx * dz);
                        if (x >= -dx && x <= dx)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBY) = this->B0 * areaFactor;
                        break;
                     case BZCASE:
                        cell->at(fsgrids::bfield::PERBZ) = 0.1 * this->B0 * areaFactor;
                        //areaFactor = (CellParams::DX * CellParams::DY) / (dx * dy);
                        if (x >= -dx && x <= dx)
                           if (y >= -dy && y <= dy)
                              cell->at(fsgrids::bfield::PERBZ) = this->B0 * areaFactor;
                        break;
                     case BALLCASE:
                        cell->at(fsgrids::bfield::PERBX) = 0.1 * this->B0 * areaFactor;
                        cell->at(fsgrids::bfield::PERBY) = 0.1 * this->B0 * areaFactor;
                        cell->at(fsgrids::bfield::PERBZ) = 0.1 * this->B0 * areaFactor;

                        //areaFactor = (CellParams::DX * CellParams::DY) / (dx * dy);

                        if (y >= -dy && y <= dy)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBX) = this->B0 * areaFactor;
                        if (x >= -dx && x <= dx)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBY) = this->B0 * areaFactor;
                        if (x >= -dx && x <= dx)
                           if (y >= -dy && y <= dy)
                              cell->at(fsgrids::bfield::PERBZ) = this->B0 * areaFactor;
                        break;
                  }
               }
            }
         }
      }
   }

   void test_fp::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {

   }

   vector<std::array<Real, 3>> test_fp::getV0(
      creal x,
      creal y,
      creal z,
      creal dx,
      creal dy,
      creal dz,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> centerPoints;

      Real VX=0.0,VY=0.0,VZ=0.0;
      if (this->shear == true)
      {
         //Real ksi;
         Real eta;
         switch (this->CASE) {
            case BXCASE:
               //ksi = ((y + 0.5 * dy)  * cos(this->ALPHA) + (z + 0.5 * dz) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(y + 0.5 * dy)  * sin(this->ALPHA) + (z + 0.5 * dz) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = 0.0;
               VY = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               VZ = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               break;
            case BYCASE:
               //ksi = ((z + 0.5 * dz)  * cos(this->ALPHA) + (x + 0.5 * dx) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(z + 0.5 * dz)  * sin(this->ALPHA) + (x + 0.5 * dx) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               VY = 0.0;
               VZ = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               break;
            case BZCASE:
               //ksi = ((x + 0.5 * dx)  * cos(this->ALPHA) + (y + 0.5 * dy) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(x + 0.5 * dx)  * sin(this->ALPHA) + (y + 0.5 * dy) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               VY = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               VZ = 0.0;
            break;
          case BALLCASE:
            std::cerr << "not implemented in " << __FILE__ << ":" << __LINE__ << std::endl;
            exit(1);
            break;
         }
      } else {
         switch (this->CASE) {
          case BXCASE:
            VX = 0.0;
            VY = cos(this->ALPHA) * 0.5;
            VZ = sin(this->ALPHA) * 0.5;
            break;
          case BYCASE:
            VX = sin(this->ALPHA) * 0.5;
            VY = 0.0;
            VZ = cos(this->ALPHA) * 0.5;
            break;
          case BZCASE:
            VX = cos(this->ALPHA) * 0.5;
            VY = sin(this->ALPHA) * 0.5;
            VZ = 0.0;
            break;
          case BALLCASE:
            VX = 0.5 / sqrt(3.0);
            VY = 0.5 / sqrt(3.0);
            VZ = 0.5 / sqrt(3.0);
            break;
         }
      }

      VX *= this->V0 * 2.0;
      VY *= this->V0 * 2.0;
      VZ *= this->V0 * 2.0;

      std::array<Real, 3> point {{VX, VY, VZ}};
      centerPoints.push_back(point);
      return centerPoints;
   }

   vector<std::array<Real, 3>> test_fp::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> centerPoints;

      creal dx = 0.0;
      creal dy = 0.0;
      creal dz = 0.0;

      return this->getV0(x,y,z,dx,dy,dz,popID);
   }

}// namespace projects
