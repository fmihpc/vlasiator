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
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../object_wrapper.h"

#include "Distributions.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Distributions::Distributions(): TriAxisSearch() { }
   Distributions::~Distributions() { }


   bool Distributions::initialize(void) {return Project::initialize();}

   void Distributions::addParameters(){
      typedef Readparameters RP;
      RP::add("Distributions.rho1", "Number density, first peak (m^-3)", 0.0);
      RP::add("Distributions.rho2", "Number density, second peak (m^-3)", 0.0);
      RP::add("Distributions.Tx1", "Temperature, first peak (K)", 0.0);
      RP::add("Distributions.Tx2", "Temperature, second peak (K)", 0.0);
      RP::add("Distributions.Ty1", "Temperature, first peak (K)", 0.0);
      RP::add("Distributions.Ty2", "Temperature, second peak (K)", 0.0);
      RP::add("Distributions.Tz1", "Temperature, first peak (K)", 0.0);
      RP::add("Distributions.Tz2", "Temperature, second peak (K)", 0.0);
      RP::add("Distributions.Vx1", "Bulk velocity x component, first peak (m/s)", 0.0);
      RP::add("Distributions.Vx2", "Bulk velocity x component, second peak (m/s)", 0.0);
      RP::add("Distributions.Vy1", "Bulk velocity y component, first peak (m/s)", 0.0);
      RP::add("Distributions.Vy2", "Bulk velocity y component, second peak (m/s)", 0.0);
      RP::add("Distributions.Vz1", "Bulk velocity z component, first peak (m/s)", 0.0);
      RP::add("Distributions.Vz2", "Bulk velocity z component, second peak (m/s)", 0.0);
      RP::add("Distributions.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("Distributions.By", "Magnetic field y component (T)", 0.0);
      RP::add("Distributions.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("Distributions.dBx", "Magnetic field x component cosine perturbation amplitude (T)", 0.0);
      RP::add("Distributions.dBy", "Magnetic field y component cosine perturbation amplitude (T)", 0.0);
      RP::add("Distributions.dBz", "Magnetic field z component cosine perturbation amplitude (T)", 0.0);
      RP::add("Distributions.magXPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along x (T)", 1.0e-9);
      RP::add("Distributions.magYPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along y (T)", 1.0e-9);
      RP::add("Distributions.magZPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along z (T)", 1.0e-9);
      RP::add("Distributions.rho1PertAbsAmp", "Absolute amplitude of the density perturbation, first peak", 0.1);
      RP::add("Distributions.rho2PertAbsAmp", "Absolute amplitude of the density perturbation, second peak", 0.1);
//       RP::add("Distributions.Vx1PertAbsAmp", "Absolute amplitude of the Vx perturbation, first peak", 1.0e6);
//       RP::add("Distributions.Vy1PertAbsAmp", "Absolute amplitude of the Vy perturbation, first peak", 1.0e6);
//       RP::add("Distributions.Vz1PertAbsAmp", "Absolute amplitude of the Vz perturbation, first peak", 1.0e6);
//       RP::add("Distributions.Vx2PertAbsAmp", "Absolute amplitude of the Vx perturbation, second peak", 1.0e6);
//       RP::add("Distributions.Vy2PertAbsAmp", "Absolute amplitude of the Vy perturbation, second peak", 1.0e6);
//       RP::add("Distributions.Vz2PertAbsAmp", "Absolute amplitude of the Vz perturbation, second peak", 1.0e6);
      RP::add("Distributions.lambda", "B cosine perturbation wavelength (m)", 0.0);
   }

   void Distributions::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;
      Project::getParameters();

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }

      RP::get("Distributions.rho1", this->rho[0]);
      RP::get("Distributions.rho2", this->rho[1]);
      RP::get("Distributions.Tx1", this->Tx[0]);
      RP::get("Distributions.Tx2", this->Tx[1]);
      RP::get("Distributions.Ty1", this->Ty[0]);
      RP::get("Distributions.Ty2", this->Ty[1]);
      RP::get("Distributions.Tz1", this->Tz[0]);
      RP::get("Distributions.Tz2", this->Tz[1]);
      RP::get("Distributions.Vx1", this->Vx[0]);
      RP::get("Distributions.Vx2", this->Vx[1]);
      RP::get("Distributions.Vy1", this->Vy[0]);
      RP::get("Distributions.Vy2", this->Vy[1]);
      RP::get("Distributions.Vz1", this->Vz[0]);
      RP::get("Distributions.Vz2", this->Vz[1]);
      RP::get("Distributions.Bx", this->Bx);
      RP::get("Distributions.By", this->By);
      RP::get("Distributions.Bz", this->Bz);
      RP::get("Distributions.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("Distributions.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("Distributions.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("Distributions.rho1PertAbsAmp", this->rhoPertAbsAmp[0]);
      RP::get("Distributions.rho2PertAbsAmp", this->rhoPertAbsAmp[1]);
//       RP::get("Distributions.Vx1PertAbsAmp", this->Vx1PertAbsAmp);
//       RP::get("Distributions.Vy1PertAbsAmp", this->Vy1PertAbsAmp);
//       RP::get("Distributions.Vz1PertAbsAmp", this->Vz1PertAbsAmp);
//       RP::get("Distributions.Vx2PertAbsAmp", this->Vx2PertAbsAmp);
//       RP::get("Distributions.Vy2PertAbsAmp", this->Vy2PertAbsAmp);
//       RP::get("Distributions.Vz2PertAbsAmp", this->Vz2PertAbsAmp);
      RP::get("Distributions.dBx", this->dBx);
      RP::get("Distributions.dBy", this->dBy);
      RP::get("Distributions.dBz", this->dBz);
      RP::get("Distributions.lambda", this->lambda);
   }

   Realf Distributions::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];
      creal relx = x/(Parameters::xmax - Parameters::xmin);
      creal rely = y/(Parameters::ymax - Parameters::ymin);
      creal relz = z/(Parameters::zmax - Parameters::zmin);
      creal scaledVx1 = this->Vx[1] * relx;
      creal scaledVy1 = this->Vy[1] * rely;
      creal scaledVz1 = this->Vz[1] * relz;

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real initRho0 = this->rhoRnd[0];
      const Real initRho1 = this->rhoRnd[1];
      const Real initT0x = this->Tx[0];
      const Real initT0y = this->Ty[0];
      const Real initT0z = this->Tz[0];
      const Real initT1x = this->Tx[1];
      const Real initT1y = this->Ty[1];
      const Real initT1z = this->Tz[1];
      const Real initV0X = this->Vx[0];
      const Real initV0Y = this->Vy[0];
      const Real initV0Z = this->Vz[0];
      const Real initV1X = this->Vx[1];
      const Real initV1Y = this->Vy[1];
      const Real initV1Z = this->Vz[1];

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
               Real vx = vxBlock + (i+0.5)*dvxCell - initV0X;
               Real vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
               Real vz = vzBlock + (k+0.5)*dvzCell - initV0Z;
               Realf value = TriMaxwellianPhaseSpaceDensity(vx,vy,vz,initT0x,initT0y,initT0z,initRho0,mass);
               vx = vxBlock + (i+0.5)*dvxCell - initV1X;
               vy = vyBlock + (j+0.5)*dvyCell - initV1Y;
               vz = vzBlock + (k+0.5)*dvzCell - initV1Z;
               value += TriMaxwellianPhaseSpaceDensity(vx,vy,vz,initT1x,initT1y,initT1z,initRho1,mass);
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
   Realf Distributions::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];
      creal relx = x/(Parameters::xmax - Parameters::xmin);
      creal rely = y/(Parameters::ymax - Parameters::ymin);
      creal relz = z/(Parameters::zmax - Parameters::zmin);
      creal scaledVx1 = this->Vx[1] * relx;
      creal scaledVy1 = this->Vy[1] * rely;
      creal scaledVz1 = this->Vz[1] * relz;

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real initRho0 = this->rhoRnd[0];
      const Real initRho1 = this->rhoRnd[1];
      const Real initT0x = this->Tx[0];
      const Real initT0y = this->Ty[0];
      const Real initT0z = this->Tz[0];
      const Real initT1x = this->Tx[1];
      const Real initT1y = this->Ty[1];
      const Real initT1z = this->Tz[1];
      const Real initV0X = this->Vx[0];
      const Real initV0Y = this->Vy[0];
      const Real initV0Z = this->Vz[0];
      const Real initV1X = this->Vx[1];
      const Real initV1Y = this->Vy[1];
      const Real initV1Z = this->Vz[1];

      Realf value = 0;
      Real vx = vx_in - initV0X;
      Real vy = vy_in - initV0Y;
      Real vz = vz_in - initV0Z;
      value = TriMaxwellianPhaseSpaceDensity(vx,vy,vz,initT0x,initT0y,initT0z,initRho0,mass);
      vx = vx_in - initV1X;
      vy = vy_in - initV1Y;
      vz = vz_in - initV1Z;
      value += TriMaxwellianPhaseSpaceDensity(vx,vy,vz,initT1x,initT1y,initT1z,initRho1,mass);
      return value;
   }

   void Distributions::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      std::default_random_engine rndState;
      setRandomCellSeed(cell,rndState);
      for (uint i=0; i<2; i++) {
         this->rhoRnd[i] = this->rho[i] + this->rhoPertAbsAmp[i] * (0.5 - getRandomNumber(rndState));
      }
   }

   void Distributions::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);

      setBackgroundField(bgField, BgBGrid);

      if(!P::isRestart) {
         const auto localSize = BgBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);

                  std::default_random_engine rndState;
                  setRandomSeed(cellid,rndState);

                  if (this->lambda != 0.0) {
                     cell->at(fsgrids::bfield::PERBX) = this->dBx*cos(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBY) = this->dBy*sin(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBZ) = this->dBz*cos(2.0 * M_PI * xyz[0] / this->lambda);
                  }

                  cell->at(fsgrids::bfield::PERBX) += this->magXPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBY) += this->magYPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBZ) += this->magZPertAbsAmp * (0.5 - getRandomNumber(rndState));
               }
            }
         }
      }
   }

   vector<std::array<Real, 3>> Distributions::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> centerPoints;
      creal relx = x/(Parameters::xmax - Parameters::xmin);
      creal rely = y/(Parameters::ymax - Parameters::ymin);
      creal relz = z/(Parameters::zmax - Parameters::zmin);
      creal scaledVx1 = this->Vx[1] * relx;
      creal scaledVy1 = this->Vy[1] * rely;
      creal scaledVz1 = this->Vz[1] * relz;
      std::array<Real, 3> point0 {{this->Vx[0], this->Vy[0], this->Vz[0]}};
      std::array<Real, 3> point1 {{scaledVx1, scaledVy1, scaledVz1}};
      centerPoints.push_back(point0);
      centerPoints.push_back(point1);
      return centerPoints;
   }

}// namespace projects
