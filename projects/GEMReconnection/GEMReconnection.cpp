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

#include "GEMReconnection.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   GEMReconnection::GEMReconnection(): TriAxisSearch() { }
   GEMReconnection::~GEMReconnection() { }
   
   bool GEMReconnection::initialize(void) {return Project::initialize();}
   
   void GEMReconnection::addParameters(){
      typedef Readparameters RP;
      RP::add<Real>("GEMReconnection.Scale_size", "GEM reconnection challenge current sheet scale size (m)", this->SCA_LAMBDA,150000.0);
      RP::add<Real>("GEMReconnection.BX0", "Magnetic field at infinity (T)", this->BX0,1e-8);
      RP::add<Real>("GEMReconnection.BY0", "Magnetic field at infinity (T)", this->BY0,0.0);
      RP::add<Real>("GEMReconnection.BZ0", "Magnetic field at infinity (T)", this->BZ0,0.0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         GEMReconnectionSpeciesParameters* sP=new GEMReconnectionSpeciesParameters();
         speciesParamsRead.push_back(sP);

         RP::add<Real>(pop + "_GEMReconnection.Temperature", "Temperature (K)", sP->TEMPERATURE,2.0e6);
         RP::add<Real>(pop + "_GEMReconnection.rho", "Number density at infinity (m^-3)", sP->DENSITY,1.0e7);
      }
   }
   
   void GEMReconnection::getParameters(){
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        this->speciesParams.push_back(*this->speciesParamsRead.at(i));
      }
   }

   Realf GEMReconnection::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const GEMReconnectionSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY / pow(cosh(z / (this->SCA_LAMBDA)), 2.0) + sP.DENSITY * 0.2;
      Real initT = sP.TEMPERATURE;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

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
   Realf GEMReconnection::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      const GEMReconnectionSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY / pow(cosh(z / (this->SCA_LAMBDA)), 2.0) + sP.DENSITY * 0.2;
      Real initT = sP.TEMPERATURE;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
      return value;
   }

   void GEMReconnection::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   vector<std::array<Real, 3>> GEMReconnection::getV0(
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

   void GEMReconnection::setProjectBField(fsgrids::perbspan perb,
                                 fsgrids::bgbspan bgb,
                                 fsgrids::technicalspan technical, FieldSolverGrid &fsgrid) {
      setBackgroundFieldToZero(fsgrid, technical, bgb);
      const GEMReconnectionSpeciesParameters& sP = speciesParams[0];

      creal Lx = Parameters::xmax - Parameters::xmin;
      creal Ly = Parameters::ymax - Parameters::ymin;
      creal Lz = Parameters::zmax - Parameters::zmin;

      creal rho0 = sP.DENSITY;

      creal di = 299792458.0/sqrt(rho0*physicalconstants::CHARGE*physicalconstants::CHARGE/physicalconstants::MASS_PROTON/physicalconstants::EPS_0);

      if(!P::isRestart) {
         // local copies for lambda capture
         const auto BX0_l = this->BX0;
         const auto BY0_l = this->BY0;
         const auto BZ0_l = this->BZ0;
         const auto SCA_LAMBDA_l = this->SCA_LAMBDA;

         fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                             phiprof::initializeTimer("setProjectBField-loop"), technical,
                             [=](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
            const std::array<Real, 3> xyz = coordinates.getPhysicalCoords(stencil.i, stencil.j, stencil.k);
            const std::array<Real, 3> gridSpacing = coordinates.physicalGridSpacing;
            auto& cell = perb[stencil.ooo()];

            Real Bx_island, By_island, Bz_island;

            Bx_island = -M_PI * di * BX0_l * 0.1 * cos(2.0 * M_PI * (xyz[0] + 0.5 * gridSpacing[0]) / Lx) * sin(M_PI * (xyz[2] + 0.5 * gridSpacing[2]) / Lz) / Lz;
            Bz_island = 2.0 * M_PI * di * BX0_l * 0.1 * sin(2.0 * M_PI * (xyz[0] + 0.5 * gridSpacing[0]) / Lx) * cos(M_PI * (xyz[2] + 0.5 * gridSpacing[2]) / Lz) / Lx;

            cell[fsgrids::bfield::PERBX] = BX0_l * tanh((xyz[2] + 0.5 * gridSpacing[2]) / SCA_LAMBDA_l) + Bx_island;
            cell[fsgrids::bfield::PERBY] = 0.0;
            cell[fsgrids::bfield::PERBZ] = Bz_island;
         });
      }
   }

} // namespace projects
