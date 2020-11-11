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
         RP::add(pop + "_Harris.nSpaceSamples", "Number of sampling points per spatial dimension.", 2);
         RP::add(pop + "_Harris.nVelocitySamples", "Number of sampling points per velocity dimension.", 2);
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
         RP::get(pop + "_Harris.nSpaceSamples", sP.nSpaceSamples);
         RP::get(pop + "_Harris.nVelocitySamples", sP.nVelocitySamples);

         speciesParams.push_back(sP);
      }
   }
   
   Real Harris::getDistribValue(
      creal& x,creal& y, creal& z,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz,
      const uint popID
   ) const {

      const HarrisSpeciesParameters& sP = speciesParams[popID];
      Real mass = getObjectWrapper().particleSpecies[popID].mass;

      return sP.DENSITY * pow(mass / (2.0 * M_PI * physicalconstants::K_B * sP.TEMPERATURE), 1.5) * (
         5.0 / pow(cosh(x / (this->SCA_LAMBDA)), 2.0) * exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * physicalconstants::K_B * sP.TEMPERATURE))
         +
         exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * physicalconstants::K_B * sP.TEMPERATURE)));
   }
   
   Real Harris::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz,const uint popID
   ) const {
      const HarrisSpeciesParameters& sP = speciesParams[popID];
      if((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
         creal d_x = dx / (sP.nSpaceSamples-1);
         creal d_y = dy / (sP.nSpaceSamples-1);
         creal d_z = dz / (sP.nSpaceSamples-1);
         creal d_vx = dvx / (sP.nVelocitySamples-1);
         creal d_vy = dvy / (sP.nVelocitySamples-1);
         creal d_vz = dvz / (sP.nVelocitySamples-1);
         
         Real avg = 0.0;
         // #pragma omp parallel for collapse(6) reduction(+:avg)
         // WARNING No threading here if calling functions are already threaded
         for (uint i=0; i<sP.nSpaceSamples; ++i)
            for (uint j=0; j<sP.nSpaceSamples; ++j)
               for (uint k=0; k<sP.nSpaceSamples; ++k)
                  for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                     for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                        for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
                           avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz, popID);
                        }
         return avg /
         (sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples) /
         (sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz, popID);
      }
      
      
      
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
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);
      
      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();
         
         #pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
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
