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
#include "../../backgroundfield/backgroundfield.h"

#include "Shock.h"

namespace projects {
   Shock::Shock(): Project() { }
   Shock::~Shock() { }

   bool Shock::initialize(void) {return Project::initialize();}

   void Shock::addParameters() {
      typedef Readparameters RP;
      RP::add("Shock.BX0", "Background field value (T)", 1.0e-9);
      RP::add("Shock.BY0", "Background field value (T)", 2.0e-9);
      RP::add("Shock.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Shock.EX0", "Background electric field", 0.0);
      RP::add("Shock.VX0", "Bulk velocity in x", 0.0);
      RP::add("Shock.VY0", "Bulk velocity in y", 0.0);
      RP::add("Shock.VZ0", "Bulk velocuty in z", 0.0);
      RP::add("Shock.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Shock.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Shock.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
      RP::add("Shock.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
      RP::add("Shock.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
      RP::add("Shock.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Shock.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      RP::add("Shock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("Shock.Scale_x", "Scale length in x (m)", 2.0e6);
      RP::add("Shock.Scale_y", "Scale length in y (m)", 2.0e6);
      RP::add("Shock.Sharp_Y", "Sharpness of tannh", 0.1);
   }

   void Shock::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("Shock.BX0", this->BX0);
      RP::get("Shock.BY0", this->BY0);
      RP::get("Shock.BZ0", this->BZ0);
      RP::get("Shock.EX0", this->EX0);
      RP::get("Shock.VX0", this->VX0);
      RP::get("Shock.VY0", this->VY0);
      RP::get("Shock.VZ0", this->VZ0);
      RP::get("Shock.rho", this->DENSITY);
      RP::get("Shock.Temperature", this->TEMPERATURE);
      RP::get("Shock.magPertAmp", this->magPertAmp);
      RP::get("Shock.densityPertAmp", this->densityPertAmp);
      RP::get("Shock.velocityPertAmp", this->velocityPertAmp);
      RP::get("Shock.nSpaceSamples", this->nSpaceSamples);
      RP::get("Shock.nVelocitySamples", this->nVelocitySamples);
      RP::get("Shock.maxwCutoff", this->maxwCutoff);
      RP::get("Shock.Scale_x", this->SCA_X);
      RP::get("Shock.Scale_y", this->SCA_Y);
      RP::get("Shock.Sharp_Y", this->Sharp_Y);
   }

   Real Shock::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, const uint popID) const {
      creal kb = physicalconstants::K_B;
      creal mass = physicalconstants::MASS_PROTON;
      return exp(- mass * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0)+ (vz-this->VZ0)*(vz-this->VZ0)) / (2.0 * kb * this->TEMPERATURE));
      //*exp(-pow(x-Parameters::xmax/2.0, 2.0)/pow(this->SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/4.0, 2.0)/pow(this->SCA_Y, 2.0));
   }

   Real Shock::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz,
           creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
      vmesh::MeshParameters& meshParams = getObjectWrapper().velocityMeshes[meshID];
      if (vx < meshParams.meshMinLimits[0] + 0.5*dvx ||
          vy < meshParams.meshMinLimits[1] + 0.5*dvy ||
          vz < meshParams.meshMinLimits[2] + 0.5*dvz ||
          vx > meshParams.meshMaxLimits[0] - 1.5*dvx ||
          vy > meshParams.meshMaxLimits[1] - 1.5*dvy ||
          vz > meshParams.meshMaxLimits[2] - 1.5*dvz) {
         return 0.0;
      }
      
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      creal d_x = dx / (this->nSpaceSamples-1);
      creal d_y = dy / (this->nSpaceSamples-1);
      creal d_z = dz / (this->nSpaceSamples-1);
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      
      for (uint i=0; i<this->nSpaceSamples; ++i)
      for (uint j=0; j<this->nSpaceSamples; ++j)
      for (uint k=0; k<this->nSpaceSamples; ++k)      
      for (uint vi=0; vi<this->nVelocitySamples; ++vi)
         for (uint vj=0; vj<this->nVelocitySamples; ++vj)
      for (uint vk=0; vk<this->nVelocitySamples; ++vk)
            {
         avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, popID);
         }
      
      creal result = avg *this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) /
                     (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) / 
   //                (Parameters::vzmax - Parameters::vzmin) / 
                     (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
               
               
      if(result < this->maxwCutoff) {
         return 0.0;
      } else {
         return result;
      }
   }

   void Shock::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }
   
   void Shock::setProjectBField(
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
                  
                  cell->at(fsgrids::bfield::PERBX) = 0.0;
                  cell->at(fsgrids::bfield::PERBY) = 0.0;
                  cell->at(fsgrids::bfield::PERBZ) = this->BZ0*(3.0 + 2.0*tanh((xyz[1] - Parameters::ymax/2.0)/(this->Sharp_Y*Parameters::ymax)));
               }
            }
         }
      }
   }
}//namespace projects
