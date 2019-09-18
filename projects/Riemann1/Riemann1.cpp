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
#include "../../object_wrapper.h"

#include "Riemann1.h"

using namespace std;

namespace projects {
   Riemann1::Riemann1(): Project() { }
   Riemann1::~Riemann1() { }

   bool Riemann1::initialize(void) {return Project::initialize();}

   void Riemann1::addParameters(){
      typedef Readparameters RP;
      RP::add("Riemann.rho1", "Number density, left state (m^-3)", 0.0);
      RP::add("Riemann.rho2", "Number density, right state (m^-3)", 0.0);
      RP::add("Riemann.T1", "Temperature, left state (K)", 0.0);
      RP::add("Riemann.T2", "Temperature, right state (K)", 0.0);
      RP::add("Riemann.Vx1", "Bulk velocity x component, left state (m/s)", 0.0);
      RP::add("Riemann.Vx2", "Bulk velocity x component, right state (m/s)", 0.0);
      RP::add("Riemann.Vy1", "Bulk velocity y component, left state (m/s)", 0.0);
      RP::add("Riemann.Vy2", "Bulk velocity y component, right state (m/s)", 0.0);
      RP::add("Riemann.Vz1", "Bulk velocity z component, left state (m/s)", 0.0);
      RP::add("Riemann.Vz2", "Bulk velocity z component, right state (m/s)", 0.0);
      RP::add("Riemann.Bx1", "Magnetic field x component, left state (T)", 0.0);
      RP::add("Riemann.Bx2", "Magnetic field x component, right state (T)", 0.0);
      RP::add("Riemann.By1", "Magnetic field y component, left state (T)", 0.0);
      RP::add("Riemann.By2", "Magnetic field y component, right state (T)", 0.0);
      RP::add("Riemann.Bz1", "Magnetic field z component, left state (T)", 0.0);
      RP::add("Riemann.Bz2", "Magnetic field z component, right state (T)", 0.0);
      RP::add("Riemann.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Riemann.nVelocitySamples", "Number of sampling points per velocity dimension", 5);      
   }

   void Riemann1::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("Riemann.rho1", this->rho[this->LEFT]);
      RP::get("Riemann.rho2", this->rho[this->RIGHT]);
      RP::get("Riemann.T1", this->T[this->LEFT]);
      RP::get("Riemann.T2", this->T[this->RIGHT]);
      RP::get("Riemann.Vx1", this->Vx[this->LEFT]);
      RP::get("Riemann.Vx2", this->Vx[this->RIGHT]);
      RP::get("Riemann.Vy1", this->Vy[this->LEFT]);
      RP::get("Riemann.Vy2", this->Vy[this->RIGHT]);
      RP::get("Riemann.Vz1", this->Vz[this->LEFT]);
      RP::get("Riemann.Vz2", this->Vz[this->RIGHT]);
      RP::get("Riemann.Bx1", this->Bx[this->LEFT]);
      RP::get("Riemann.Bx2", this->Bx[this->RIGHT]);
      RP::get("Riemann.By1", this->By[this->LEFT]);
      RP::get("Riemann.By2", this->By[this->RIGHT]);
      RP::get("Riemann.Bz1", this->Bz[this->LEFT]);
      RP::get("Riemann.Bz2", this->Bz[this->RIGHT]);
      RP::get("Riemann.nSpaceSamples", this->nSpaceSamples);
      RP::get("Riemann.nVelocitySamples", this->nVelocitySamples);
   }

   Real Riemann1::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const uint popID) const {
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;
      
      return this->rho[side] * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T[side]), 1.5) *
      exp(- physicalconstants::MASS_PROTON * (pow(vx - this->Vx[side], 2.0) + pow(vy - this->Vy[side], 2.0) + pow(vz - this->Vz[side], 2.0)) / (2.0 * physicalconstants::K_B * this->T[side]));
   }

   Real Riemann1::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
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
                     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz, popID);
                  }
      return avg / pow(this->nSpaceSamples, 3.0) / pow(this->nVelocitySamples, 3.0);
   }

   void Riemann1::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }
   
   void Riemann1::setProjectBField(
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
                  
                  Real Bxavg, Byavg, Bzavg;
                  Bxavg = Byavg = Bzavg = 0.0;
                  if(this->nSpaceSamples > 1) {
                     Real d_x = perBGrid.DX / (this->nSpaceSamples - 1);
                     Real d_z = perBGrid.DZ / (this->nSpaceSamples - 1);
                     for (uint i=0; i<this->nSpaceSamples; ++i) {
                        for (uint k=0; k<this->nSpaceSamples; ++k) {
                           Bxavg += ((xyz[0] + i * d_x) < 0.0) ? this->Bx[this->LEFT] : this->Bx[this->RIGHT];
                           Byavg += ((xyz[0] + i * d_x) < 0.0) ? this->By[this->LEFT] : this->By[this->RIGHT];
                           Bzavg += ((xyz[0] + i * d_x) < 0.0) ? this->Bz[this->LEFT] : this->Bz[this->RIGHT];
                        }
                     }
                     cuint nPts = pow(this->nSpaceSamples, 3.0);
                     
                     cell->at(fsgrids::bfield::PERBX) = Bxavg / nPts;
                     cell->at(fsgrids::bfield::PERBY) = Byavg / nPts;
                     cell->at(fsgrids::bfield::PERBZ) = Bzavg / nPts;
                  } else {
                     cell->at(fsgrids::bfield::PERBX) = (xyz[0] < 0.0) ? this->Bx[this->LEFT] : this->Bx[this->RIGHT];
                     cell->at(fsgrids::bfield::PERBY) = (xyz[0] < 0.0) ? this->By[this->LEFT] : this->By[this->RIGHT];
                     cell->at(fsgrids::bfield::PERBZ) = (xyz[0] < 0.0) ? this->Bz[this->LEFT] : this->Bz[this->RIGHT];
                  }
               }
            }
         }
      }
   }
}
