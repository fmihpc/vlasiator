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

   Real Distributions::getDistribValue(
      creal& x, creal& y, creal& z,
      creal& vx, creal& vy, creal& vz,
      const uint popID
   ) const {
      Real value = 0.0;
      creal relx = x/(Parameters::xmax - Parameters::xmin);
      creal rely = y/(Parameters::ymax - Parameters::ymin);
      creal relz = z/(Parameters::zmax - Parameters::zmin);
      creal scaledVx1 = this->Vx[1] * relx;
      creal scaledVy1 = this->Vy[1] * rely;
      creal scaledVz1 = this->Vz[1] * relz;
      
      value += this->rhoRnd[0] * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B ), 1.5) * 1.0 / sqrt(this->Tx[0]*this->Ty[0]*this->Tz[0]) *
      exp(-physicalconstants::MASS_PROTON * (pow(vx - this->Vx[0], 2.0) / (2.0 * physicalconstants::K_B * this->Tx[0]) + pow(vy - this->Vy[0], 2.0) / (2.0 * physicalconstants::K_B * this->Ty[0]) + pow(vz - this->Vz[0], 2.0) / (2.0 * physicalconstants::K_B * this->Tz[0])));
      
      value += this->rhoRnd[1] * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B ), 1.5) * 1.0 / sqrt(this->Tx[1]*this->Ty[1]*this->Tz[1]) *
      exp(-physicalconstants::MASS_PROTON * (pow(vx - scaledVx1, 2.0) / (2.0 * physicalconstants::K_B * this->Tx[1]) + pow(vy - scaledVy1, 2.0) / (2.0 * physicalconstants::K_B * this->Ty[1]) + pow(vz - scaledVz1, 2.0) / (2.0 * physicalconstants::K_B * this->Tz[1])));
      
      return value;
   }

   Real Distributions::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {   
      return getDistribValue(x+0.5*dx, y+0.5*dy,z+0.5*dz,vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, popID);
   }

   void Distributions::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      setRandomCellSeed(cell);
      for (uint i=0; i<2; i++) {
         this->rhoRnd[i] = this->rho[i] + this->rhoPertAbsAmp[i] * (0.5 - getRandomNumber());
      }
   }

   void Distributions::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField, BgBGrid);
      
      if(!P::isRestart) {
         const auto localSize = BgBGrid.getLocalSize().data();
         
#pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  
                  setRandomSeed(cellid);
                  
                  if (this->lambda != 0.0) {
                     cell->at(fsgrids::bfield::PERBX) = this->dBx*cos(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBY) = this->dBy*sin(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBZ) = this->dBz*cos(2.0 * M_PI * xyz[0] / this->lambda);
                  }

                  cell->at(fsgrids::bfield::PERBX) += this->magXPertAbsAmp * (0.5 - getRandomNumber());
                  cell->at(fsgrids::bfield::PERBY) += this->magYPertAbsAmp * (0.5 - getRandomNumber());
                  cell->at(fsgrids::bfield::PERBZ) += this->magZPertAbsAmp * (0.5 - getRandomNumber());
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
