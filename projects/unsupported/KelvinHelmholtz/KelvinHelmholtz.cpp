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

#include "KelvinHelmholtz.h"

namespace projects {
   using namespace std;
   KelvinHelmholtz::KelvinHelmholtz(): Project() { }
   KelvinHelmholtz::~KelvinHelmholtz() { }
   
   bool KelvinHelmholtz::initialize(void) {return true;}
   
   void KelvinHelmholtz::addParameters() {
      typedef Readparameters RP;
      RP::add("KelvinHelmholtz.rho1", "Number density, this->TOP state (m^-3)", 0.0);
      RP::add("KelvinHelmholtz.rho2", "Number density, this->BOTTOM state (m^-3)", 0.0);
      RP::add("KelvinHelmholtz.T1", "Temperature, this->TOP state (K)", 0.0);
      RP::add("KelvinHelmholtz.T2", "Temperature, this->BOTTOM state (K)", 0.0);
      RP::add("KelvinHelmholtz.Vx1", "Bulk velocity x component, this->TOP state (m/s)", 0.0);
      RP::add("KelvinHelmholtz.Vx2", "Bulk velocity x component, this->BOTTOM state (m/s)", 0.0);
      RP::add("KelvinHelmholtz.Vy1", "Bulk velocity y component, this->TOP state (m/s)", 0.0);
      RP::add("KelvinHelmholtz.Vy2", "Bulk velocity y component, this->BOTTOM state (m/s)", 0.0);
      RP::add("KelvinHelmholtz.Vz1", "Bulk velocity z component, this->TOP state (m/s)", 0.0);
      RP::add("KelvinHelmholtz.Vz2", "Bulk velocity z component, this->BOTTOM state (m/s)", 0.0);
      RP::add("KelvinHelmholtz.Bx1", "Magnetic field x component, this->TOP state (T)", 0.0);
      RP::add("KelvinHelmholtz.Bx2", "Magnetic field x component, this->BOTTOM state (T)", 0.0);
      RP::add("KelvinHelmholtz.By1", "Magnetic field y component, this->TOP state (T)", 0.0);
      RP::add("KelvinHelmholtz.By2", "Magnetic field y component, this->BOTTOM state (T)", 0.0);
      RP::add("KelvinHelmholtz.Bz1", "Magnetic field z component, this->TOP state (T)", 0.0);
      RP::add("KelvinHelmholtz.Bz2", "Magnetic field z component, this->BOTTOM state (T)", 0.0);
      RP::add("KelvinHelmholtz.lambda", "Initial perturbation wavelength (m)", 0.0);
      RP::add("KelvinHelmholtz.amp", "Initial perturbation amplitude (m)", 0.0);
      RP::add("KelvinHelmholtz.offset", "Boundaries offset from 0 (m)", 0.0);
      RP::add("KelvinHelmholtz.transitionWidth", "Width of tanh transition for all changing values", 0.0);
      RP::add("KelvinHelmholtz.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("KelvinHelmholtz.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }

   void KelvinHelmholtz::getParameters() {
      typedef Readparameters RP;
      RP::get("KelvinHelmholtz.rho1", this->rho[this->TOP]);
      RP::get("KelvinHelmholtz.rho2", this->rho[this->BOTTOM]);
      RP::get("KelvinHelmholtz.T1", this->T[this->TOP]);
      RP::get("KelvinHelmholtz.T2", this->T[this->BOTTOM]);
      RP::get("KelvinHelmholtz.Vx1", this->Vx[this->TOP]);
      RP::get("KelvinHelmholtz.Vx2", this->Vx[this->BOTTOM]);
      RP::get("KelvinHelmholtz.Vy1", this->Vy[this->TOP]);
      RP::get("KelvinHelmholtz.Vy2", this->Vy[this->BOTTOM]);
      RP::get("KelvinHelmholtz.Vz1", this->Vz[this->TOP]);
      RP::get("KelvinHelmholtz.Vz2", this->Vz[this->BOTTOM]);
      RP::get("KelvinHelmholtz.Bx1", this->Bx[this->TOP]);
      RP::get("KelvinHelmholtz.Bx2", this->Bx[this->BOTTOM]);
      RP::get("KelvinHelmholtz.By1", this->By[this->TOP]);
      RP::get("KelvinHelmholtz.By2", this->By[this->BOTTOM]);
      RP::get("KelvinHelmholtz.Bz1", this->Bz[this->TOP]);
      RP::get("KelvinHelmholtz.Bz2", this->Bz[this->BOTTOM]);
      RP::get("KelvinHelmholtz.lambda", this->lambda);
      RP::get("KelvinHelmholtz.amp", this->amp);
      RP::get("KelvinHelmholtz.offset", this->offset);
      RP::get("KelvinHelmholtz.transitionWidth", this->transitionWidth);
      RP::get("KelvinHelmholtz.nSpaceSamples", this->nSpaceSamples);
      RP::get("KelvinHelmholtz.nVelocitySamples", this->nVelocitySamples);
   }
   
   
   Real KelvinHelmholtz::profile(creal top, creal bottom, creal x, creal z) {
      if(top == bottom) {
         return top;
      }
      if(this->offset != 0.0) {
         return 0.5 * ((top-bottom) * (
         tanh((z + this->offset + this->amp * cos(2.0*M_PI*x/this->lambda))/this->transitionWidth) -
         tanh((z-(this->offset + this->amp * cos(2.0*M_PI*x/this->lambda)))/this->transitionWidth) -1) + top+bottom);
      } else {
         return 0.5 * ((top-bottom) * tanh(z/this->transitionWidth) + top+bottom);
      }
   }
   
   Real KelvinHelmholtz::getDistribValue(creal& x, creal& z, creal& vx, creal& vy, creal& vz){
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      Real rho = profile(this->rho[this->BOTTOM], this->rho[this->TOP], x, z);
      Real T = profile(this->T[this->BOTTOM], this->T[this->TOP], x, z);
      Real Vx = profile(this->Vx[this->BOTTOM], this->Vx[this->TOP], x, z);
      Real Vy = profile(this->Vy[this->BOTTOM], this->Vy[this->TOP], x, z);
      Real Vz = profile(this->Vz[this->BOTTOM], this->Vz[this->TOP], x, z);
      
      return rho * pow(mass / (2.0 * M_PI * kb * T), 1.5) *
      exp(- mass * (pow(vx - Vx, 2.0) + pow(vy - Vy, 2.0) + pow(vz - Vz, 2.0)) / (2.0 * kb * T));
   }

   Real KelvinHelmholtz::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
      creal d_x = dx / (this->nSpaceSamples-1);
      creal d_z = dz / (this->nSpaceSamples-1);
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      uint samples=0;

      Real middleValue=getDistribValue(x+0.5*dx, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz);
      #warning TODO: add SpatialCell::velocity_block_threshold() in place of sparseMinValue
      if(middleValue<0.000001*Parameters::sparseMinValue){
         return middleValue; //abort, this will not be accepted anyway
      }
      
   //#pragma omp parallel for collapse(6) reduction(+:avg)
      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint k=0; k<this->nSpaceSamples; ++k)
            for (uint vi=0; vi<this->nVelocitySamples; ++vi)
               for (uint vj=0; vj<this->nVelocitySamples; ++vj)
                  for (uint vk=0; vk<this->nVelocitySamples; ++vk){
                     avg +=getDistribValue(x+i*d_x, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
                     samples++;
                  }
      return avg / samples;
   }
   

   void KelvinHelmholtz::calcCellParameters(Real* cellParams,creal& t) {
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
   }
   
   void KelvinHelmholtz::setCellBackgroundField(SpatialCell* cell) {
      creal x = cell->parameters[CellParams::XCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dz = cell->parameters[CellParams::DZ];
      
      Real Bxavg, Byavg, Bzavg;
      Bxavg = Byavg = Bzavg = 0.0;
      Real d_x = dx / (this->nSpaceSamples - 1);
      Real d_z = dz / (this->nSpaceSamples - 1);
      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint k=0; k<this->nSpaceSamples; ++k) {
            Bxavg += profile(this->Bx[this->BOTTOM], this->Bx[this->TOP], x+i*d_x, z+k*d_z);
            Byavg += profile(this->By[this->BOTTOM], this->By[this->TOP], x+i*d_x, z+k*d_z);
            Bzavg += profile(this->Bz[this->BOTTOM], this->Bz[this->TOP], x+i*d_x, z+k*d_z);
         }
      cuint nPts = pow(this->nSpaceSamples, 2.0);
      
      cell->parameters[CellParams::BGBX   ] = Bxavg / nPts;
      cell->parameters[CellParams::BGBY   ] = Byavg / nPts;
      cell->parameters[CellParams::BGBZ   ] = Bzavg / nPts;
   }
} // namespace projects
