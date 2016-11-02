/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

#include "Firehose.h"

using namespace std;

namespace projects {
   Firehose::Firehose(): Project() { }
   Firehose::~Firehose() { }
   
   bool Firehose::initialize(void) {return Project::initialize();}

   void Firehose::addParameters(){
      typedef Readparameters RP;
      RP::add("Firehose.rho1", "Number density, first peak (m^-3)", 0.0);
      RP::add("Firehose.rho2", "Number density, second peak (m^-3)", 0.0);
      RP::add("Firehose.Tx1", "Temperature x, first peak (K)", 0.0);
      RP::add("Firehose.Tx2", "Temperature x, second peak (K)", 0.0);
      RP::add("Firehose.Ty1", "Temperature y, first peak (K)", 0.0);
      RP::add("Firehose.Ty2", "Temperature y, second peak (K)", 0.0);
      RP::add("Firehose.Tz1", "Temperature z, first peak (K)", 0.0);
      RP::add("Firehose.Tz2", "Temperature z, second peak (K)", 0.0);
      RP::add("Firehose.Vx1", "Bulk velocity x component, first peak (m/s)", 0.0);
      RP::add("Firehose.Vx2", "Bulk velocity x component, second peak (m/s)", 0.0);
      RP::add("Firehose.Vy1", "Bulk velocity y component, first peak (m/s)", 0.0);
      RP::add("Firehose.Vy2", "Bulk velocity y component, second peak (m/s)", 0.0);
      RP::add("Firehose.Vz1", "Bulk velocity z component, first peak (m/s)", 0.0);
      RP::add("Firehose.Vz2", "Bulk velocity z component, second peak (m/s)", 0.0);
      RP::add("Firehose.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("Firehose.By", "Magnetic field y component (T)", 0.0);
      RP::add("Firehose.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("Firehose.lambda", "Initial perturbation wavelength (m)", 0.0);
      RP::add("Firehose.amp", "Initial perturbation amplitude (m)", 0.0);
      RP::add("Firehose.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Firehose.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }

   void Firehose::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;
      RP::get("Firehose.rho1", this->rho[1]);
      RP::get("Firehose.rho2", this->rho[2]);
      RP::get("Firehose.Tx1", this->Tx[1]);
      RP::get("Firehose.Tx2", this->Tx[2]);
      RP::get("Firehose.Ty1", this->Ty[1]);
      RP::get("Firehose.Ty2", this->Ty[2]);
      RP::get("Firehose.Tz1", this->Tz[1]);
      RP::get("Firehose.Tz2", this->Tz[2]);
      RP::get("Firehose.Vx1", this->Vx[1]);
      RP::get("Firehose.Vx2", this->Vx[2]);
      RP::get("Firehose.Vy1", this->Vy[1]);
      RP::get("Firehose.Vy2", this->Vy[2]);
      RP::get("Firehose.Vz1", this->Vz[1]);
      RP::get("Firehose.Vz2", this->Vz[2]);
      RP::get("Firehose.Bx", this->Bx);
      RP::get("Firehose.By", this->By);
      RP::get("Firehose.Bz", this->Bz);
      RP::get("Firehose.lambda", this->lambda);
      RP::get("Firehose.amp", this->amp);
      RP::get("Firehose.nSpaceSamples", this->nSpaceSamples);
      RP::get("Firehose.nVelocitySamples", this->nVelocitySamples);
   }

   Real Firehose::profile(creal top, creal bottom, creal x) {
      return top * (1.0 + this->amp*cos(2.0*M_PI*x/this->lambda));
   }

   Real Firehose::getDistribValue(
      creal& x, creal& y,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz) {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      Real Vx = profile(this->Vx[1],this->Vx[1], x);
      
      return
      this->rho[1] * pow(mass / (2.0 * M_PI * kb * this->Tx[1]), 1.5) *
      exp(- mass * (pow(vx - Vx, 2.0) / (2.0 * kb * this->Tx[1]) +
                  pow(vy - this->Vy[1], 2.0) / (2.0 * kb * this->Ty[1]) +
               pow(vz - this->Vz[1], 2.0) / (2.0 * kb * this->Tz[1])));
   //   this->rho[2] * pow(mass / (2.0 * M_PI * kb * this->Tx[2]), 1.5) *
   //   exp(- mass * (pow(vx - this->Vx[2], 2.0) / (2.0 * kb * this->Tx[2]) + 
   //                 pow(vy - this->Vy[2], 2.0) / (2.0 * kb * this->Ty[2]) + 
   //           pow(vz - this->Vz[2], 2.0) / (2.0 * kb * this->Tz[2]))); 
   }

   Real Firehose::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) {
      creal d_x = dx / (this->nSpaceSamples-1);
      creal d_y = dy / (this->nSpaceSamples-1);
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
   //#pragma omp parallel for collapse(6) reduction(+:avg)
      for (uint i=0; i<this->nSpaceSamples; ++i)
      for (uint j=0; j<this->nSpaceSamples; ++j)
         for (uint vi=0; vi<this->nVelocitySamples; ++vi)
         for (uint vj=0; vj<this->nVelocitySamples; ++vj)
            for (uint vk=0; vk<this->nVelocitySamples; ++vk)
         {
            avg += getDistribValue(x+i*d_x, y+j*d_y, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
         }
      return avg / pow(this->nSpaceSamples, 2.0) /  pow(this->nVelocitySamples, 3.0);
   }

   void Firehose::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      cellParams[CellParams::PERBX   ] = this->Bx;
      cellParams[CellParams::PERBY   ] = this->By;
      cellParams[CellParams::PERBZ   ] = this->Bz;
   }
} // namespace projects
