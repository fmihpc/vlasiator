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
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../object_wrapper.h"


#include "MultiPeak.h"

using namespace std;
using namespace spatial_cell;


vector<Real> projects::MultiPeak::rhoRnd;

namespace projects {
   MultiPeak::MultiPeak(): TriAxisSearch() { }
   
   MultiPeak::~MultiPeak() { }

   bool MultiPeak::initialize(void) {
      return Project::initialize();
   }

   void MultiPeak::addParameters(){
      typedef Readparameters RP;
      RP::add("MultiPeak.n", "Number of populations to use", 0);
      RP::addComposing("MultiPeak.rho", "Number density (m^-3)");
      RP::addComposing("MultiPeak.Tx", "Temperature (K)");
      RP::addComposing("MultiPeak.Ty", "Temperature");
      RP::addComposing("MultiPeak.Tz", "Temperature");
      RP::addComposing("MultiPeak.Vx", "Bulk velocity x component (m/s)");
      RP::addComposing("MultiPeak.Vy", "Bulk velocity y component (m/s)");
      RP::addComposing("MultiPeak.Vz", "Bulk velocity z component (m/s)");
      RP::add("MultiPeak.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("MultiPeak.By", "Magnetic field y component (T)", 0.0);
      RP::add("MultiPeak.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("MultiPeak.dBx", "Magnetic field x component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBy", "Magnetic field y component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBz", "Magnetic field z component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.magXPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along x (T)", 1.0e-9);
      RP::add("MultiPeak.magYPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along y (T)", 1.0e-9);
      RP::add("MultiPeak.magZPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along z (T)", 1.0e-9);
      RP::addComposing("MultiPeak.rhoPertAbsAmp", "Absolute amplitude of the density perturbation");
      RP::add("MultiPeak.lambda", "B cosine perturbation wavelength (m)", 1.0);
      RP::add("MultiPeak.nVelocitySamples", "Number of sampling points per velocity dimension", 2);
      RP::add("MultiPeak.useMultipleSpecies","Is each peak a separate particle species",false);
      RP::add("MultiPeak.densityModel","Which spatial density model is used?",string("uniform"));
   }

   void MultiPeak::getParameters(){
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("MultiPeak.n", this->numberOfPopulations);
      if(this->numberOfPopulations < 1) {
         cerr << "You should set MultiPeak.n to more than 0 populations." << endl;
         abort();
      }
      RP::get("MultiPeak.rho", this->rho);
      RP::get("MultiPeak.Tx", this->Tx);
      RP::get("MultiPeak.Ty", this->Ty);
      RP::get("MultiPeak.Tz", this->Tz);
      RP::get("MultiPeak.Vx", this->Vx);
      RP::get("MultiPeak.Vy", this->Vy);
      RP::get("MultiPeak.Vz", this->Vz);
      RP::get("MultiPeak.Bx", this->Bx);
      RP::get("MultiPeak.By", this->By);
      RP::get("MultiPeak.Bz", this->Bz);
      RP::get("MultiPeak.rhoPertAbsAmp", this->rhoPertAbsAmp);
      if(!(this->rho.size() == this->Tx.size() &&
         this->Tx.size() == this->Ty.size() &&
         this->Ty.size() == this->Tz.size() &&
         this->Tz.size() == this->Vx.size() &&
         this->Vx.size() == this->Vy.size() &&
         this->Vy.size() == this->Vz.size() &&
         this->Vz.size() == this->rhoPertAbsAmp.size() &&
         this->rhoPertAbsAmp.size() == this->rho.size()
      )) {
         cerr << "You should define all parameters (MultiPeak.rho, MultiPeak.Tx, MultiPeak.Ty, MultiPeak.Tz, MultiPeak.Vx, MultiPeak.Vy, MultiPeak.Vz, MultiPeak.rhoPertAbsAmp) for each of the populations." << endl;
         abort();
      }
      if(this->numberOfPopulations > this->rho.size()) {
         cerr << "You are requesting more populations than are currently defined. Change MultiPeak.n or define more populations." << endl;
         abort();
      }
      
      RP::get("MultiPeak.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("MultiPeak.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("MultiPeak.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("MultiPeak.dBx", this->dBx);
      RP::get("MultiPeak.dBy", this->dBy);
      RP::get("MultiPeak.dBz", this->dBz);
      RP::get("MultiPeak.lambda", this->lambda);
      RP::get("MultiPeak.nVelocitySamples", this->nVelocitySamples);
      RP::get("MultiPeak.useMultipleSpecies", useMultipleSpecies);

      
      string densModelString;
      RP::get("MultiPeak.densityModel",densModelString);
      
      if (densModelString == "uniform") densityModel = Uniform;
      else if (densModelString == "testcase") densityModel = TestCase;
   }

   Real MultiPeak::getDistribValue(creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;

      Real value = 0.0;

      if (useMultipleSpecies == false) { // one species, multiple peaks
         if (popID != 0) return 0.0;

         for (uint i=0; i<this->numberOfPopulations; ++i) {
            value += this->rhoRnd[i]
                  * pow(mass / (2.0 * M_PI * kb ), 1.5) * 1.0
                  / sqrt(Tx[i]*Ty[i]*Tz[i]) 
                  * exp(- mass * (pow(vx - Vx[i], 2.0) / (2.0 * kb * Tx[i]) 
                                + pow(vy - Vy[i], 2.0) / (2.0 * kb * Ty[i]) 
                                + pow(vz - Vz[i], 2.0) / (2.0 * kb * Tz[i])));
         }
      } else { // multiple species, one peak each
         if (this->numberOfPopulations != getObjectWrapper().particleSpecies.size()) {
            cerr << "error number of peaks and populations do not match" << endl;
            exit(1);
         }

         value += this->rhoRnd[popID]
               * pow(mass / (2.0 * M_PI * kb ), 1.5)
               * 1.0 / sqrt(this->Tx[popID]*this->Ty[popID]*this->Tz[popID])
               * exp(- mass * (pow(vx - this->Vx[popID], 2.0) / (2.0 * kb * this->Tx[popID]) + pow(vy - this->Vy[popID], 2.0) / (2.0 * kb * this->Ty[popID]) + pow(vz - this->Vz[popID], 2.0) / (2.0 * kb * this->Tz[popID])));
      }
      return value;
   }

   Real MultiPeak::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, 
                                         creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                         const uint popID) const {
      // Iterative sampling of the distribution function. Keep track of the 
      // accumulated volume average over the iterations. When the next 
      // iteration improves the average by less than 1%, return the value.
      Real avgTotal = 0.0;
      bool ok = false;
      uint N = nVelocitySamples; // Start by using nVelocitySamples
      int N3_sum = 0;           // Sum of sampling points used so far
                                            
      #warning TODO: Replace getObjectWrapper().particleSpecies[popID].sparseMinValue with SpatialCell::velocity_block_threshold(?)
      const Real avgLimit = 0.01*getObjectWrapper().particleSpecies[popID].sparseMinValue;
      do {
         Real avg = 0.0;        // Volume average obtained during this sampling
         creal DVX = dvx / N; 
         creal DVY = dvy / N;
         creal DVZ = dvz / N;

         Real rhoFactor = 1.0;
         switch (densityModel) {
            case Uniform:
               rhoFactor = 1.0;
               break;
            case TestCase:
               rhoFactor = 1.0;
               if ((x >= 3.9e5 && x <= 6.1e5) && (y >= 3.9e5 && y <= 6.1e5)) {
                  rhoFactor = 1.5;
               }
               break;
            default:
               rhoFactor = 1.0;
               break;
         }
         
         // Sample the distribution using N*N*N points
         for (uint vi=0; vi<N; ++vi) {
            for (uint vj=0; vj<N; ++vj) {
               for (uint vk=0; vk<N; ++vk) {
                  creal VX = vx + 0.5*DVX + vi*DVX;
                  creal VY = vy + 0.5*DVY + vj*DVY;
                  creal VZ = vz + 0.5*DVZ + vk*DVZ;
                  avg += getDistribValue(VX,VY,VZ,DVX,DVY,DVZ,popID);
               }
            }
         }
         avg *= rhoFactor;
         
         // Compare the current and accumulated volume averages:
         Real eps = max(numeric_limits<creal>::min(),avg * static_cast<Real>(1e-6));
         Real avgAccum   = avgTotal / (avg + N3_sum);
         Real avgCurrent = avg / (N*N*N);
         if (fabs(avgCurrent-avgAccum)/(avgAccum+eps) < 0.01) ok = true;
         else if (avg < avgLimit) ok = true;
         else if (N > 10) {
            ok = true;
         }

         avgTotal += avg;
         N3_sum += N*N*N;
         ++N;
      } while (ok == false);

      return avgTotal / N3_sum;
   }

   void MultiPeak::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      setRandomCellSeed(cell,cellParams);

      if (this->lambda != 0.0) {
         cellParams[CellParams::PERBX] = this->dBx*cos(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
         cellParams[CellParams::PERBY] = this->dBy*sin(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
         cellParams[CellParams::PERBZ] = this->dBz*cos(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
      }

      cellParams[CellParams::PERBX] += this->magXPertAbsAmp * (0.5 - getRandomNumber(cell));
      cellParams[CellParams::PERBY] += this->magYPertAbsAmp * (0.5 - getRandomNumber(cell));
      cellParams[CellParams::PERBZ] += this->magZPertAbsAmp * (0.5 - getRandomNumber(cell));

      rhoRnd.clear();
      for (uint i=0; i<numberOfPopulations; ++i) {
         rhoRnd.push_back(rho[i] + rhoPertAbsAmp[i] * (0.5 - getRandomNumber(cell)));
      }
   }

   void MultiPeak::setActivePopulation(const uint popID) {
      this->popID = popID;
   }

   void MultiPeak::setCellBackgroundField(SpatialCell* cell) const {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
   
   std::vector<std::array<Real, 3> > MultiPeak::getV0(
                                                creal x,
                                                creal y,
                                                creal z,
                                                const uint popID
                                               ) const {
      vector<std::array<Real, 3> > centerPoints;
      for(uint i=0; i<this->numberOfPopulations; i++) {
         array<Real, 3> point {{this->Vx[i], this->Vy[i], this->Vz[i]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }
   
}// namespace projects
