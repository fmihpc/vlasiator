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


#include "MultiPeak.h"

using namespace std;
using namespace spatial_cell;


Real projects::MultiPeak::rhoRnd;

namespace projects {
   MultiPeak::MultiPeak(): TriAxisSearch() { }
   
   MultiPeak::~MultiPeak() { }

   bool MultiPeak::initialize(void) {
      return Project::initialize();
   }

   void MultiPeak::addParameters(){
      typedef Readparameters RP;

      RP::add("MultiPeak.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("MultiPeak.By", "Magnetic field y component (T)", 0.0);
      RP::add("MultiPeak.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("MultiPeak.dBx", "Magnetic field x component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBy", "Magnetic field y component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBz", "Magnetic field z component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.magXPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along x (T)", 1.0e-9);
      RP::add("MultiPeak.magYPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along y (T)", 1.0e-9);
      RP::add("MultiPeak.magZPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along z (T)", 1.0e-9);
      RP::add("MultiPeak.lambda", "B cosine perturbation wavelength (m)", 1.0);
      RP::add("MultiPeak.nVelocitySamples", "Number of sampling points per velocity dimension", 2);
      RP::add("MultiPeak.densityModel","Which spatial density model is used?",string("uniform"));

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         RP::add(pop+"_MultiPeak.n", "Number of peaks to create", 0);
         RP::addComposing(pop+"_MultiPeak.rho", "Number density (m^-3)");
         RP::addComposing(pop+"_MultiPeak.Tx", "Temperature (K)");
         RP::addComposing(pop+"_MultiPeak.Ty", "Temperature");
         RP::addComposing(pop+"_MultiPeak.Tz", "Temperature");
         RP::addComposing(pop+"_MultiPeak.Vx", "Bulk velocity x component (m/s)");
         RP::addComposing(pop+"_MultiPeak.Vy", "Bulk velocity y component (m/s)");
         RP::addComposing(pop+"_MultiPeak.Vz", "Bulk velocity z component (m/s)");
         RP::addComposing(pop+"_MultiPeak.rhoPertAbsAmp", "Absolute amplitude of the density perturbation");
      }
   }

   void MultiPeak::getParameters(){

      typedef Readparameters RP;
      Project::getParameters();
      RP::get("MultiPeak.Bx", this->Bx);
      RP::get("MultiPeak.By", this->By);
      RP::get("MultiPeak.Bz", this->Bz);
      RP::get("MultiPeak.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("MultiPeak.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("MultiPeak.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("MultiPeak.dBx", this->dBx);
      RP::get("MultiPeak.dBy", this->dBy);
      RP::get("MultiPeak.dBz", this->dBz);
      RP::get("MultiPeak.lambda", this->lambda);
      RP::get("MultiPeak.nVelocitySamples", this->nVelocitySamples);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         MultiPeakSpeciesParameters sP;
         RP::get(pop + "_MultiPeak.n", sP.numberOfPeaks);
         RP::get(pop + "_MultiPeak.rho",sP.rho);
         RP::get(pop + "_MultiPeak.Tx", sP.Tx);
         RP::get(pop + "_MultiPeak.Ty", sP.Ty);
         RP::get(pop + "_MultiPeak.Tz", sP.Tz);
         RP::get(pop + "_MultiPeak.Vx", sP.Vx);
         RP::get(pop + "_MultiPeak.Vy", sP.Vy);
         RP::get(pop + "_MultiPeak.Vz", sP.Vz);

         RP::get(pop + "_MultiPeak.rhoPertAbsAmp", sP.rhoPertAbsAmp);
         if(!sP.isConsistent()) {
            cerr << "You should define all parameters (MultiPeak.rho, MultiPeak.Tx, MultiPeak.Ty, MultiPeak.Tz, MultiPeak.Vx, MultiPeak.Vy, MultiPeak.Vz, MultiPeak.rhoPertAbsAmp) for all " << sP.numberOfPeaks << " peaks of population " << pop << "." << endl;
            abort();
         }

         speciesParams.push_back(sP);
      }
      
      string densModelString;
      RP::get("MultiPeak.densityModel",densModelString);
      
      if (densModelString == "uniform") densityModel = Uniform;
      else if (densModelString == "testcase") densityModel = TestCase;
   }

   Real MultiPeak::getDistribValue(creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      const MultiPeakSpeciesParameters& sP = speciesParams[popID];
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;

      Real value = 0.0;

      for (uint i=0; i<sP.numberOfPeaks; ++i) {
         value += (sP.rho[i] + sP.rhoPertAbsAmp[i] * rhoRnd)
               * pow(mass / (2.0 * M_PI * kb ), 1.5) * 1.0
               / sqrt(sP.Tx[i]*sP.Ty[i]*sP.Tz[i]) 
               * exp(- mass * (pow(vx - sP.Vx[i], 2.0) / (2.0 * kb * sP.Tx[i]) 
                             + pow(vy - sP.Vy[i], 2.0) / (2.0 * kb * sP.Ty[i]) 
                             + pow(vz - sP.Vz[i], 2.0) / (2.0 * kb * sP.Tz[i])));
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

      const MultiPeakSpeciesParameters& sP = speciesParams[popID];
                                            
      #warning TODO: Replace getObjectWrapper().particleSpecies[popID].sparseMinValue with SpatialCell::getVelocityBlockMinValue(popID)
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
      setRandomCellSeed(cell);
      rhoRnd = 0.5 - getRandomNumber();
   }

   void MultiPeak::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField, BgBGrid);
      
      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();
         
#pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
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
   
   std::vector<std::array<Real, 3> > MultiPeak::getV0(
                                                creal x,
                                                creal y,
                                                creal z,
                                                const uint popID
                                               ) const {
      const MultiPeakSpeciesParameters& sP = speciesParams[popID];
      vector<std::array<Real, 3> > centerPoints;
      for(uint i=0; i<sP.numberOfPeaks; i++) {
         array<Real, 3> point {{sP.Vx[i], sP.Vy[i], sP.Vz[i]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }
   
}// namespace projects
