/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "MultiPeak.h"

using namespace std;

vector<Real> projects::MultiPeak::rhoRnd;

namespace projects {
   MultiPeak::MultiPeak(): TriAxisSearch() { }
   MultiPeak::~MultiPeak() { }


   bool MultiPeak::initialize(void) {return true;}

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
   }

   void MultiPeak::getParameters(){
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("MultiPeak.n", this->numberOfPopulations);
      if(this->numberOfPopulations < 1) {
         std::cerr << "You should set MultiPeak.n to more than 0 populations." << std::endl;
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
         std::cerr << "You should define all parameters (MultiPeak.rho, MultiPeak.Tx, MultiPeak.Ty, MultiPeak.Tz, MultiPeak.Vx, MultiPeak.Vy, MultiPeak.Vz, MultiPeak.rhoPertAbsAmp) for each of the populations." << std::endl;
         abort();
      }
      if(this->numberOfPopulations > this->rho.size()) {
         std::cerr << "You are requesting more populations than are currently defined. Change MultiPeak.n or define more populations." << std::endl;
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
      
   }

   Real MultiPeak::getDistribValue(creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;

      Real value = 0.0;
      for(uint i=0; i<this->numberOfPopulations; i++) {
         value += this->rhoRnd[i] * pow(mass / (2.0 * M_PI * kb ), 1.5) * 1.0 / sqrt(this->Tx[i]*this->Ty[i]*this->Tz[i]) *
      exp(- mass * (pow(vx - this->Vx[i], 2.0) / (2.0 * kb * this->Tx[i]) + pow(vy - this->Vy[i], 2.0) / (2.0 * kb * this->Ty[i]) + pow(vz - this->Vz[i], 2.0) / (2.0 * kb * this->Tz[i])));
      }
      return value;
   }

   Real MultiPeak::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      // Iterative sampling of the distribution function. Keep track of the 
      // accumulated volume average over the iterations. When the next 
      // iteration improves the average by less than 1%, return the value.
      Real avgTotal = 0.0;
      bool ok = false;
      int N = nVelocitySamples; // Start by using nVelocitySamples
      int N3_sum = 0;           // Sum of sampling points used so far
      do {
         Real avg = 0.0;        // Volume average obtained during this sampling
         creal DVX = dvx / N; 
         creal DVY = dvy / N;
         creal DVZ = dvz / N;

         // Sample the distribution using N*N*N points
         for (uint vi=0; vi<N; ++vi) {
            for (uint vj=0; vj<N; ++vj) {
               for (uint vk=0; vk<N; ++vk) {
                  creal VX = vx + 0.5*DVX + vi*DVX;
                  creal VY = vy + 0.5*DVY + vj*DVX;
                  creal VZ = vz + 0.5*DVZ + vk*DVX;
                  avg += getDistribValue(VX,VY,VZ,DVX,DVY,DVZ);
               }
            }
         }
         
         // Compare the current and accumulated volume averages:
         //ok = true;
         Real eps = max(numeric_limits<creal>::min(),avg * static_cast<Real>(1e-6));
         Real avgAccum   = avgTotal / (avg + N3_sum);
         Real avgCurrent = avg / (N*N*N);
         if (fabs(avgCurrent-avgAccum)/(avgAccum+eps) < 0.01) ok = true;
         else if (avg < Parameters::sparseMinValue*0.01) ok = true;
         else if (N > 10) {
            ok = true;
         }

         avgTotal += avg;
         N3_sum += N*N*N;
         ++N;
      } while (ok == false);

      return avgTotal / N3_sum;
   }

   void MultiPeak::calcCellParameters(Real* cellParams,creal& t) {
      setRandomCellSeed(cellParams);

      if (this->lambda != 0.0) {
         cellParams[CellParams::PERBX] = this->dBx*cos(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
         cellParams[CellParams::PERBY] = this->dBy*sin(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
         cellParams[CellParams::PERBZ] = this->dBz*cos(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
      }
      
      cellParams[CellParams::PERBX] += this->magXPertAbsAmp * (0.5 - getRandomNumber());
      cellParams[CellParams::PERBY] += this->magYPertAbsAmp * (0.5 - getRandomNumber());
      cellParams[CellParams::PERBZ] += this->magZPertAbsAmp * (0.5 - getRandomNumber());


      rhoRnd.clear();
      for(uint i=0; i<this->numberOfPopulations; i++) {
         this->rhoRnd.push_back(this->rho[i] + this->rhoPertAbsAmp[i] * (0.5 - getRandomNumber()));
      }
   }

   void MultiPeak::setCellBackgroundField(SpatialCell* cell) {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
   
   vector<std::array<Real, 3>> MultiPeak::getV0(
                                                creal x,
                                                creal y,
                                                creal z
                                               ) {
      vector<std::array<Real, 3>> centerPoints;
      for(uint i=0; i<this->numberOfPopulations; i++) {
         std::array<Real, 3> point {{this->Vx[i], this->Vy[i], this->Vz[i]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }
   
}// namespace projects
