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
      
      // initialize that vector to avoid segmentation faults
      for(uint i=0; i<this->numberOfPopulations; i++) {
         this->rhoRnd.push_back(0.0);
      }
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
// FIXME allow only 1 sample
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
   //#pragma omp parallel for collapse(6) reduction(+:avg)
      for (uint vi=0; vi<this->nVelocitySamples; ++vi)
      for (uint vj=0; vj<this->nVelocitySamples; ++vj)
         for (uint vk=0; vk<this->nVelocitySamples; ++vk)
            {
               avg += getDistribValue(vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
            }
            return avg / (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
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
      
      for(uint i=0; i<this->numberOfPopulations; i++) {
         this->rhoRnd[i] = this->rho[i] + this->rhoPertAbsAmp[i] * (0.5 - getRandomNumber());
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
