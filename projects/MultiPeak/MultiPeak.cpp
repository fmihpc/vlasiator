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
      RP::add("MultiPeak.rho1", "Number density, first peak (m^-3)", 0.0);
      RP::add("MultiPeak.rho2", "Number density, second peak (m^-3)", 0.0);
      RP::add("MultiPeak.Tx1", "Temperature, first peak (K)", 0.0);
      RP::add("MultiPeak.Tx2", "Temperature, second peak (K)", 0.0);
      RP::add("MultiPeak.Ty1", "Temperature, first peak (K)", 0.0);
      RP::add("MultiPeak.Ty2", "Temperature, second peak (K)", 0.0);
      RP::add("MultiPeak.Tz1", "Temperature, first peak (K)", 0.0);
      RP::add("MultiPeak.Tz2", "Temperature, second peak (K)", 0.0);
      RP::add("MultiPeak.Vx1", "Bulk velocity x component, first peak (m/s)", 0.0);
      RP::add("MultiPeak.Vx2", "Bulk velocity x component, second peak (m/s)", 0.0);
      RP::add("MultiPeak.Vy1", "Bulk velocity y component, first peak (m/s)", 0.0);
      RP::add("MultiPeak.Vy2", "Bulk velocity y component, second peak (m/s)", 0.0);
      RP::add("MultiPeak.Vz1", "Bulk velocity z component, first peak (m/s)", 0.0);
      RP::add("MultiPeak.Vz2", "Bulk velocity z component, second peak (m/s)", 0.0);
      RP::add("MultiPeak.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("MultiPeak.By", "Magnetic field y component (T)", 0.0);
      RP::add("MultiPeak.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("MultiPeak.dBx", "Magnetic field x component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBy", "Magnetic field y component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.dBz", "Magnetic field z component cosine perturbation amplitude (T)", 0.0);
      RP::add("MultiPeak.magXPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along x (T)", 1.0e-9);
      RP::add("MultiPeak.magYPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along y (T)", 1.0e-9);
      RP::add("MultiPeak.magZPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along z (T)", 1.0e-9);
      RP::add("MultiPeak.rho1PertAbsAmp", "Absolute amplitude of the density perturbation, first peak", 0.1);
      RP::add("MultiPeak.rho2PertAbsAmp", "Absolute amplitude of the density perturbation, second peak", 0.1);
//       RP::add("MultiPeak.Vx1PertAbsAmp", "Absolute amplitude of the Vx perturbation, first peak", 1.0e6);
//       RP::add("MultiPeak.Vy1PertAbsAmp", "Absolute amplitude of the Vy perturbation, first peak", 1.0e6);
//       RP::add("MultiPeak.Vz1PertAbsAmp", "Absolute amplitude of the Vz perturbation, first peak", 1.0e6);
//       RP::add("MultiPeak.Vx2PertAbsAmp", "Absolute amplitude of the Vx perturbation, second peak", 1.0e6);
//       RP::add("MultiPeak.Vy2PertAbsAmp", "Absolute amplitude of the Vy perturbation, second peak", 1.0e6);
//       RP::add("MultiPeak.Vz2PertAbsAmp", "Absolute amplitude of the Vz perturbation, second peak", 1.0e6);
      RP::add("MultiPeak.lambda", "B cosine perturbation wavelength (m)", 0.0);
      RP::add("MultiPeak.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }

   void MultiPeak::getParameters(){
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("MultiPeak.rho1", this->rho[0]);
      RP::get("MultiPeak.rho2", this->rho[1]);
      RP::get("MultiPeak.Tx1", this->Tx[0]);
      RP::get("MultiPeak.Tx2", this->Tx[1]);
      RP::get("MultiPeak.Ty1", this->Ty[0]);
      RP::get("MultiPeak.Ty2", this->Ty[1]);
      RP::get("MultiPeak.Tz1", this->Tz[0]);
      RP::get("MultiPeak.Tz2", this->Tz[1]);
      RP::get("MultiPeak.Vx1", this->Vx[0]);
      RP::get("MultiPeak.Vx2", this->Vx[1]);
      RP::get("MultiPeak.Vy1", this->Vy[0]);
      RP::get("MultiPeak.Vy2", this->Vy[1]);
      RP::get("MultiPeak.Vz1", this->Vz[0]);
      RP::get("MultiPeak.Vz2", this->Vz[1]);
      RP::get("MultiPeak.Bx", this->Bx);
      RP::get("MultiPeak.By", this->By);
      RP::get("MultiPeak.Bz", this->Bz);
      RP::get("MultiPeak.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("MultiPeak.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("MultiPeak.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("MultiPeak.rho1PertAbsAmp", this->rhoPertAbsAmp[0]);
      RP::get("MultiPeak.rho2PertAbsAmp", this->rhoPertAbsAmp[1]);
//       RP::get("MultiPeak.Vx1PertAbsAmp", this->Vx1PertAbsAmp);
//       RP::get("MultiPeak.Vy1PertAbsAmp", this->Vy1PertAbsAmp);
//       RP::get("MultiPeak.Vz1PertAbsAmp", this->Vz1PertAbsAmp);
//       RP::get("MultiPeak.Vx2PertAbsAmp", this->Vx2PertAbsAmp);
//       RP::get("MultiPeak.Vy2PertAbsAmp", this->Vy2PertAbsAmp);
//       RP::get("MultiPeak.Vz2PertAbsAmp", this->Vz2PertAbsAmp);
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
      for(uint i=0; i<2; i++) {
         value += this->rhoRnd[i] * pow(mass / (2.0 * M_PI * kb ), 1.5) * 1.0 / sqrt(this->Tx[i]*this->Ty[i]*this->Tz[i]) *
      exp(- mass * (pow(vx - this->Vx[i], 2.0) / (2.0 * kb * this->Tx[i]) + pow(vy - this->Vy[i], 2.0) / (2.0 * kb * this->Ty[i]) + pow(vz - this->Vz[i], 2.0) / (2.0 * kb * this->Tz[i])));
      }
      return value;
   }

   Real MultiPeak::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
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
      if(this->lambda != 0.0) {
         cellParams[CellParams::PERBX] = this->dBx*cos(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
         cellParams[CellParams::PERBY] = this->dBy*sin(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
         cellParams[CellParams::PERBZ] = this->dBz*cos(2.0 * M_PI * cellParams[CellParams::XCRD] / this->lambda);
      }
      
      cellParams[CellParams::PERBX] += this->magXPertAbsAmp * (0.5 - getRandomNumber());
      cellParams[CellParams::PERBY] += this->magYPertAbsAmp * (0.5 - getRandomNumber());
      cellParams[CellParams::PERBZ] += this->magZPertAbsAmp * (0.5 - getRandomNumber());
      
      for(uint i=0; i<2; i++) {
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
      for(uint i=0; i<2; i++) {
         std::array<Real, 3> point {{this->Vx[i], this->Vy[i], this->Vz[i]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }
   
}// namespace projects
