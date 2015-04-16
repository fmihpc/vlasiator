/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"

#include "Harris.h"

using namespace std;

namespace projects {
   Harris::Harris(): TriAxisSearch() { }
   Harris::~Harris() { }
   
   bool Harris::initialize(void) {return true;}
   
   void Harris::addParameters(){
      typedef Readparameters RP;
      RP::add("Harris.Scale_size", "Harris sheet scale size (m)", 150000.0);
      RP::add("Harris.BX0", "Magnetic field at infinity (T)", 8.33061003094e-8);
      RP::add("Harris.BY0", "Magnetic field at infinity (T)", 8.33061003094e-8);
      RP::add("Harris.BZ0", "Magnetic field at infinity (T)", 8.33061003094e-8);
      RP::add("Harris.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Harris.rho", "Number density at infinity (m^-3)", 1.0e7);
      RP::add("Harris.nSpaceSamples", "Number of sampling points per spatial dimension.", 2);
      RP::add("Harris.nVelocitySamples", "Number of sampling points per velocity dimension.", 2);
   }
   
   void Harris::getParameters(){
      typedef Readparameters RP;
      RP::get("Harris.Scale_size", this->SCA_LAMBDA);
      RP::get("Harris.BX0", this->BX0);
      RP::get("Harris.BY0", this->BY0);
      RP::get("Harris.BZ0", this->BZ0);
      RP::get("Harris.Temperature", this->TEMPERATURE);
      RP::get("Harris.rho", this->DENSITY);
      RP::get("Harris.nSpaceSamples", this->nSpaceSamples);
      RP::get("Harris.nVelocitySamples", this->nVelocitySamples);
   }
   
   Real Harris::getDistribValue(creal& x,creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      return this->DENSITY * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->TEMPERATURE), 1.5) * (
         5.0 / pow(cosh(x / (this->SCA_LAMBDA)), 2.0) * exp(- physicalconstants::MASS_PROTON * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * physicalconstants::K_B * this->TEMPERATURE))
         +
         exp(- physicalconstants::MASS_PROTON * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * physicalconstants::K_B * this->TEMPERATURE)));
   }
   
   Real Harris::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz
   ) {
      if((this->nSpaceSamples > 1) && (this->nVelocitySamples > 1)) {
         creal d_x = dx / (this->nSpaceSamples-1);
         creal d_y = dy / (this->nSpaceSamples-1);
         creal d_z = dz / (this->nSpaceSamples-1);
         creal d_vx = dvx / (this->nVelocitySamples-1);
         creal d_vy = dvy / (this->nVelocitySamples-1);
         creal d_vz = dvz / (this->nVelocitySamples-1);
         
         Real avg = 0.0;
         // #pragma omp parallel for collapse(6) reduction(+:avg)
         // WARNING No threading here if calling functions are already threaded
         for (uint i=0; i<this->nSpaceSamples; ++i)
            for (uint j=0; j<this->nSpaceSamples; ++j)
               for (uint k=0; k<this->nSpaceSamples; ++k)
                  for (uint vi=0; vi<this->nVelocitySamples; ++vi)
                     for (uint vj=0; vj<this->nVelocitySamples; ++vj)
                        for (uint vk=0; vk<this->nVelocitySamples; ++vk) {
                           avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
                        }
                        return avg /
                        (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) /
                        (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz);
      }
      
      
      
   }
   
   void Harris::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = this->BX0 * tanh((y + 0.5 * dy) / this->SCA_LAMBDA);
      cellParams[CellParams::PERBY   ] = this->BY0 * tanh((z + 0.5 * dz) / this->SCA_LAMBDA);
      cellParams[CellParams::PERBZ   ] = this->BZ0 * tanh((x + 0.5 * dx) / this->SCA_LAMBDA);
      cellParams[CellParams::BGBX   ] = 0.0;
      cellParams[CellParams::BGBY   ] = 0.0;
      cellParams[CellParams::BGBZ   ] = 0.0;
   }
   
   vector<std::array<Real, 3>> Harris::getV0(
      creal x,
      creal y,
      creal z
   ) {
      vector<std::array<Real, 3>> V0;
      std::array<Real, 3> v = {{0.0, 0.0, 0.0 }};
      V0.push_back(v);
      return V0;
   }
} // namespace projects
