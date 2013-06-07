/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License

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
   Harris::Harris(): Project() { }
   Harris::~Harris() { }
   
   bool Harris::initialize(void) {return true;}
   
   void Harris::addParameters(){
      typedef Readparameters RP;
      RP::add("Harris.Scale_size", "Harris sheet scale size (m)", 150000.0);
      RP::add("Harris.B0", "Magnetic field at infinity (T)", 8.33061003094e-8);
      RP::add("Harris.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Harris.rho", "Number density at infinity (m^-3)", 1.0e7);
   }
   
   void Harris::getParameters(){
      typedef Readparameters RP;
      RP::get("Harris.Scale_size", this->SCA_LAMBDA);
      RP::get("Harris.B0", this->B0);
      RP::get("Harris.Temperature", this->TEMPERATURE);
      RP::get("Harris.rho", this->DENSITY);
   }
   
   Real Harris::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz
   ) {
      creal mass = 1.67262171e-27; // m_p in kg
      creal k = 1.3806505e-23; // Boltzmann
      creal mu0 = 1.25663706144e-6; // mu_0
      creal q = 1.60217653e-19; // q_i
      
      creal Vy0 = 0.0;
      
      return this->DENSITY * pow(mass / (2.0 * M_PI * k * this->TEMPERATURE), 1.5) * (
         5.0 / pow(cosh((x + 0.5 * dx) / (this->SCA_LAMBDA)), 2.0) * exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy - Vy0, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * this->TEMPERATURE))
         +
         exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * this->TEMPERATURE)));
   }
   
   void Harris::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
      cellParams[CellParams::BGBX   ] = 0.0;
      cellParams[CellParams::BGBY   ] = 0.0;
      cellParams[CellParams::BGBZ   ] = this->B0 * tanh((x + 0.5 * dx) / this->SCA_LAMBDA);
   }
} // namespace projects
