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
   Harris::Harris(): Project() { }
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
   }
   
   void Harris::getParameters(){
      typedef Readparameters RP;
      RP::get("Harris.Scale_size", this->SCA_LAMBDA);
      RP::get("Harris.BX0", this->BX0);
      RP::get("Harris.BY0", this->BY0);
      RP::get("Harris.BZ0", this->BZ0);
      RP::get("Harris.Temperature", this->TEMPERATURE);
      RP::get("Harris.rho", this->DENSITY);
   }
   
   Real Harris::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz
   ) {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      creal Vy0 = 0.0;
      
      return this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) * (
         5.0 / pow(cosh((x + 0.5 * dx) / (this->SCA_LAMBDA)), 2.0) * exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy - Vy0, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * kb * this->TEMPERATURE))
         +
         exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * kb * this->TEMPERATURE)));
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
} // namespace projects
