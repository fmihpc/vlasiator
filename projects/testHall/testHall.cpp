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
#include "../../backgroundfield/dipole.hpp"

#include "testHall.h"

using namespace std;

namespace projects {
   TestHall::TestHall(): Project() { }
   TestHall::~TestHall() { }
   
   bool TestHall::initialize(void) {
      bool success = Project::initialize();
      this->constBgB[0] = 0.0;
      this->constBgB[1] = 0.0;
      this->constBgB[2] = 0.0;
      this->dipoleScalingFactor = 1.0;
      this->dipoleTilt = 0.0;
      this->noDipoleInSW = 0;
      return success;
   }
   
   void TestHall::addParameters(){
      typedef Readparameters RP;
      RP::add("TestHall.BX0", "Magnetic field x (T)", 1.0e-9);
      RP::add("TestHall.BY0", "Magnetic field y (T)", 1.0e-9);
      RP::add("TestHall.BZ0", "Magnetic field z (T)", 1.0e-9);
      RP::add("TestHall.VX0", "velocity x (m/s)", -1.0e3);
      RP::add("TestHall.VY0", "velocity y (m/s)", 1.0e3);
      RP::add("TestHall.VZ0", "velocity z (m/s)", 1.0e3);
      RP::add("TestHall.Temperature", "Temperature (K)", 1.0e6);
      RP::add("TestHall.rho", "Number density (m^-3)", 1.0e6);
   }
   
   void TestHall::getParameters(){
      Project::getParameters();
      typedef Readparameters RP;
      RP::get("TestHall.BX0", this->BX0);
      RP::get("TestHall.BY0", this->BY0);
      RP::get("TestHall.BZ0", this->BZ0);
      RP::get("TestHall.VX0", this->VX0);
      RP::get("TestHall.VY0", this->VY0);
      RP::get("TestHall.VZ0", this->VZ0);
      RP::get("TestHall.Temperature", this->TEMPERATURE);
      RP::get("TestHall.rho", this->DENSITY);
   }
   
   Real TestHall::calcPhaseSpaceDensity(
      creal& x,creal& y,creal& z,
      creal& dx,creal& dy,creal& dz,
      creal& vx,creal& vy,creal& vz,
      creal& dvx,creal& dvy,creal& dvz,const int& popID
   ) {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      return this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) * (
         exp(- mass * (pow(vx + 0.5 * dvx - this->VX0, 2.0) + pow(vy + 0.5 * dvy - this->VY0, 2.0) + pow(vz + 0.5 * dvz - this->VZ0, 2.0)) / (2.0 * kb * this->TEMPERATURE)));
   }
   
//    void TestHall::setCellBackgroundField(SpatialCell *cell){
//       Dipole bgField;
//       bgField.initialize(8e15 *this->dipoleScalingFactor,this->dipoleTilt); //set dipole moment
//       if(cell->sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN && this->noDipoleInSW) {
//          setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
//       } else {
//          setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
//       }
//       
//       cell->parameters[CellParams::EX   ] = 0.0;
//       cell->parameters[CellParams::EY   ] = 0.0;
//       cell->parameters[CellParams::EZ   ] = 0.0;
//       cell->parameters[CellParams::PERBX  ] = cell->parameters[CellParams::BGBX];
//       cell->parameters[CellParams::BGBX   ] = 0.0;
//       cell->parameters[CellParams::BGBXVOL] = 0.0;
//       cell->parameters[CellParams::PERBY  ] = cell->parameters[CellParams::BGBY];
//       cell->parameters[CellParams::BGBY   ] = 0.0;
//       cell->parameters[CellParams::BGBYVOL] = 0.0;
//       cell->parameters[CellParams::PERBZ  ] = cell->parameters[CellParams::BGBZ];
//       cell->parameters[CellParams::BGBZ   ] = 0.0;
//       cell->parameters[CellParams::BGBZVOL] = 0.0;
//       
//       cell->derivatives[fieldsolver::dBGBydx]=0.0;
//       cell->derivatives[fieldsolver::dBGBzdx]=0.0;
//       cell->derivatives[fieldsolver::dBGBxdy]=0.0;
//       cell->derivatives[fieldsolver::dBGBxdz]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBYVOLdx]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBZVOLdx]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBXVOLdy]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBXVOLdz]=0.0;
//       cell->derivatives[fieldsolver::dBGBxdy]=0.0;
//       cell->derivatives[fieldsolver::dBGBzdy]=0.0;
//       cell->derivatives[fieldsolver::dBGBydx]=0.0;
//       cell->derivatives[fieldsolver::dBGBydz]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBXVOLdy]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBZVOLdy]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBYVOLdx]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBYVOLdz]=0.0;
//       cell->derivatives[fieldsolver::dBGBxdy]=0.0;
//       cell->derivatives[fieldsolver::dBGBxdz]=0.0;
//       cell->derivatives[fieldsolver::dBGBydx]=0.0;
//       cell->derivatives[fieldsolver::dBGBydz]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBXVOLdy]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBXVOLdz]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBYVOLdx]=0.0;
//       cell->derivativesBVOL[bvolderivatives::dBGBYVOLdz]=0.0;
//       
//       for(uint component=0; component<3; component++) {
//          if(this->constBgB[component] != 0.0) {
//             cell->parameters[CellParams::BGBX+component] += this->constBgB[component];
//             cell->parameters[CellParams::BGBXVOL+component] += this->constBgB[component];
//          }
//       }
//    }
//    
   void TestHall::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      creal Dx = 0.5*dx;
      creal Dy = 0.5*dy;
      creal Dz = 0.5*dz;
      
//       creal r = sqrt((x+Dx)*(x+Dx) + (y+Dy)*(y+Dy));
//       creal theta = atan2(y+Dy, x+Dx);
      
//       creal I = 1.0e6; // current
//       creal B = physicalconstants::MU_0 * I / (2.0 * 3.1415927);
      
//       cellParams[CellParams::PERBX] = this->BX0;
//       cellParams[CellParams::PERBY] = this->BY0;
//       cellParams[CellParams::PERBZ] = this->BZ0;
      
//       cellParams[CellParams::PERBX] = this->BX0 * y;
//       cellParams[CellParams::PERBY] = this->BY0 * z;
//       cellParams[CellParams::PERBZ] = this->BZ0 * x;
      
      cellParams[CellParams::PERBX] = this->BX0 * cos(2.0*M_PI * 1.0 * x / (P::xmax - P::xmin)) * cos(2.0*M_PI * 1.0 * y / (P::ymax - P::ymin)) * cos(2.0*M_PI * 1.0 * z / (P::zmax - P::zmin));
      cellParams[CellParams::PERBY] = this->BY0 * cos(2.0*M_PI * 1.0 * x / (P::xmax - P::xmin)) * cos(2.0*M_PI * 1.0 * y / (P::ymax - P::ymin)) * cos(2.0*M_PI * 1.0 * z / (P::zmax - P::zmin));
      cellParams[CellParams::PERBZ] = this->BZ0 * cos(2.0*M_PI * 1.0 * x / (P::xmax - P::xmin)) * cos(2.0*M_PI * 1.0 * y / (P::ymax - P::ymin)) * cos(2.0*M_PI * 1.0 * z / (P::zmax - P::zmin));
      
//       cellParams[CellParams::PERBX] = -1.0*(y+Dy) / ((x+Dx)*(x+Dx) + (y+Dy)*(y+Dy));
//       cellParams[CellParams::PERBY] = (x+Dx) / ((x+Dx)*(x+Dx) + (y+Dy)*(y+Dy));
//       cellParams[CellParams::PERBZ] = 0.0;
      
//       cellParams[CellParams::PERBX] = this->BX0 * tanh((y + 0.5 * dy) / (15.0 * dy));
//       cellParams[CellParams::PERBY] = this->BY0 * tanh((z + 0.5 * dz) / (15.0 * dz));
//       cellParams[CellParams::PERBZ] = this->BZ0 * tanh((x + 0.5 * dx) / (15.0 * dx));
      
//       cellParams[CellParams::PERBX   ] = this->BX0 * (x+0.5*Dx + y+0.5*Dy + (z+0.5*Dz));
//       cellParams[CellParams::PERBY   ] = this->BY0 * ((x+0.5*Dx)*(x+0.5*Dx) + (y+0.5*Dy)*(y+0.5*Dy) + (z+0.5*Dz)*(z+0.5*Dz));
//       cellParams[CellParams::PERBX   ] = this->BX0 * ((x+0.5*Dx)*(x+0.5*Dx)*(x+0.5*Dx)/ pow(Parameters::xmax - Parameters::xmin, 3.0) + (y+0.5*Dy)*(y+0.5*Dy)*(y+0.5*Dy)/ pow(Parameters::ymax - Parameters::ymin, 3.0) + (z+0.5*Dz)*(z+0.5*Dz)*(z+0.5*Dz)/ pow(Parameters::zmax - Parameters::zmin, 3.0))   ;
//       cellParams[CellParams::PERBY   ] = this->BY0 * ((x+0.5*Dx)*(x+0.5*Dx)*(x+0.5*Dx)/ pow(Parameters::xmax - Parameters::xmin, 3.0) + (y+0.5*Dy)*(y+0.5*Dy)*(y+0.5*Dy)/ pow(Parameters::ymax - Parameters::ymin, 3.0) + (z+0.5*Dz)*(z+0.5*Dz)*(z+0.5*Dz)/ pow(Parameters::zmax - Parameters::zmin, 3.0));
//       cellParams[CellParams::PERBZ   ] = this->BZ0 * ((x+0.5*Dx)*(x+0.5*Dx)*(x+0.5*Dx)/ pow(Parameters::xmax - Parameters::xmin, 3.0) + (y+0.5*Dy)*(y+0.5*Dy)*(y+0.5*Dy)/ pow(Parameters::ymax - Parameters::ymin, 3.0) + (z+0.5*Dz)*(z+0.5*Dz)*(z+0.5*Dz)/ pow(Parameters::zmax - Parameters::zmin, 3.0));
      
//       cellParams[CellParams::PERBX   ] = this->BX0 * (x+0.5*Dx)*(y+0.5*Dy)*(z+0.5*Dz);
//       cellParams[CellParams::PERBY   ] = this->BY0 * (x+0.5*Dx)*(y+0.5*Dy)*(z+0.5*Dz)*(x+0.5*Dx)*(y+0.5*Dy)*(z+0.5*Dz);
//       cellParams[CellParams::PERBZ   ] = this->BZ0 * (x+0.5*Dx)*(y+0.5*Dy)*(z+0.5*Dz)*(x+0.5*Dx)*(y+0.5*Dy)*(z+0.5*Dz)*(x+0.5*Dx)*(y+0.5*Dy)*(z+0.5*Dz);
   }
} // namespace projects
