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

#include "Template.h"

using namespace std;

namespace projects {
   Template::Template(): IsotropicMaxwellian() { }
   Template::~Template() { }
   
   void Template::addParameters() {
      typedef Readparameters RP;
      RP::add("Template.param", "This is my project's parameter. Default is 0.0", 0.0);
   }
   
   void Template::getParameters(){
      typedef Readparameters RP;
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!RP::get("Template.param", this->param)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
   }
   
   bool Template::initialize() {
      this->param += 1.0;
      return true;
   }

   Real Template::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      creal rho = 1.0;
      creal T = 1.0;
      return rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
      exp(- physicalconstants::MASS_PROTON * ((vx-this->getV0(x,y,z,0))*(vx-this->getV0(x,y,z,0)) + (vy-this->getV0(x,y,z,1))*(vy-this->getV0(x,y,z,1)) + (vz-this->getV0(x,y,z,2))*(vz-this->getV0(x,y,z,2))) / (2.0 * physicalconstants::K_B * T));
   }
   
   void Template::setCellBackgroundField(SpatialCell *cell){
      Dipole bgField;
      bgField.initialize(8e15,0.0); //set dipole moment
      if(cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
      } else {
         setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
      }
   }
   
   Real Template::getV0(
      creal x,
      creal y,
      creal z,
      cuint component
   ) {
      Real V0 = 0.0;
      switch(component) {
         case 0:
            V0 = 0.0;
            break;
         case 1:
            if(x < 0.0) V0 = 1.0;
            else V0 = 0.0;
            break;
         case 2:
            V0 = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ": Wrong component. This should not happen." << endl;
      }
      return V0;
   }
   
} // namespace projects

