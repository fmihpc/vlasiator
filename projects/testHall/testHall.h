/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef TESTHALL_H
#define TESTHALL_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class TestHall: public Project {
      public:
         TestHall();
         virtual ~TestHall();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
//          virtual void setCellBackgroundField(SpatialCell* cell);
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const int& popID
         );
         
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         
         bool noDipoleInSW;
         Real constBgB[3];
         Real dipoleScalingFactor;
         Real dipoleTilt;
         Real BX0;
         Real BY0;
         Real BZ0;
         Real VX0;
         Real VY0;
         Real VZ0;
         Real TEMPERATURE;
         Real DENSITY;
   }; // class TestHall
} // namespace TestHall

#endif
