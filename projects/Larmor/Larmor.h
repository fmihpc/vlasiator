/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef LARMOR_H
#define LARMOR_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"


namespace projects {
   class Larmor: public Project {
      public:
         Larmor();
         virtual ~Larmor();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         virtual void setCellBackgroundField(SpatialCell* cell);
      protected:
         Real getDistribValue(
         creal& x,creal& y, creal& z,
         creal& vx, creal& vy, creal& vz
         );
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz,const int& popID
         );
         
         Real BX0;
         Real BY0;
         Real BZ0;
         Real VX0;
         Real VY0;
         Real VZ0;
         Real DENSITY;
         Real TEMPERATURE;
         Real maxwCutoff; 
         uint nSpaceSamples;
         uint nVelocitySamples;
         Real SCA_X;
         Real SCA_Y;
         
   }; //Class Larmor
} // namespace projects
   

#endif
