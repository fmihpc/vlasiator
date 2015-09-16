/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef SHOCK_H
#define SHOCK_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Shock: public Project {
      public:
         Shock();
         virtual ~Shock();
      
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
      
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz
         ); 
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const int& popID
         );


         Real BX0;
         Real BY0;
         Real BZ0;
         Real EX0;
         Real VX0;
         Real VY0;
         Real VZ0;
         Real DENSITY;
         Real TEMPERATURE;
         Real magPertAmp;
         Real densityPertAmp;
         Real velocityPertAmp;
         Real maxwCutoff;
         uint nSpaceSamples;
         uint nVelocitySamples;
         Real SCA_X;
         Real SCA_Y;
         Real Sharp_Y;
   } ; //class Shock
} // namespace projects
#endif
