/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef FLOWTHROUGH_H
#define FLOWTHROUGH_H

#include "../../definitions.h"
#include "../projectIsotropicMaxwellian.h"


namespace projects {
   class Flowthrough: public IsotropicMaxwellian {
      public:
         Flowthrough();
         virtual ~Flowthrough();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         
         
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         virtual Real getV0(
            creal x,
            creal y,
            creal z,
            cuint component
         );
         
         Real rho;
         Real T;
         Real V0[3];
         Real Bx;
         Real By;
         Real Bz;
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class Flowthrough
} // namespace projects

#endif
