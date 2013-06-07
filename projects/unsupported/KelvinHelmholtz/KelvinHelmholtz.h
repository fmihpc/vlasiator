/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef KELVINHELMHOLTZ_H
#define KELVINHELMHOLTZ_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class KelvinHelmholtz: public Project {
      public:
         KelvinHelmholtz();
         virtual ~KelvinHelmholtz();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual void setCellBackgroundField(SpatialCell* cell);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
      protected:
         Real getDistribValue(
            creal& x, creal& z,
            creal& vx, creal& vy, creal& vz);
         Real profile(creal top, creal bottom, creal x, creal z);
         
         enum {
            TOP,
            BOTTOM
         };
         Real rho[2];
         Real T[2];
         Real Vx[2];
         Real Vy[2];
         Real Vz[2];
         Real Bx[2];
         Real By[2];
         Real Bz[2];
         Real lambda;
         Real amp;
         Real offset;
         Real transitionWidth;
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class KelvinHelmholtz
} // namespace
#endif
