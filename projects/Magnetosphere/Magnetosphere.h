/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef MAGNETOSPHERE_H
#define MAGNETOSPHERE_H

#include "../../definitions.h"
#include "../projectIsotropicMaxwellian.h"

namespace projects {
   class Magnetosphere: public IsotropicMaxwellian {
      public:
         Magnetosphere();
         virtual ~Magnetosphere();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         virtual void setCellBackgroundField(SpatialCell* cell);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         
         virtual Real getV0(
            creal x,
            creal y,
            creal z,
            cuint component
         );
         
         Real T;
         Real tailRho;
         Real V0[3];
         Real ionosphereV0[3];
         Real constBgB[3];
         bool noDipoleInSW;
         Real rhoTransitionCenter;
         Real rhoTransitionWidth;
         Real ionosphereRho;
         Real ionosphereRadius;
         Real ionosphereTaperRadius;
         Real dipoleScalingFactor;
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class Magnetosphere
} // namespace projects

#endif

