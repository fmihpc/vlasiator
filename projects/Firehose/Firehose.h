/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef FIREHOSE_H
#define FIREHOSE_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Firehose: public Project {
      public:
         Firehose();
         virtual ~Firehose();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
            
      protected:
         Real getDistribValue(
            creal& x,creal& y,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         Real profile(creal top, creal bottom, creal x);
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         
         Real rho[2];
         Real Tx[2];
         Real Ty[2];
         Real Tz[2];
         Real Vx[2];
         Real Vy[2];
         Real Vz[2];
         Real Bx;
         Real By;
         Real Bz;   
         Real lambda;
         Real amp;
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class Firehose
} // namespace projects

#endif
