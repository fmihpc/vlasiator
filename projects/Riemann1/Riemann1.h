/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef RIEMANN_H
#define RIEMANN_H

#include "../../definitions.h"
#include "../project.h"


namespace projects {
   class Riemann1: public Project {
      public:
         Riemann1();
         virtual ~Riemann1();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
      
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const int& popID
         );

      
         enum {
            LEFT,
            RIGHT
         };
         Real rho[2];
         Real T[2];
         Real Vx[2];
         Real Vy[2];
         Real Vz[2];
         Real Bx[2];
         Real By[2];
         Real Bz[2];
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class Riemann1
}  // namespace projects

#endif
