/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

*/


#ifndef VELOCITYBOX_H
#define VELOCITYBOX_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class VelocityBox: public Project {
      public:
         VelocityBox();
         virtual ~VelocityBox();
         
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
            creal& dvx, creal& dvy, creal& dvz
            );
         
         Real rho;
         Real Vx[2];
         Real Vy[2];
         Real Vz[2];
         Real Bx;
         Real By;
         Real Bz;
   }; // class VelocityBox
} //  namespace projects



#endif

