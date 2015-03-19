/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef VERIFICATIONLARMOR_H
#define VERIFICATIONLARMOR_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class verificationLarmor: public Project {
   public:
      verificationLarmor();
      virtual ~verificationLarmor();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell);
   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz);
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
      Real X0;
      Real Y0;
      Real Z0;
      Real DENSITY;

   }; // class verificationLarmor
} // namespace projects

#endif

