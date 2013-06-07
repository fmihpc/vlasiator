/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef HARM1D_H
#define HARM1D_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class harm1D: public Project {
      public:
      harm1D();
      virtual ~harm1D();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      
      protected:
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
   }; // class harm1D
} // namespace projects

#endif
