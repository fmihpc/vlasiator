/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef TEST_FP_H
#define TEST_FP_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class test_fp: public Project {
   public:
      test_fp();
      virtual ~test_fp();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      
      
   protected:
      Real sign(creal value);
      Real getDistribValue(creal& vx, creal& vy, creal& vz);
      virtual void calcCellParameters(Real* cellParams,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz
      );

      Real B0;
      Real DENSITY;
      Real TEMPERATURE;
      Real ALPHA;
      int CASE;
      bool shear;
   }; // class test_fp
}// namespace projects

#endif
