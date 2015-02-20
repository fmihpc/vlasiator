/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef TEST_FP_H
#define TEST_FP_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class test_fp: public TriAxisSearch {
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
         creal& dvx, creal& dvy, creal& dvz,const int& popID
      );
      
      virtual vector<std::array<Real, 3>> getV0(
         creal x,
         creal y,
         creal z
      );
      
      virtual vector<std::array<Real, 3>> getV0(
         creal x,
         creal y,
         creal z,
         creal dx,
         creal dy,
         creal dz
      );
      
      Real V0;
      Real B0;
      Real DENSITY;
      Real TEMPERATURE;
      Real ALPHA;
      int CASE;
      bool shear;
   }; // class test_fp
}// namespace projects

#endif
