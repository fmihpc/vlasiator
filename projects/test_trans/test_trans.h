/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef TEST_TRANS_H
#define TEST_TRANS_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class test_trans: public Project {
   public:
      test_trans();
      virtual ~test_trans();

            virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      void setCellBackgroundField(SpatialCell* cell);

   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz);
      virtual void calcCellParameters(Real* cellParams,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz
      );
      Real cellPosition;
   }; // class test_trans
} // namespace projects

#endif

