/*
This file is part of Vlasiator.

Copyright 2011-2012,2015 Finnish Meteorological Institute

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
      void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;

   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz);
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz,
         const int& popID
      ) const;
      Real cellPosition;
      Real peakValue;
   }; // class test_trans
} // namespace projects

#endif

