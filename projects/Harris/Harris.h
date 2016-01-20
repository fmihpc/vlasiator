/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute



*/

#ifndef HARRIS_H
#define HARRIS_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Harris: public Project {
    public:
      Harris();
      virtual ~Harris();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,const int& popID
                                        );
      
    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz
                          );
      
      Real SCA_LAMBDA;
      Real BX0, BY0, BZ0;
      Real TEMPERATURE;
      Real DENSITY;
   }; // class Harris
} // namespace Harris

#endif
