/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef ALFVEN_H
#define ALFVEN_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Alfven: public Project {
    public:
      Alfven();
      virtual ~Alfven();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      
    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,const int& popID
                          );
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,const int& popID
                                        );
      
      Real B0;
      Real Bx_guiding;
      Real By_guiding;
      Real Bz_guiding;
      Real DENSITY;
      Real ALPHA;
      Real WAVELENGTH;
      Real TEMPERATURE;
      Real A_VEL;
      Real A_MAG;
      uint nSpaceSamples;
      uint nVelocitySamples;
   } ; // class Alfven
} // namespace projects

#endif
