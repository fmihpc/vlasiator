
/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef DISPERSION_H
#define DISPERSION_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Dispersion: public Project {
    public:
      Dispersion();
      virtual ~Dispersion();
      
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
                                         creal& dvx, creal& dvy, creal& dvz
                                        );

      Real B0;
      Real VX0;
      Real VY0;
      Real VZ0;
      Real angleXY;
      Real angleXZ;
      Real DENSITY;
      Real TEMPERATURE;
      Real magXPertAbsAmp;
      Real magYPertAbsAmp;
      Real magZPertAbsAmp;
      Real densityPertRelAmp;
      Real velocityPertAbsAmp;
      Real maxwCutoff;
      uint seed;
      uint nSpaceSamples;
      uint nVelocitySamples;
      
      char rngStateBuffer[256];
      random_data rngDataBuffer;
      
      static Real rndRho, rndVel[3];
      #pragma omp threadprivate(rndRho,rndVel)
   } ; // class Dispersion
} // namespace projects

#endif
