/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef FLUCTUATIONS_H
#define FLUCTUATIONS_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class Fluctuations: public TriAxisSearch {
   public:
      Fluctuations();
      virtual ~Fluctuations();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
      virtual std::vector<std::array<Real, 3> > getV0(
         creal x,
         creal y,
         creal z
      ) const;
   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz,const int& popID
      ) const;
      
      Real BX0;
      Real BY0;
      Real BZ0;
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

      static Real rndRho, rndVel[3];
      #pragma omp threadprivate(rndRho,rndVel)
   } ; // class Fluctuations
} // namespace projects
#endif
