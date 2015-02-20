/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Diffusion: public Project {
      public:
         Diffusion();
         virtual ~Diffusion();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         /*! set background field, should set it for all cells */
         virtual void setCellBackgroundField(SpatialCell* cell);
      
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz
         );
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const int& popID
         );
         
         Real B0;
         Real DENSITY;
         Real TEMPERATURE;
         Real SCA_X;
         Real SCA_Y;
         uint nSpaceSamples;
         uint nVelocitySamples;
   } ; // class Diffusion
} // namespace projects

#endif
