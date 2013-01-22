/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef KHB_H
#define KHB_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class KHB: public Project {
      public:
         KHB();
         virtual ~KHB();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
      protected:
         Real getDistribValue(
            creal& x, creal& z,
            creal& vx, creal& vy, creal& vz);
         Real profile(creal top, creal bottom, creal x, creal z);
         
         enum {
            TOP,
            BOTTOM
         };
         Real rho[2];
         Real T[2];
         Real Vx[2];
         Real Vy[2];
         Real Vz[2];
         Real Bx[2];
         Real By[2];
         Real Bz[2];
         Real lambda;
         Real amp;
         Real offset;
         Real transitionWidth;
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class KHB
} // namespace
#endif
