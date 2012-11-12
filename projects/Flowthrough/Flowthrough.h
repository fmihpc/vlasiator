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

#ifndef FLOWTHROUGH_H
#define FLOWTHROUGH_H

#include "../../definitions.h"
#include "../project.h"


namespace projects {
   class Flowthrough: public Project {
      public:
         Flowthrough();
         virtual ~Flowthrough();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         
         
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         
         Real rho;
         Real T;
         Real Bx;
         Real By;
         Real Bz;
         Real V0[3];
         uint nSpaceSamples;
         uint nVelocitySamples;
   }; // class Flowthrough
} // namespace projects

#endif
