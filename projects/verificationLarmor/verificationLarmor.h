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

#ifndef VERIFICATIONLARMOR_H
#define VERIFICATIONLARMOR_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class verificationLarmor: public Project {
   public:
      verificationLarmor();
      virtual ~verificationLarmor();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      
   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz);
      virtual void calcCellParameters(Real* cellParams,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz
      );

      Real BX0;
      Real BY0;
      Real BZ0;
      Real VX0;
      Real VY0;
      Real VZ0;
      Real X0;
      Real Y0;
      Real Z0;
      Real DENSITY;

   }; // class verificationLarmor
} // namespace projects

#endif

