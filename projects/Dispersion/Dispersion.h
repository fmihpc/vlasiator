/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const int& popID
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
