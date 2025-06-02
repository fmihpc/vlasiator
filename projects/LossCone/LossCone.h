/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * 2017-2025 University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
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

#ifndef LOSSCONE_H
#define LOSSCONE_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {

   struct LossConeSpeciesParameters {
      Real DENSITY;
      Real TEMPERATUREX;
      Real TEMPERATUREY;
      Real TEMPERATUREZ;
      Real densityPertRelAmp;
      Real velocityPertAbsAmp;
      Real V0[3];
      Real muLimit;
   };

   class LossCone: public TriAxisSearch {
   public:
      LossCone();
      virtual ~LossCone();
      
      virtual bool initialize(void) override;
      static void addParameters(void);
      virtual void getParameters(void) override;
      virtual void setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
         FsGrid< fsgrids::technical, 2>& technicalGrid
      ) override;
      virtual std::vector<std::array<Real, 3> > getV0(
         creal x,
         creal y,
         creal z,
         const uint popID
      ) const override;

      virtual Realf fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                  const uint popID,
                                  const uint nRequested) const override;
      virtual Realf probePhaseSpace(spatial_cell::SpatialCell *cell,
                                    const uint popID,
                                    Real vx_in, Real vy_in, Real vz_in) const override;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) override;
      
      Real BX0;
      Real BY0;
      Real BZ0;
      Real magXPertAbsAmp;
      Real magYPertAbsAmp;
      Real magZPertAbsAmp;
      uint seed;
      std::vector<LossConeSpeciesParameters> speciesParams;

      static Real rndRho, rndVel[3];
      #pragma omp threadprivate(rndRho,rndVel)
   } ; // class LossCone
} // namespace projects
#endif
