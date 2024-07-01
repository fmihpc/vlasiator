/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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


#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class Distributions: public TriAxisSearch {
    public:
      Distributions();
      virtual ~Distributions();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
      );
    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           const uint popID
                          ) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const uint popID
                                        ) const;
      virtual std::vector<std::array<Real, 3> > getV0(
                                                      creal x,
                                                      creal y,
                                                      creal z,
                                                      const uint popID
                                                     ) const;

      Real rho[2];
      Real rhoRnd[2];
      Real Tx[2];
      Real Ty[2];
      Real Tz[2];
      Real Vx[2];
      Real Vy[2];
      Real Vz[2];
      Real Bx;
      Real By;
      Real Bz;
      Real dBx;
      Real dBy;
      Real dBz;
      Real magXPertAbsAmp;
      Real magYPertAbsAmp;
      Real magZPertAbsAmp;
      Real rhoPertAbsAmp[2];
      //          Real Vx1PertAbsAmp;
      //          Real Vy1PertAbsAmp;
      //          Real Vz1PertAbsAmp;
      //          Real Vx2PertAbsAmp;
      //          Real Vy2PertAbsAmp;
      //          Real Vz2PertAbsAmp;
      Real lambda;
   }; // class Distributions
} //  namespace projects



#endif

