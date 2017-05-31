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

#ifndef FLOWTHROUGH_H
#define FLOWTHROUGH_H

#include <array>
#include <vector>

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {

   struct FlowthroughSpeciesParameters {
      Real rho;
      Real T;
      Real V0[3];
      uint nSpaceSamples;
      uint nVelocitySamples;
   };

   class Flowthrough: public TriAxisSearch {
    public:
      Flowthrough();
      virtual ~Flowthrough();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;

    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,
                           const uint popID
                          ) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,const uint popID
                                        ) const;
      virtual std::vector<std::array<Real, 3> > getV0(
                                                      creal x,
                                                      creal y,
                                                      creal z,
                                                      const uint popID
                                                     ) const;

      bool emptyBox;               /**< If true, then the simulation domain is empty initially 
                                    * and matter will flow in only through the boundaries.*/

      Real Bx;
      Real By;
      Real Bz;
      std::vector<FlowthroughSpeciesParameters> speciesParams;
   }; // class Flowthrough
} // namespace projects

#endif
