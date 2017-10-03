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
 * 
 * File:   electric_sail.h
 * Author: sandroos
 *
 * Created on March 3, 2015
 */

#ifndef ELECTRIC_SAIL_H
#define	ELECTRIC_SAIL_H

#include <vector>
#include "../project.h"
#include "../projectTriAxisSearch.h"

namespace projects {

   struct Population {
      Real rho;
      Real T[3];
      Real V[3];
      
      Population(const Real& rho,const Real& Tx,const Real& Ty,const Real& Tz,
                 const Real& Vx,const Real& Vy,const Real& Vz);
   };

   class ElectricSail: public TriAxisSearch {
   public:
      ElectricSail();
      virtual ~ElectricSail();
        
      static void addParameters();
      Real getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const uint popID) const;
      virtual void getParameters();
      virtual bool initialize();
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
                
   protected:
      uint popID;
      std::vector<Population> populations;

      bool addParticleCloud;     /**< If true, a charge-neutralising particle cloud is added around the tether.*/
      Real particleCloudRadius;  /**< Radius of the particle cloud.*/
      
      Real tether_x;             /**< Electric sail tether x-position.*/
      Real tether_y;             /**< Electric sail tether y-position.*/
      Real tetherChargeRiseTime; /**< Time when tether charge reaches its maximum value.
                                  * Only has an effect if ElectricSail::timeDependentCharge is true.*/
      Real tetherUnitCharge;     /**< Unit charge per meter of the tether in Coulombs.*/
      bool timeDependentCharge;  /**< If true, tether charge is time-dependent.*/
      bool useBackgroundField;   /**< If true, then tether electric field is calculated as a constant 
                                  * background field. If false, a charge placed near the tether position
                                  * is used instead.*/

      void tetherElectricField(Real* x,Real* E) const;

      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);

      virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const uint popID) const;
      
      Real getDistribValue(creal& vx,creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,const uint popID) const;

      std::vector<std::array<Real,3>> getV0(creal x,creal y,creal z, const uint popID) const;

      virtual bool rescalesDensity(const uint popID) const;

    }; // class PoissonTest

} // namespace projects

#endif	// POISSON_TEST_H

