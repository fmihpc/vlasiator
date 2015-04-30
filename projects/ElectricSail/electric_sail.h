/* This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
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
      Real getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const int& popID) const;
      virtual void getParameters();
      virtual bool initialize();
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
                
   protected:
      int popID;
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
            creal& dvx, creal& dvy, creal& dvz,const int& popID) const;
      
      Real getDistribValue(creal& vx,creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,const int& popID) const;

      std::vector<std::array<Real,3>> getV0(creal x,creal y,creal z) const;

      virtual bool rescalesDensity(const int& popID) const;

    }; // class PoissonTest

} // namespace projects

#endif	// POISSON_TEST_H

