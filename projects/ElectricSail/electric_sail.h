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

      Real ionCloudRadius;
      
      Real tether_x;          /**< Electric sail tether x-position.*/
      Real tether_y;          /**< Electric sail tether y-position.*/
      Real tetherUnitCharge;  /**< Unit charge per meter of the tether in Coulombs.*/
      Real tetherVoltage;     /**< Electric sail tether voltage.*/

      void tetherElectricField(Real* x,Real* E) const;

      virtual void calcCellParameters(Real* cellParams,creal& t);

      virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const int& popID) const;
      
      Real getDistribValue(creal& vx,creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,const int& popID) const;

      std::vector<std::array<Real,3>> getV0(creal x,creal y,creal z) const;

    }; // class PoissonTest

} // namespace projects

#endif	// POISSON_TEST_H

