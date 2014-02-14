/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef SHOCKTEST_H
#define SHOCKTEST_H

#include "../../definitions.h"
#include "../../spatial_cell.hpp"
//#include "../projects_common.h"
#include "../project.h"
#include "../projectIsotropicMaxwellian.h"

//#include "../projects_vlasov_acceleration.h"

#include "dccrg.hpp"

namespace projects {
   class Shocktest: public IsotropicMaxwellian {
      public:
         Shocktest(); // Constructor
         virtual ~Shocktest(); // Destructor

         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         virtual void setCell(SpatialCell* cell);

      protected:
         enum {
            LEFT,
            RIGHT
         };
         Real rho[2];
         Real T[2];
         Real Vx[2];
         Real Vy[2];
         Real Vz[2];
         Real Bx[2];
         Real By[2];
         Real Bz[2];
         uint nSpaceSamples;
         uint nVelocitySamples;

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
         //virtual std::vector<uint> findBlocksToInitialize(SpatialCell* cell);
         virtual Real getV0(
            creal x,
            creal y,
            creal z,
            cuint component
         );

         // Couldn't find an explanation for this        
//         template<typename UINT,typename REAL> void calcAccFaceX(
//            REAL& ax, REAL& ay, REAL& az,
//            const UINT& I, const UINT& J, const UINT& K,
//            const REAL* const cellParams,
//            const REAL* const blockParams,
//            const REAL* const cellBVOLDerivatives
//         ) {
//            lorentzForceFaceX(ax,ay,az,I,J,K,cellParams,blockParams,cellBVOLDerivatives);
//         }
//         
//         template<typename UINT,typename REAL> void calcAccFaceY(
//            REAL& ax, REAL& ay, REAL& az,
//            const UINT& I, const UINT& J, const UINT& K,
//            const REAL* const cellParams,
//            const REAL* const blockParams,
//            const REAL* const cellBVOLDerivatives
//         ) {
//            lorentzForceFaceY(ax,ay,az,I,J,K,cellParams,blockParams,cellBVOLDerivatives);
//         }
//         
//         template<typename UINT,typename REAL> void calcAccFaceZ(
//            REAL& ax, REAL& ay, REAL& az,
//            const UINT& I, const UINT& J, const UINT& K,
//            const REAL* const cellParams,
//            const REAL* const blockParams,
//            const REAL* const cellBVOLDerivatives
//         ) {
//            lorentzForceFaceZ(ax,ay,az,I,J,K,cellParams,blockParams,cellBVOLDerivatives);
//         }
   }; // Class Shocktest

} // Namespace projects
#endif

