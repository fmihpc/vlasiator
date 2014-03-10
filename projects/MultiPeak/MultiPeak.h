/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/


#ifndef MULTIPEAK_H
#define MULTIPEAK_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class MultiPeak: public TriAxisSearch {
      public:
         MultiPeak();
         virtual ~MultiPeak();
         
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
         virtual void setCellBackgroundField(SpatialCell* cell);
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz
         );
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );
         virtual vector<std::array<Real, 3>> getV0(
            creal x,
            creal y,
            creal z
         );
         
         Real rho[2];
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
//          Real rho1PertRelAmp;
//          Real rho2PertRelAmp;
//          Real Vx1PertAbsAmp;
//          Real Vy1PertAbsAmp;
//          Real Vz1PertAbsAmp;
//          Real Vx2PertAbsAmp;
//          Real Vy2PertAbsAmp;
//          Real Vz2PertAbsAmp;
         Real lambda;
         uint nVelocitySamples;
         uint seed;
         
         char rngStateBuffer[256];
         random_data rngDataBuffer;
         
//          #ifdef _AIX
//          static int64_t rndRho, rndVel1[3], rndVel2[3];
//          #else
//          static int32_t rndRho, rndVel1[3], rndVel2[3];
//          #endif
//          #pragma omp threadprivate(rndRho,rndVel)
         
   }; // class MultiPeak
} //  namespace projects



#endif

