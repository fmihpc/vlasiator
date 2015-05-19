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
         int numberOfPopulations;
         vector<Real> rho;
         static vector<Real> rhoRnd; //static as it has to be threadprivate
#pragma omp threadprivate(rhoRnd)       
         vector<Real> Tx;
         vector<Real> Ty;
         vector<Real> Tz;
         vector<Real> Vx;
         vector<Real> Vy;
         vector<Real> Vz;
         Real Bx;
         Real By;
         Real Bz;
         Real dBx;
         Real dBy;
         Real dBz;
         Real magXPertAbsAmp;
         Real magYPertAbsAmp;
         Real magZPertAbsAmp;
         vector<Real> rhoPertAbsAmp;
         Real lambda;
         uint nVelocitySamples;
   }; // class MultiPeak
} //  namespace projects



#endif

