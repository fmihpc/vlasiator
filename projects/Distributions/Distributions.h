/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












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
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell);
    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz
                          );
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const int& popID
                                        );
      virtual vector<std::array<Real, 3>> getV0(
                                                creal x,
                                                creal y,
                                                creal z
                                               );
      
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
      uint nVelocitySamples;
   }; // class Distributions
} //  namespace projects



#endif

