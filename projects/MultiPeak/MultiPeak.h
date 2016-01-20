/*
This file is part of Vlasiator.

Copyright 2011-2015 Finnish Meteorological Institute

*/

#ifndef MULTIPEAK_H
#define MULTIPEAK_H

#include <vector>

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
      virtual void setActivePopulation(const int& popID);
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                          const int& popID) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const int& popID) const;
      virtual std::vector<std::array<Real, 3> > getV0(
                                                      creal x,
                                                      creal y,
                                                      creal z
                                                     ) const;
      int popID;
      int numberOfPopulations;
      std::vector<Real> rho;
      static std::vector<Real> rhoRnd; //static as it has to be threadprivate
      #pragma omp threadprivate(rhoRnd)       
      std::vector<Real> Tx;
      std::vector<Real> Ty;
      std::vector<Real> Tz;
      std::vector<Real> Vx;
      std::vector<Real> Vy;
      std::vector<Real> Vz;
      Real Bx;
      Real By;
      Real Bz;
      Real dBx;
      Real dBy;
      Real dBz;
      Real magXPertAbsAmp;
      Real magYPertAbsAmp;
      Real magZPertAbsAmp;
      std::vector<Real> rhoPertAbsAmp;
      Real lambda;
      uint nVelocitySamples;
      bool useMultipleSpecies;      /**< If true, then each peak is a separate particle species.
                                     * Defaults to false.*/
      
      enum densitymodel {
         Uniform,
         TestCase
      } densityModel;

   }; // class MultiPeak
} //  namespace projects

#endif

