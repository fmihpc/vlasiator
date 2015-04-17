/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef SHOCKTEST_H
#define SHOCKTEST_H

#include <vector>

#include "../../definitions.h"
#include "../../spatial_cell.hpp"
#include "../project.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class Shocktest: public TriAxisSearch {
    public:
      Shocktest(); // Constructor
      virtual ~Shocktest(); // Destructor
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      
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
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell);
      
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const int& popID
                                        );
         
      virtual std::vector<std::array<Real, 3>> getV0(
                                                     creal x,
                                                     creal y,
                                                     creal z
                                                    );
         
   }; // Class Shocktest

} // Namespace projects
#endif

