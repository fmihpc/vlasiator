/*
This file is part of Vlasiator.

Copyright 2011, 2012, 2015 Finnish Meteorological Institute

*/

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class Template: public TriAxisSearch {
    public:
      Template();
      virtual ~Template();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const int& popID
                                        );
      
    protected:
      virtual std::vector<std::array<Real, 3>> getV0(
                                                     creal x,
                                                     creal y,
                                                     creal z
                                                    ) const;
      
      Real param;
   }; // class Template
} // namespace projects

#endif

