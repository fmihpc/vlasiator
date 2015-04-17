/* This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * File:   poisson_test.h
 * Author: sandroos
 *
 * Created on January 14, 2015, 3:17 PM
 */

#ifndef POISSON_TEST_H
#define POISSON_TEST_H

#include <vector>
#include "../project.h"

namespace projects {

    class PoissonTest: public Project {
    public:
       PoissonTest();
       virtual ~PoissonTest();
       
       static void addParameters();
       virtual void getParameters();
       virtual bool initialize();
       virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell);

     protected:
       virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
       
       virtual Real calcPhaseSpaceDensity(
                                          creal& x, creal& y, creal& z,
                                          creal& dx, creal& dy, creal& dz,
                                          creal& vx, creal& vy, creal& vz,
                                          creal& dvx, creal& dvy, creal& dvz);
       
       virtual std::vector<uint> findBlocksToInitialize(spatial_cell::SpatialCell* cell);
       
    }; // class PoissonTest

} // namespace projects

#endif	// POISSON_TEST_H

