/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
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
                                          creal& dvx, creal& dvy, creal& dvz,
                                          const uint popID) const;
       
       virtual std::vector<uint> findBlocksToInitialize(spatial_cell::SpatialCell* cell);
       
    }; // class PoissonTest

} // namespace projects

#endif	// POISSON_TEST_H

