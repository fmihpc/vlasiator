/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * File:   poisson_solver_jacobi.h
 * Author: sandroos
 *
 * Created on January 16, 2015.
 */

#ifndef POISSON_SOLVER_JACOBI_H
#define	POISSON_SOLVER_JACOBI_H

#include <vector>

#include "poisson_solver.h"

namespace poisson {

    class PoissonSolverJacobi: public PoissonSolver {
    public:
        PoissonSolverJacobi();
        ~PoissonSolverJacobi();
        
        bool calculateElectrostaticField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        bool initialize();
        bool finalize();
        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        
    private:
        void evaluate(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const std::vector<CellID>& cells);
        bool iterate(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
    };
   
    PoissonSolver* makeJacobi();
    
} // namespace poisson

#endif	// POISSON_SOLVER_JACOBI_H

