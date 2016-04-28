/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
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

