/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver_sor.h
 * Author: sandroos
 *
 * Created on January 15, 2015, 12:45 PM
 */

#ifndef POISSON_SOLVER_SOR_H
#define	POISSON_SOLVER_SOR_H

#include <vector>

#include "poisson_solver.h"

namespace poisson {

    class PoissonSolverSOR: public PoissonSolver {
    public:
        PoissonSolverSOR();
        ~PoissonSolverSOR();
        
        bool initialize();
        bool finalize();
        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        
    private:
        void evaluate(std::vector<poisson::CellCache3D>& cellPointers,const int& cellColor);
       
        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                   const int& oddness);
    };

    PoissonSolver* makeSOR();
   
} // namespace poisson

#endif	// POISSON_SOLVER_SOR_H

