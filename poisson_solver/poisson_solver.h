/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver.h
 * Author: sandroos
 *
 * Created on January 14, 2015, 1:42 PM
 */

#ifndef POISSON_SOLVER_H
#define	POISSON_SOLVER_H

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../object_factory.h"
#include "../spatial_cell.hpp"

namespace poisson {

    class PoissonSolver {
    public:
        PoissonSolver();
        virtual ~PoissonSolver();

        virtual Real error(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        virtual bool initialize();
        virtual bool finalize();
        virtual bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) = 0;

    protected:
        
    };

    struct Poisson {
        static int RHOQ_TOT;
        static int PHI;

        static ObjectFactory<PoissonSolver> solvers;
        static PoissonSolver* solver;
        static std::string solverName;
    };

    bool initialize(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
    bool finalize();
    bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

} // namespace poisson

#endif	// POISSON_SOLVER_H

