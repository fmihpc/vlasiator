/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver_cg.h
 * Author: sandroos
 *
 * Created on January 15, 2015, 12:45 PM
 */

#ifndef POISSON_SOLVER_CG_H
#define	POISSON_SOLVER_CG_H

#include <vector>

#include "poisson_solver.h"

namespace poisson {

    class PoissonSolverCG: public PoissonSolver {
    public:
        PoissonSolverCG();
        ~PoissonSolverCG();
        
        bool calculateElectrostaticField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        bool initialize();
        bool finalize();
        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        
    private:

        Real bndryCellParams[CellParams::N_SPATIAL_CELL_PARAMS];

        void cachePointers2D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,std::vector<poisson::CellCache3D>& redCache,
                             std::vector<poisson::CellCache3D>& blackCache);
        void cachePointers3D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,std::vector<poisson::CellCache3D>& redCache,
                             std::vector<poisson::CellCache3D>& blackCache);
        void evaluate2D(std::vector<poisson::CellCache3D>& cellPointers,const int& cellColor);
        void evaluate3D(std::vector<poisson::CellCache3D>& cellPointers,const int& cellColor);
       
        void (*evaluator)(std::vector<poisson::CellCache3D>& cellPointers,const int& cellColor);
        
        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                   const int& oddness);
    };

    PoissonSolver* makeSOR();
    
   
} // namespace poisson

#endif	// POISSON_SOLVER_CG_H

