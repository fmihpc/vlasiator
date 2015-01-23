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

    struct CellCache2D {
       Real* parameters[5];
       Real*& operator[](const int& i) {return parameters[i];}
    };

    struct CellCache3D {
       Real* parameters[7];
       Real*& operator[](const int& i) {return parameters[i];}
    };

    class PoissonSolver {
    public:
       PoissonSolver();
       virtual ~PoissonSolver();

       virtual Real error(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
       virtual Real error3D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
       virtual bool initialize();
       virtual bool finalize();
       virtual bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) = 0;

    protected:
        
    };

   /** Wrapper for all variables needed by Poisson solvers.*/
   struct Poisson {
      static int RHOQ_TOT;
      static int PHI;

      static ObjectFactory<PoissonSolver> solvers; /**< Container for all existing Poisson solvers.*/
      static PoissonSolver* solver;                /**< Poisson solver used in the simulation.*/
      static std::string solverName;               /**< Name of the Poisson solver in use.*/
      
      static uint maxIterations;                   /**< Maximum number of iterations allowed, only 
                                                    * has an effect on iterative solvers.*/
      static Real minRelativePotentialChange;      /**< Iterative solvers keep on iterating the solution 
                                                    * until the change in potential during successive 
                                                    * iterations is less than this value.*/
      
      static std::vector<Real*> localCellParams;   /**< Pointers to spatial cell parameters, order 
						    * is the same as in getLocalCells() vector.*/
      
      static void cacheCellParameters(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
				      const std::vector<CellID>& cells);
   };

   bool initialize(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
   bool finalize();
   bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

} // namespace poisson

#endif	// POISSON_SOLVER_H

