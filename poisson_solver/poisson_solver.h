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
       spatial_cell::SpatialCell* cell;
       Real* parameters[5];
       Real*& operator[](const int& i) {return parameters[i];}
    };

    struct CellCache3D {
#warning TEMP remove
       CellID cellID;
       spatial_cell::SpatialCell* cell;
       Real* parameters[7];
       Real*& operator[](const int& i) {return parameters[i];}
    };

    class PoissonSolver {
    public:
       PoissonSolver();
       virtual ~PoissonSolver();

       // ***** DECLARATIONS OF VIRTUAL MEMBER FUNCTIONS ***** //

       virtual bool calculateBackgroundField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                             const std::vector<CellID>& cells);
       virtual bool calculateChargeDensity(spatial_cell::SpatialCell* cell);
       virtual bool calculateElectrostaticField2D(const std::vector<poisson::CellCache3D>& cells);
       virtual bool calculateElectrostaticField3D(const std::vector<poisson::CellCache3D>& cells);
       virtual Real error(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
       virtual Real error3D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
       virtual bool initialize();
       virtual bool finalize();
       
       // ***** DECLARATIONS OF PURE VIRTUAL MEMBER FUNCTIONS ***** //

       /** Calculate electric field on all non-system boundary spatial cells.
        * @param mpiGrid Parallel grid library.
        * @return If true, electric field was calculated successfully.*/
       virtual bool calculateElectrostaticField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) = 0;

       /** Solve Poisson equation on all non-system boundary spatial cells.
        * @param mpiGrid Parallel grid library.
        * @return If true, Poisson equation was successfully solved.*/
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

      static bool clearPotential;                  /**< If true, then potential is cleared each timestep 
                                                    * before solving Poisson's equation. Otherwise the old 
                                                    * potential is used as an initial guess.*/
      static bool is2D;                            /**< If true, then system is two-dimensional, i.e., 
                                                    * electrostatic potential and electric field is solved 
                                                    * in xy-plane.*/
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

