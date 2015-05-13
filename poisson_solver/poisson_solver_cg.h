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

   namespace cgvar {
      enum Variable {
         PHI,
         B,          /**< Charge density multiplied by dx2/epsilon0.*/
         R,
         A_TIMES_P,  /**< Matrix A times P.*/
         SIZE
      };
   }
   
   namespace cgglobal {
      enum Variable {
         ALPHA,
         R_T_R,
         P_T_A_P,
         BETA,
         R_MAX,
         SIZE
      };
   }

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

        void bvalue(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                    std::vector<CellCache3D<cgvar::SIZE> >& cells);
        
        void cachePointers2D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,
                             std::vector<poisson::CellCache3D<cgvar::SIZE> >& cellCache);
        /*void cachePointers3D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,std::vector<poisson::CellCache3D>& redCache,
                             std::vector<poisson::CellCache3D>& blackCache);*/

        bool calculateAlpha();
        void calculateAlpha(CellCache3D<cgvar::SIZE>& cell,Real& mySum0,Real& mySum1);
        bool startIteration(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        bool startIteration(std::vector<CellCache3D<cgvar::SIZE> >& cells);
        bool update_p(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        bool update_x_r();

        Real globalVariables[cgglobal::SIZE];
        int iterations;
   };

   PoissonSolver* makeCG();

} // namespace poisson

#endif	// POISSON_SOLVER_CG_H

