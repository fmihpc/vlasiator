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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
        uint iterations;
   };

   PoissonSolver* makeCG();

} // namespace poisson

#endif	// POISSON_SOLVER_CG_H

