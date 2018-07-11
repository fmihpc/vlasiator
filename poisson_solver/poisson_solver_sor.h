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

   const unsigned int SOR_VARS = 0;
   
    class PoissonSolverSOR: public PoissonSolver {
    public:
        PoissonSolverSOR();
        ~PoissonSolverSOR();
        
        bool calculateElectrostaticField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        bool initialize();
        bool finalize();
        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        
    private:

        Real bndryCellParams[CellParams::N_SPATIAL_CELL_PARAMS];

        bool boundaryConds(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
        void boundaryConds(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                        std::vector<poisson::CellCache3D<SOR_VARS> >& cells);
        
        void cachePointers2D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,std::vector<poisson::CellCache3D<SOR_VARS> >& redCache,
                             std::vector<poisson::CellCache3D<SOR_VARS> >& blackCache);
        void cachePointers3D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,std::vector<poisson::CellCache3D<SOR_VARS> >& redCache,
                             std::vector<poisson::CellCache3D<SOR_VARS> >& blackCache);
        void evaluate2D(std::vector<poisson::CellCache3D<SOR_VARS> >& cellPointers,const int& cellColor);
        void evaluate3D(std::vector<poisson::CellCache3D<SOR_VARS> >& cellPointers,const int& cellColor);

        void (*evaluator)(std::vector<poisson::CellCache3D<SOR_VARS> >& cellPointers,const int& cellColor);

        bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                   const int& oddness);

    };

    PoissonSolver* makeSOR();
   
} // namespace poisson

#endif	// POISSON_SOLVER_SOR_H

