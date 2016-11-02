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
 * File:   poisson_solver_jacobi.cpp
 * Author: sandroos
 *
 * Created on January 16, 2015.
 */

#include <cstdlib>
#include <iostream>

#include "../grid.h"

#include "poisson_solver_jacobi.h"

using namespace std;

namespace poisson {

   static const int RED   = 0;
   static const int BLACK = 1;

   PoissonSolver* makeJacobi() {
      return new PoissonSolverJacobi();
   }

   PoissonSolverJacobi::PoissonSolverJacobi(): PoissonSolver() { }
   
   PoissonSolverJacobi::~PoissonSolverJacobi() { }
   
   bool PoissonSolverJacobi::initialize() {
      bool success = true;
      return success;
   }
   
   bool PoissonSolverJacobi::finalize() {
      bool success = true;
      return success;
   }
   
   bool PoissonSolverJacobi::calculateElectrostaticField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      #warning Jacobi solver does not calculate electric field
      return false;
   }
   
   void PoissonSolverJacobi::evaluate(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   const std::vector<CellID>& cells) {      
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         
         // Skip cells on domain boundaries:
         if (mpiGrid[cellID]->sysBoundaryFlag != 1) continue;
         
         // Fetch data
         const Real rho_q = mpiGrid[cellID]->parameters[CellParams::RHOQ_TOT];
         Real phi_111 = mpiGrid[cellID]->parameters[CellParams::PHI_TMP];
         
         // Calculate cell i/j/k indices
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cellID);
         CellID nbrID;

         // +/- x face neighbor potential
         indices[0] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);
         Real phi_011 = mpiGrid[nbrID]->parameters[CellParams::PHI_TMP];
         indices[0] += 2; nbrID = mpiGrid.mapping.get_cell_from_indices(indices,0);
         Real phi_211 = mpiGrid[nbrID]->parameters[CellParams::PHI_TMP];
         indices[0] -= 1;

         // +/- y face neighbor potential
         indices[1] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_101 = mpiGrid[nbrID]->parameters[CellParams::PHI_TMP];
         indices[1] += 2; nbrID = mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_121 = mpiGrid[nbrID]->parameters[CellParams::PHI_TMP];
         indices[1] -= 1;

         // +/- z face neighbor potential
         /*indices[2] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_110 = mpiGrid[nbrID]->parameters[CellParams::PHI_TMP];
         indices[2] += 2; mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_112 = mpiGrid[nbrID]->parameters[CellParams::PHI_TMP];
         indices[2] -= 1;*/
         Real phi_110 = phi_111; Real phi_112 = phi_111;

         Real DX2 = mpiGrid[cellID]->parameters[CellParams::DX]*mpiGrid[cellID]->parameters[CellParams::DX];
         Real DY2 = mpiGrid[cellID]->parameters[CellParams::DY]*mpiGrid[cellID]->parameters[CellParams::DY];
         Real DZ2 = mpiGrid[cellID]->parameters[CellParams::DZ]*mpiGrid[cellID]->parameters[CellParams::DZ];
         Real factor = 2*(1/DX2 + 1/DY2 + 1/DZ2);
         Real rhs = ((phi_011+phi_211)/DX2 + (phi_101+phi_121)/DY2 + (phi_110+phi_112)/DZ2 + rho_q)/factor;
         mpiGrid[cellID]->parameters[CellParams::PHI] = rhs;
      }
   }

   bool PoissonSolverJacobi::solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;
      phiprof::start("Poisson Solver");
      
      // Update charge density
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_RHOQ_TOT,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);

      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
      
      do {
         if (iterate(mpiGrid) == false) success = false;

         // Evaluate the error in potential solution and reiterate if necessary         
         break;
      } while (true);

      phiprof::stop("Poisson Solver");

      return success;
   }

   bool PoissonSolverJacobi::iterate(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

      bool success = true;
      
      // Swap new and temporary potentials:
      vector<CellID> cells = mpiGrid.get_cells();
      for (size_t c=0; c<cells.size(); ++c) {
         mpiGrid[cells[c]]->parameters[CellParams::PHI_TMP] = mpiGrid[cells[c]]->parameters[CellParams::PHI];
      }
      cells = mpiGrid.get_remote_cells_on_process_boundary(POISSON_NEIGHBORHOOD_ID);
      for (size_t c=0; c<cells.size(); ++c) {
         mpiGrid[cells[c]]->parameters[CellParams::PHI_TMP] = mpiGrid[cells[c]]->parameters[CellParams::PHI];
      }
      
      // Compute new potential on process boundary cells
      cells = mpiGrid.get_local_cells_on_process_boundary(POISSON_NEIGHBORHOOD_ID);
      evaluate(mpiGrid,cells);

      // Exchange new potential values on process boundaries
      mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
         
      // Compute new potential on inner cells
      cells = mpiGrid.get_local_cells_not_on_process_boundary(POISSON_NEIGHBORHOOD_ID);
      evaluate(mpiGrid,cells);

      // Wait for MPI transfers to complete
      mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);

      return success;
   }
   
} // namespace poisson
