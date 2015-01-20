/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver_sor.h
 * Author: sandroos
 *
 * Created on January 15, 2015, 12:45 PM
 */

#include <cstdlib>
#include <iostream>

#include "../grid.h"

#include "poisson_solver_sor.h"

using namespace std;

namespace poisson {

   static const int RED   = 0;
   static const int BLACK = 1;

   // Static initializer for the solver that adds it to solver object factory
   
   PoissonSolver* makeSOR() {
      return new PoissonSolverSOR();
   }
   
   class InitSOR {
   public:
      InitSOR() {
         Poisson::solvers.add("SOR",makeSOR);
      }
   };
   static InitSOR initSor;

   PoissonSolverSOR::PoissonSolverSOR(): PoissonSolver() { }
   
   PoissonSolverSOR::~PoissonSolverSOR() { }
   
   bool PoissonSolverSOR::initialize() {
      bool success = true;
      return success;
   }
   
   bool PoissonSolverSOR::finalize() {
      bool success = true;
      return success;
   }
   
   void PoissonSolverSOR::evaluate(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   const std::vector<CellID>& cells,const int& oddness) {
      
      const Real weight = 1.5;

      #pragma omp for
      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         
         // Skip cells on domain boundaries:
         if (mpiGrid[cellID]->sysBoundaryFlag != 1) continue;

         // Calculate cell i/j/k indices
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cellID);
         
         // Check that the cell has the correct coloring (red or black)
         if ((indices[0] + indices[1]%2 + indices[2]%2) % 2 != oddness) {
            continue;
         }

         // Fetch data
         const Real rho_q = mpiGrid[cellID]->parameters[CellParams::RHOQ_TOT];
         Real phi_111 = mpiGrid[cellID]->parameters[CellParams::PHI];

         CellID nbrID;

         // +/- x face neighbor potential
         indices[0] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_011 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[0] += 2; nbrID = mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_211 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[0] -= 1;
         
         // +/- y face neighbor potential
         indices[1] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_101 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[1] += 2; nbrID = mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_121 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[1] -= 1;
         
         // +/- z face neighbor potential
         indices[2] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);
         Real phi_110 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[2] += 2; nbrID = mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_112 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[2] -= 1;
//         Real phi_110 = phi_111; Real phi_112 = phi_111;
         
         Real DX2 = mpiGrid[cellID]->parameters[CellParams::DX]*mpiGrid[cellID]->parameters[CellParams::DX];
         Real DY2 = mpiGrid[cellID]->parameters[CellParams::DY]*mpiGrid[cellID]->parameters[CellParams::DY];
         Real DZ2 = mpiGrid[cellID]->parameters[CellParams::DZ]*mpiGrid[cellID]->parameters[CellParams::DZ];
         Real factor = 2*(1/DX2 + 1/DY2 + 1/DZ2);
         
         Real rhs = ((phi_011+phi_211)/DX2 + (phi_101+phi_121)/DY2 + (phi_110+phi_112)/DZ2 + rho_q)/factor;         
         Real correction = rhs - phi_111;
         mpiGrid[cellID]->parameters[CellParams::PHI] = phi_111 + weight*correction;
      }      
   }
   
   bool PoissonSolverSOR::solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;
      
      // Update charge density
      phiprof::start("MPI");
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_RHOQ_TOT,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
      phiprof::stop("MPI");

      do {
         #pragma omp parallel
           {
              // Solve red cells first, the black cells
              if (solve(mpiGrid,RED  ) == false) success = false;
              if (solve(mpiGrid,BLACK) == false) success = false;
           }

         // Evaluate the error in potential solution and reiterate if necessary
         error(mpiGrid);

         break;
      } while (true);

      return success;
   }

   bool PoissonSolverSOR::solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                const int& oddness) {
      bool success = true;
      
      // Compute new potential on process boundary cells
      phiprof::start("Evaluate potential");
      vector<CellID> cells = mpiGrid.get_local_cells_on_process_boundary(POISSON_NEIGHBORHOOD_ID);
      evaluate(mpiGrid,cells,oddness);
      phiprof::stop("Evaluate potential");

      // Exchange new potential values on process boundaries
      phiprof::start("MPI");
      //mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      phiprof::stop("MPI");

      // Compute new potential on inner cells
      phiprof::start("Evaluate potential");
      cells = mpiGrid.get_local_cells_not_on_process_boundary(POISSON_NEIGHBORHOOD_ID);
      evaluate(mpiGrid,cells,oddness);
      phiprof::stop("Evaluate potential");

      // Wait for MPI transfers to complete
      phiprof::start("MPI");
      //mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      phiprof::stop("MPI");

      return success;
   }
   
} // namespace poisson
