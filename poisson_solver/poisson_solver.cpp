/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver.h
 * Author: sandroos
 *
 * Created on January 14, 2015, 1:42 PM
 */

#include <cstdlib>
#include <iostream>
#include <vector>

#include "../common.h"
#include "../logger.h"
#include "../mpiconversion.h"
#include "../grid.h"
#include "../spatial_cell.hpp"       

#include "poisson_solver.h"

using namespace std;

extern Logger logFile;

namespace poisson {

   // ***** INITIALIZE STATIC VARIABLES ***** //   
   int Poisson::RHOQ_TOT = CellParams::RHOQ_TOT;
   int Poisson::PHI = CellParams::PHI;
   ObjectFactory<PoissonSolver> Poisson::solvers;
   PoissonSolver* Poisson::solver = NULL;
   string Poisson::solverName;

   // ***** DEFINITION OF POISSON SOLVER BASE CLASS ***** //
   
   PoissonSolver::PoissonSolver() { }
   
   PoissonSolver::~PoissonSolver() { }
   
   bool PoissonSolver::initialize() {return true;}
   
   bool PoissonSolver::finalize() {return true;}

   Real PoissonSolver::error(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::start("Evaluate Error");
      
      // DEBUG: Make sure values are up to date
      SpatialCell::set_mpi_transfer_type(spatial_cell::Transfer::CELL_RHOQ_TOT,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(spatial_cell::Transfer::CELL_PHI,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      
      Real totalError2 = 0;
      const vector<CellID>& cells = getLocalCells();
      for (size_t c=0; c<cells.size(); ++c) {
         CellID cellID = cells[c];
         
         // Skip cells on domain boundaries:
         if (mpiGrid[cellID]->sysBoundaryFlag != 1) continue;
         
         // Fetch data:
         // Fetch data
         const Real rho_q = mpiGrid[cellID]->parameters[CellParams::RHOQ_TOT];
         Real phi_111 = mpiGrid[cellID]->parameters[CellParams::PHI];
         
         // Calculate cell i/j/k indices
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cellID);
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
         /*indices[2] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_110 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[2] += 2; mpiGrid.mapping.get_cell_from_indices(indices,0);         
         Real phi_112 = mpiGrid[nbrID]->parameters[CellParams::PHI];
         indices[2] -= 1;*/
         Real phi_110 = phi_111; Real phi_112 = phi_111;

         Real DX2 = mpiGrid[cellID]->parameters[CellParams::DX]*mpiGrid[cellID]->parameters[CellParams::DX];
         Real DY2 = mpiGrid[cellID]->parameters[CellParams::DY]*mpiGrid[cellID]->parameters[CellParams::DY];
         Real DZ2 = mpiGrid[cellID]->parameters[CellParams::DZ]*mpiGrid[cellID]->parameters[CellParams::DZ];
         Real factor = 2*(1/DX2 + 1/DY2 + 1/DZ2);
         Real rhs = ((phi_011+phi_211)/DX2 + (phi_101+phi_121)/DY2 + (phi_110+phi_112)/DZ2 + rho_q)/factor;
         
         Real cellError = rhs - phi_111;
         totalError2 += cellError*cellError;
         
         mpiGrid[cellID]->parameters[CellParams::PHI_TMP] = fabs(cellError);
      }
      
      Real globalError;
      MPI_Reduce(&totalError2,&globalError,1,MPI_Type<Real>(),MPI_SUM,0,MPI_COMM_WORLD);
      cerr << Parameters::tstep << '\t' << sqrt(totalError2) << endl;
      
      phiprof::stop("Evaluate Error");
      return sqrt(totalError2);
   }

   // ***** DEFINITIONS OF HIGH-LEVEL DRIVER FUNCTIONS ***** //

   bool initialize(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;
      
      // Create and initialize the Poisson solver
      Poisson::solver = Poisson::solvers.create(Poisson::solverName);
      if (Poisson::solver == NULL) {
         logFile << "Failed to create Poisson solver '" << Poisson::solverName << "'" << endl << write;
         return false;
      } else {
         logFile << "Successfully initialized Poisson solver '" << Poisson::solverName << "'" << endl << write;
      }

      // Set up the initial state unless the simulation was restarted
      vector<CellID> local_cells = mpiGrid.get_cells();
      for (size_t c=0; c<local_cells.size(); ++c) {
         spatial_cell::SpatialCell* cell = mpiGrid[local_cells[c]];
         cell->parameters[CellParams::PHI] = 0;
      }

      return success;
   }

   bool finalize() {
      bool success = true;
      if (Poisson::solver != NULL) {
         if (Poisson::solver->finalize() == false) success = false;
         delete Poisson::solver;
         Poisson::solver = NULL;
      }
      return success;
   }

   bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::start("Poisson Solver");

      bool success = true;
      if (Poisson::solver != NULL) {
         if (Poisson::solver->solve(mpiGrid) == false) success = false;
      }
      
      phiprof::stop("Poisson Solver");
      return success;
   }

} // namespace poisson
