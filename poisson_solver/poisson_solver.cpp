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
#include <omp.h>

#include "../common.h"
#include "../logger.h"
#include "../mpiconversion.h"
#include "../grid.h"
#include "../spatial_cell.hpp"       

#include "poisson_solver.h"
#include "poisson_solver_jacobi.h"
#include "poisson_solver_sor.h"

using namespace std;

extern Logger logFile;

namespace poisson {

   // ***** INITIALIZE STATIC VARIABLES ***** //   
   int Poisson::RHOQ_TOT = CellParams::RHOQ_TOT;
   int Poisson::PHI = CellParams::PHI;
   ObjectFactory<PoissonSolver> Poisson::solvers;
   PoissonSolver* Poisson::solver = NULL;
   string Poisson::solverName;
   uint Poisson::maxIterations;
   Real Poisson::minRelativePotentialChange;
   vector<Real*> Poisson::localCellParams;

   void Poisson::cacheCellParameters(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
				     const std::vector<CellID>& cells) {
      // NOTE: This is surprisingly slow as compared to the 
      // similar cache-function in poisson_solver_sor.cpp
      
      // Clear old cache
      Poisson::localCellParams.clear();
      Poisson::localCellParams.resize(cells.size());

      // Fetch pointers
      for (size_t c=0; c<cells.size(); ++c) {
	 Poisson::localCellParams[c] = mpiGrid[cells[c]]->parameters;
      }
   }

   // ***** DEFINITION OF POISSON SOLVER BASE CLASS ***** //

   PoissonSolver::PoissonSolver() { }

   PoissonSolver::~PoissonSolver() { }

   bool PoissonSolver::initialize() {return true;}

   bool PoissonSolver::finalize() {return true;}

   Real PoissonSolver::error(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::start("Potential Change");

      Real* maxError = new Real[omp_get_max_threads()];
      const Real epsilon = 1e-100;

      #pragma omp parallel
        {
	   // Each thread evaluates how much the potential changed as 
	   // compared to the previous iteration and stores the value in 
	   // threadMaxError. After all cells have been processed, the 
	   // per-thread values are stored to array maxError.
           const int tid = omp_get_thread_num();
	   Real threadMaxError = 0;
	   #pragma omp for
	   for (size_t c=0; c<Poisson::localCellParams.size(); ++c) {
	      Real phi     = Poisson::localCellParams[c][CellParams::PHI];
	      Real phi_old = Poisson::localCellParams[c][CellParams::PHI_TMP];
	      
	      Real d_phi     = phi-phi_old;
	      Real d_phi_rel = d_phi / (phi + epsilon);
	      if (fabs(d_phi_rel) > threadMaxError) threadMaxError = fabs(d_phi_rel);
	   }

	   maxError[tid] = threadMaxError;
        }

      // Reduce max local error to master thread (index 0)
      for (int i=1; i<omp_get_max_threads(); ++i) {
         if (maxError[i] > maxError[0]) maxError[0] = maxError[i];
      }
      phiprof::stop("Potential Change");

      // Reduce max error to all MPI processes
      phiprof::start("MPI");
      Real globalMaxError;
      MPI_Allreduce(maxError,&globalMaxError,1,MPI_Type<Real>(),MPI_MAX,MPI_COMM_WORLD);
      delete [] maxError; maxError = NULL;
      phiprof::stop("MPI");

      return globalMaxError;
   }
   
   Real PoissonSolver::error3D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::start("Evaluate Error");

      // DEBUG: Make sure values are up to date
      SpatialCell::set_mpi_transfer_type(spatial_cell::Transfer::CELL_RHOQ_TOT,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(spatial_cell::Transfer::CELL_PHI,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);

      Real localError = 0;
      Real* maxError = new Real[omp_get_max_threads()];
      const vector<CellID>& cells = getLocalCells();

      #pragma omp parallel 
        {
           const int tid = omp_get_thread_num();
           maxError[tid] = 0;

           #pragma omp for reduction(+:localError)
           for (size_t c=0; c<cells.size(); ++c) {
              CellID cellID = cells[c];

              // Skip cells on domain boundaries:
              if (mpiGrid[cellID]->sysBoundaryFlag != 1) continue;

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
              indices[2] -= 1; nbrID =  mpiGrid.mapping.get_cell_from_indices(indices,0);         
              Real phi_110 = mpiGrid[nbrID]->parameters[CellParams::PHI];
              indices[2] += 2; mpiGrid.mapping.get_cell_from_indices(indices,0);         
              Real phi_112 = mpiGrid[nbrID]->parameters[CellParams::PHI];
              indices[2] -= 1;

              // Evaluate error
              Real DX2 = mpiGrid[cellID]->parameters[CellParams::DX]*mpiGrid[cellID]->parameters[CellParams::DX];
              Real DY2 = mpiGrid[cellID]->parameters[CellParams::DY]*mpiGrid[cellID]->parameters[CellParams::DY];
              Real DZ2 = mpiGrid[cellID]->parameters[CellParams::DZ]*mpiGrid[cellID]->parameters[CellParams::DZ];
              Real factor = 2*(1/DX2 + 1/DY2 + 1/DZ2);
              Real rhs = ((phi_011+phi_211)/DX2 + (phi_101+phi_121)/DY2 + (phi_110+phi_112)/DZ2 + rho_q)/factor;

              Real cellError = rhs - phi_111;
              localError += cellError*cellError;
              mpiGrid[cellID]->parameters[CellParams::PHI_TMP] = fabs(cellError);

              if (fabs(cellError) > maxError[tid]) maxError[tid] = fabs(cellError);
           } // for-loop over cells

        } // #pragma omp parallel

      // Compute max error (over threads)
      for (int i=1; i<omp_get_max_threads(); ++i) {
         if (maxError[i] > maxError[0]) maxError[0] = maxError[i];
      }

      Real globalError;
      Real globalMaxError;
      MPI_Allreduce(&localError,&globalError,1,MPI_Type<Real>(),MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(maxError,&globalMaxError,1,MPI_Type<Real>(),MPI_MAX,MPI_COMM_WORLD);      

      if (mpiGrid.get_rank() == 0) {
         cerr << Parameters::tstep << '\t' << sqrt(globalError) << '\t' << globalMaxError << endl;
      }

      delete [] maxError; maxError = NULL;

      phiprof::stop("Evaluate Error");

      return sqrt(globalError);
   }

   // ***** DEFINITIONS OF HIGH-LEVEL DRIVER FUNCTIONS ***** //

   bool initialize(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;

      Poisson::solvers.add("Jacobi",makeJacobi);
      Poisson::solvers.add("SOR",makeSOR);

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
      phiprof::start("Poisson Solver (Total)");

      // If mesh partitioning has changed, recalculate spatial 
      // cell parameters pointer cache:
      if (Parameters::meshRepartitioned == true) {
	 phiprof::start("Cache Cell Parameters");
	 Poisson::cacheCellParameters(mpiGrid,getLocalCells());
	 phiprof::stop("Cache Cell Parameters");
      }

      bool success = true;
      if (Poisson::solver != NULL) {
         if (Poisson::solver->solve(mpiGrid) == false) success = false;
      }
      
      phiprof::stop("Poisson Solver (Total)");
      return success;
   }

} // namespace poisson
