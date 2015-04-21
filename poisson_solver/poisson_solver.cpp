/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute.
 * 
 * File:   poisson_solver.cpp
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
#include "../object_wrapper.h"

#include "poisson_solver.h"
#include "poisson_solver_jacobi.h"
#include "poisson_solver_sor.h"
#include "poisson_solver_cg.h"

#ifndef NDEBUG
   #define DEBUG_POISSON
#endif

using namespace std;

extern Logger logFile;

namespace poisson {

   // ***** INITIALIZE STATIC VARIABLES ***** //   
   int Poisson::RHOQ_TOT = CellParams::RHOQ_TOT;
   int Poisson::PHI = CellParams::PHI;
   ObjectFactory<PoissonSolver> Poisson::solvers;
   PoissonSolver* Poisson::solver = NULL;
   bool Poisson::clearPotential = true;
   bool Poisson::is2D = false;
   string Poisson::solverName;
   Real Poisson::maxAbsoluteError = 1e-4;
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

   bool PoissonSolver::calculateBackgroundField(
            dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const std::vector<CellID>& cells) {

      phiprof::start("Background Field");
      
      if (Poisson::clearPotential == true || Parameters::tstep == 0 || Parameters::meshRepartitioned == true) {
         #pragma omp parallel for
         for (size_t c=0; c<cells.size(); ++c) {
            spatial_cell::SpatialCell* cell = mpiGrid[cells[c]];
            cell->parameters[CellParams::PHI] = 0;
            cell->parameters[CellParams::PHI_TMP] = 0;
            cell->parameters[CellParams::EXVOL] = cell->parameters[CellParams::BGEXVOL];
            cell->parameters[CellParams::EYVOL] = cell->parameters[CellParams::BGEYVOL];
            cell->parameters[CellParams::EZVOL] = cell->parameters[CellParams::BGEZVOL];
         }
      } else {
         #pragma omp parallel for
         for (size_t c=0; c<cells.size(); ++c) {
            spatial_cell::SpatialCell* cell = mpiGrid[cells[c]];
            cell->parameters[CellParams::EXVOL] = cell->parameters[CellParams::BGEXVOL];
            cell->parameters[CellParams::EYVOL] = cell->parameters[CellParams::BGEYVOL];
            cell->parameters[CellParams::EZVOL] = cell->parameters[CellParams::BGEZVOL];
         }
      }
      phiprof::stop("Background Field",cells.size(),"Spatial Cells");
      return true;
   }
   
   /** Calculate total charge density on given spatial cells.
    * @param mpiGrid Parallel grid library.
    * @param cells List of spatial cells.
    * @return If true, charge densities were successfully calculated.*/
   bool PoissonSolver::calculateChargeDensity(spatial_cell::SpatialCell* cell) {
      phiprof::start("Charge Density");
      bool success = true;

      Real rho_q = 0.0;
      #pragma omp parallel reduction (+:rho_q)
      {
         // Iterate all particle species
         for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {            
            Real rho_q_spec=0;
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
            if (blockContainer.size() == 0) continue;

            const Real charge       = getObjectWrapper().particleSpecies[popID].charge;
            const Realf* data       = blockContainer.getData();
            const Real* blockParams = blockContainer.getParameters();

            // Sum charge density over all phase-space cells
            #pragma omp for
            for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
               Real sum = 0.0;
               for (int i=0; i<WID3; ++i) sum += data[blockLID*WID3+i];

               const Real DV3 
                  = blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVX]
                  * blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVY]
                  * blockParams[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::DVZ];
               rho_q_spec += sum*DV3;
            }
            
            #ifdef DEBUG_POISSON
            bool ok = true;
            if (rho_q_spec != rho_q_spec) ok = false;
            if (rho_q != rho_q) ok = false;
            if (charge != charge) ok = false;
            if (ok == false) {
               stringstream ss;
               ss << "(POISSON SOLVER) NAN detected, rho_q: " << rho_q_spec << '\t' << rho_q << '\t';
               ss << "pop " << popID << " charge: " << charge << '\t';
               ss << endl;
               cerr << ss.str();
               exit(1);
            }
            #endif
            
            rho_q += charge*rho_q_spec;
         } // for-loop over particle species
      }
      
#warning This works for esail only
      Real t_max = 5e-3;
      Real factor = std::max((Real)0.0,(Real)1.0 + (Parameters::t - t_max));
      factor = std::min((Real)1.0,factor);
      
      cell->parameters[CellParams::RHOQ_TOT] = factor*cell->parameters[CellParams::RHOQ_EXT] + rho_q/physicalconstants::EPS_0;

      #ifdef DEBUG_POISSON
      bool ok = true;
      if (rho_q != rho_q) ok = false;
      if (ok == false) {
         stringstream ss;
         ss << "(POISSON SOLVER) NAN detected, rho_q " << rho_q << '\t';
         ss << endl;
         cerr << ss.str(); exit(1);
      }
      #endif
      
      size_t phaseSpaceCells=0;
      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
         phaseSpaceCells += cell->get_velocity_blocks(popID).size()*WID3;

      phiprof::stop("Charge Density",phaseSpaceCells,"Phase-space cells");
      return success;
   }

   /*bool PoissonSolver::calculateElectrostaticField2D(const std::vector<poisson::CellCache3D>& cells) {
      phiprof::start("Electrostatic E");
      
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const Real DX = 2*cells[c].parameters[0][CellParams::DX];
         const Real DY = 2*cells[c].parameters[0][CellParams::DY];

         const Real phi_01 = cells[c].parameters[1][CellParams::PHI]; // -x neighbor
         const Real phi_21 = cells[c].parameters[2][CellParams::PHI]; // +x neighbor
         const Real phi_10 = cells[c].parameters[3][CellParams::PHI]; // -y neighbor
         const Real phi_12 = cells[c].parameters[4][CellParams::PHI]; // +y neighbor

         cells[c].parameters[0][CellParams::EXVOL] -= (phi_21-phi_01)/DX;
         cells[c].parameters[0][CellParams::EYVOL] -= (phi_12-phi_10)/DY;
      }

      phiprof::stop("Electrostatic E",cells.size(),"Spatial Cells");
      return true;
   }*/
   
   /*bool PoissonSolver::calculateElectrostaticField3D(const std::vector<poisson::CellCache3D>& cells) {
      phiprof::start("Electrostatic E");

      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         const Real DX = 2*cells[c].parameters[0][CellParams::DX];
         const Real DY = 2*cells[c].parameters[0][CellParams::DY];
         const Real DZ = 2*cells[c].parameters[0][CellParams::DZ];

         const Real phi_011 = cells[c].parameters[1][CellParams::PHI]; // -x neighbor
         const Real phi_211 = cells[c].parameters[2][CellParams::PHI]; // +x neighbor
         const Real phi_101 = cells[c].parameters[3][CellParams::PHI]; // -y neighbor
         const Real phi_121 = cells[c].parameters[4][CellParams::PHI]; // +y neighbor
         const Real phi_110 = cells[c].parameters[5][CellParams::PHI]; // -z neighbor
         const Real phi_112 = cells[c].parameters[6][CellParams::PHI]; // +z neighbor

         cells[c].parameters[0][CellParams::EXVOL] -= (phi_211-phi_011)/DX;
         cells[c].parameters[0][CellParams::EYVOL] -= (phi_121-phi_101)/DY;
         cells[c].parameters[0][CellParams::EZVOL] -= (phi_112-phi_110)/DZ;
      }

      phiprof::stop("Electrostatic E",cells.size(),"Spatial Cells");
      return true;
   }*/

   /*bool PoissonSolver::checkGaussLaw(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     const std::vector<poisson::CellCache3D>& cells,
                                     Real& efieldFlux,Real& totalCharge) {
      bool success = true;
      Real chargeSum = 0;
      Real eFluxSum  = 0;
      
      #pragma omp parallel for reduction(+:chargeSum,eFluxSum)
      for (size_t c=0; c<cells.size(); ++c) {
         if (cells[c].cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;
         
         const Real D3 
            = cells[c][0][CellParams::DX]
            * cells[c][0][CellParams::DY]
            * cells[c][0][CellParams::DZ];
         chargeSum += cells[c][0][CellParams::RHOQ_TOT]*D3;

         spatial_cell::SpatialCell* nbr;
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cells[c].cellID);
         
         // -x neighbor
         indices[0] -= 1;
         nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
         if (nbr != NULL) if (nbr->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            const Real area   = cells[c][0][CellParams::DY]*cells[c][0][CellParams::DZ];
            const Real Ex     = cells[c][0][CellParams::EXVOL] - cells[c][0][CellParams::BGEXVOL];
            const Real Ex_nbr = nbr->parameters[CellParams::EXVOL] - nbr->parameters[CellParams::BGEXVOL];

            eFluxSum -= 0.5*(Ex+Ex_nbr)*area;
         }
         // +x neighbor
         indices[0] += 2;
         nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
         if (nbr != NULL) if (nbr->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            const Real area = cells[c][0][CellParams::DY]*cells[c][0][CellParams::DZ];
            const Real Ex   = cells[c][0][CellParams::EXVOL] - cells[c][0][CellParams::BGEXVOL];
            const Real Ex_nbr = nbr->parameters[CellParams::EXVOL] - nbr->parameters[CellParams::BGEXVOL];

            eFluxSum += 0.5*(Ex+Ex_nbr)*area;
         }
         indices[0] -= 1;

         // -y neighbor
         indices[1] -= 1;
         nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
         if (nbr != NULL) if (nbr->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            const Real area = cells[c][0][CellParams::DX]*cells[c][0][CellParams::DZ];
            const Real Ey   = cells[c][0][CellParams::EYVOL] - cells[c][0][CellParams::BGEYVOL];
            const Real Ey_nbr = nbr->parameters[CellParams::EYVOL] - nbr->parameters[CellParams::BGEYVOL];

            eFluxSum -= 0.5*(Ey+Ey_nbr)*area;
         }
         // +y neighbor
         indices[1] += 2;
         nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
         if (nbr != NULL) if (nbr->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            const Real area = cells[c][0][CellParams::DX]*cells[c][0][CellParams::DZ];
            const Real Ey   = cells[c][0][CellParams::EYVOL] - cells[c][0][CellParams::BGEYVOL];
            const Real Ey_nbr = nbr->parameters[CellParams::EYVOL] - nbr->parameters[CellParams::BGEYVOL];

            eFluxSum += 0.5*(Ey+Ey_nbr)*area;
         }
         indices[1] -= 1;
         
         // -z neighbor
         indices[2] -= 1;
         nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
         if (nbr != NULL) if (nbr->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            const Real area = cells[c][0][CellParams::DX]*cells[c][0][CellParams::DY];
            const Real Ez   = cells[c][0][CellParams::EZVOL] - cells[c][0][CellParams::BGEZVOL];
            eFluxSum -= Ez*area;
         }
         // +z neighbor
         indices[2] += 2;
         nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
         if (nbr != NULL) if (nbr->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            const Real area = cells[c][0][CellParams::DX]*cells[c][0][CellParams::DY];
            const Real Ez   = cells[c][0][CellParams::EZVOL] - cells[c][0][CellParams::BGEZVOL];
            eFluxSum += Ez*area;
         }
         indices[2] -= 1;
      }

      efieldFlux  += eFluxSum;
      totalCharge += chargeSum;

      return success;
   }*/

   /** Estimate the error in the numerical solution of the electrostatic potential.
    * The error is calculated as the difference of new and old potential, divided 
    * by the new potential value. This function works for both 2D and 3D solver.
    * This function must be called simultaneously by all MPI processes.
    * @param mpiGrid Parallel grid.
    * @return The relative error in Poisson equation solution. The return value is 
    * the same at all MPI processes.*/
   /*Real PoissonSolver::error(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
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
      phiprof::stop("Potential Change",Poisson::localCellParams.size(),"Spatial Cells");

      // Reduce max error to all MPI processes
      phiprof::start("MPI");
      Real globalMaxError;
      MPI_Allreduce(&maxError,&globalMaxError,1,MPI_Type<Real>(),MPI_MAX,MPI_COMM_WORLD);
      delete [] maxError; maxError = NULL;
      phiprof::stop("MPI");

      return globalMaxError;
   }*/
   
   Real PoissonSolver::maxError2D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::start("Evaluate Error");

      // DEBUG: Make sure values are up to date
      SpatialCell::set_mpi_transfer_type(spatial_cell::Transfer::CELL_RHOQ_TOT,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(spatial_cell::Transfer::CELL_PHI,false);
      mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);

      //Real localError = 0;
      Real* maxError = new Real[omp_get_max_threads()];
      const vector<CellID>& cells = getLocalCells();

      #pragma omp parallel 
        {
           const int tid = omp_get_thread_num();
           maxError[tid] = 0;

           #pragma omp for
           for (size_t c=0; c<cells.size(); ++c) {
              CellID cellID = cells[c];

              // Skip cells on domain boundaries:
              if (mpiGrid[cellID]->sysBoundaryFlag != 1) {
                 mpiGrid[cellID]->parameters[CellParams::PHI_TMP] = 0;
                 continue;
              }

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

              // Evaluate error
              Real DX2 = mpiGrid[cellID]->parameters[CellParams::DX]*mpiGrid[cellID]->parameters[CellParams::DX];

              Real RHS = phi_011+phi_211+phi_101+phi_121-4*phi_111;
              Real cellError = fabs(-rho_q*DX2 - RHS);
              mpiGrid[cellID]->parameters[CellParams::PHI_TMP] = cellError;
              if (fabs(cellError) > maxError[tid]) maxError[tid] = fabs(cellError);
           } // for-loop over cells

        } // #pragma omp parallel

      // Compute max error (over threads)
      for (int i=1; i<omp_get_max_threads(); ++i) {
         if (maxError[i] > maxError[0]) maxError[0] = maxError[i];
      }

      Real globalMaxError;
      MPI_Allreduce(maxError,&globalMaxError,1,MPI_Type<Real>(),MPI_MAX,MPI_COMM_WORLD);      

      delete [] maxError; maxError = NULL;

      phiprof::stop("Evaluate Error",cells.size(),"Spatial Cells");

      return globalMaxError;
   }

   // ***** DEFINITIONS OF HIGH-LEVEL DRIVER FUNCTIONS ***** //

   bool initialize(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;

      Poisson::solvers.add("Jacobi",makeJacobi);
      Poisson::solvers.add("SOR",makeSOR);
      Poisson::solvers.add("CG",makeCG);

      // Create and initialize the Poisson solver
      Poisson::solver = Poisson::solvers.create(Poisson::solverName);
      if (Poisson::solver == NULL) {
         logFile << "(POISSON SOLVER) ERROR: Failed to create Poisson solver '" << Poisson::solverName << "'" << endl << write;
         return false;
      } else {
         if (Poisson::solver->initialize() == false) success = false;
         if (success == true) 
           logFile << "(POISSON SOLVER) Successfully initialized Poisson solver '" << Poisson::solverName << "'" << endl << write;
         else {
            logFile << "(POISSON SOLVER) ERROR: Failed to initialize Poisson solver '" << Poisson::solverName << "'" << endl << write;
            return success;
         }
      }

      // Set up the initial state unless the simulation was restarted
      //if (Parameters::isRestart == true) return success;

      for (size_t c=0; c<getLocalCells().size(); ++c) {
         spatial_cell::SpatialCell* cell = mpiGrid[getLocalCells()[c]];
         if (Poisson::solver->calculateChargeDensity(cell) == false) {
            logFile << "(POISSON SOLVER) ERROR: Failed to calculate charge density in " << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }
      }

      // Force calculateBackgroundField to reset potential arrays to zero values.
      // This may not otherwise happen if the simulation was restarted.
      const bool oldValue = Poisson::clearPotential;
      Poisson::clearPotential = true;
      if (Poisson::solver->calculateBackgroundField(mpiGrid,getLocalCells()) == false) {
         logFile << "(POISSON SOLVER) ERROR: Failed to calculate background field in " << __FILE__ << ":" << __LINE__ << endl << write;
         success = false;
      }
      if (solve(mpiGrid) == false) {
         logFile << "(POISSON SOLVER) ERROR: Failed to solve potential in " << __FILE__ << ":" << __LINE__ << endl << write;
         success = false;
      }
      Poisson::clearPotential = oldValue;

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
      bool success = true;
      
      // If mesh partitioning has changed, recalculate spatial 
      // cell parameters pointer cache:
      if (Parameters::meshRepartitioned == true) {
         phiprof::start("Cache Cell Parameters");
         Poisson::cacheCellParameters(mpiGrid,getLocalCells());
         phiprof::stop("Cache Cell Parameters");
      }

      // Solve Poisson equation
      if (success == true) if (Poisson::solver != NULL) {
         if (Poisson::solver->calculateBackgroundField(mpiGrid,getLocalCells()) == false) success = false;

         SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
         mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);

         if (Poisson::solver->solve(mpiGrid) == false) success = false;
      }

      // Add electrostatic electric field to volume-averaged E
      //if (success == true) if (Poisson::solver != NULL) {
      //   if (Poisson::solver->calculateElectrostaticField(mpiGrid) == false) success = false;
      //}
      
      phiprof::stop("Poisson Solver (Total)");
      return success;
   }

} // namespace poisson
