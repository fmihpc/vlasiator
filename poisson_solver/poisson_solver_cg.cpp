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
 * File:   poisson_solver_cg.cpp
 * Author: sandroos
 *
 * Created on January 15, 2015, 12:45 PM
 */

#include <cstdlib>
#include <iostream>
#include <omp.h>

#include "../logger.h"
#include "../grid.h"
#include "../mpiconversion.h"

#include "poisson_solver_cg.h"

#ifndef NDEBUG
   #define DEBUG_POISSON_CG
#endif

using namespace std;

extern Logger logFile;

namespace poisson {

   static std::vector<CellCache3D<cgvar::SIZE> > innerCellPointers;
   static std::vector<CellCache3D<cgvar::SIZE> > bndryCellPointers;

   PoissonSolver* makeCG() {
      return new PoissonSolverCG();
   }

   PoissonSolverCG::PoissonSolverCG(): PoissonSolver() { 

   }

   PoissonSolverCG::~PoissonSolverCG() { }

   bool PoissonSolverCG::initialize() {
      bool success = true;
      bndryCellParams[CellParams::PHI] = 0;
      bndryCellParams[CellParams::PHI_TMP] = 0;
      return success;
   }

   bool PoissonSolverCG::finalize() {
      bool success = true;
      return success;
   }
   
   inline void PoissonSolverCG::calculateAlpha(CellCache3D<cgvar::SIZE>& cell,Real& mySum0,Real& mySum1) {
      // Calculate r(transpose) * r
      mySum0 += cell.variables[cgvar::R]*cell.variables[cgvar::R];

      // Calculate p(transpose) * A * p
      Real A_p = -4*cell.parameters[0][CellParams::PHI_TMP]
        + cell.parameters[1][CellParams::PHI_TMP]
        + cell.parameters[2][CellParams::PHI_TMP]
        + cell.parameters[3][CellParams::PHI_TMP]
        + cell.parameters[4][CellParams::PHI_TMP];
      cell.variables[cgvar::A_TIMES_P] = A_p;
      mySum1 += cell.parameters[0][CellParams::PHI_TMP]*A_p;
   }
   
   /** Calculate the value of alpha parameter.
    * @return If true, all processes have the same value of alpha in PoissonSolverCG::alphaGlobal.*/
   bool PoissonSolverCG::calculateAlpha() {
      phiprof::start("calculate alpha");
      Real t_start = 0;
      if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();
      
      Real sums[2];
      Real mySum0 = 0;
      Real mySum1 = 0;

      // Calculate P(transpose)*A*P and R(transpose)*R for all local cells:
      const size_t offset = bndryCellPointers.size();
      const size_t N_cells = bndryCellPointers.size()+innerCellPointers.size();

      #pragma omp parallel for reduction(+:mySum0,mySum1)
      for (size_t c=0; c<N_cells; ++c) {
         if (c >= offset) calculateAlpha(innerCellPointers[c-offset],mySum0,mySum1);
         else             calculateAlpha(bndryCellPointers[c       ],mySum0,mySum1);
      }

      sums[0] = mySum0;
      sums[1] = mySum1;
      
      // Measure computation time (if needed)
      if (Parameters::prepareForRebalance == true) {
         const size_t N = max((size_t)1,innerCellPointers.size()+bndryCellPointers.size());
         Real t_average = (MPI_Wtime() - t_start) / N;

         #pragma omp parallel
         {
            #pragma omp for nowait
            for (size_t c=0; c<bndryCellPointers.size(); ++c) {
               CellCache3D<cgvar::SIZE>& cell = bndryCellPointers[c];
               cell.parameters[0][CellParams::LBWEIGHTCOUNTER] += t_average;
            }
            #pragma omp for nowait
            for (size_t c=0; c<innerCellPointers.size(); ++c) {
               CellCache3D<cgvar::SIZE>& cell = innerCellPointers[c];
               cell.parameters[0][CellParams::LBWEIGHTCOUNTER] += t_average;
            }            
         }
      }

      phiprof::stop("calculate alpha",bndryCellPointers.size()+innerCellPointers.size(),"Spatial Cells");

      // Reduce sums to master process:
      #warning TEST if Allreduce is faster here
      phiprof::start("MPI (Alpha)");
      MPI_Reduce(sums,&(globalVariables[cgglobal::R_T_R]),2,MPI_Type<Real>(),MPI_SUM,0,MPI_COMM_WORLD);

      // Calculate alpha and broadcast to all processes:
      globalVariables[cgglobal::ALPHA]
              = globalVariables[cgglobal::R_T_R] 
              /(globalVariables[cgglobal::P_T_A_P] + 100*numeric_limits<Real>::min());
      MPI_Bcast(globalVariables,3,MPI_Type<Real>(),0,MPI_COMM_WORLD);
      phiprof::stop("MPI (Alpha)");
      return true;
   }

   bool PoissonSolverCG::calculateElectrostaticField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
      
      mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      
      // Calculate electric field on inner cells
      if (Poisson::is2D == true) {
         if (calculateElectrostaticField2D(innerCellPointers) == false) success = false;
      } else {
         if (calculateElectrostaticField3D(innerCellPointers) == false) success = false;
      }

      mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      
      // Calculate electric field on boundary cells
      if (Poisson::is2D == true) {
         if (calculateElectrostaticField2D(bndryCellPointers) == false) success = false;
      } else {
         if (calculateElectrostaticField3D(bndryCellPointers) == false) success = false;
      }

      return success;
   }

   void PoissonSolverCG::cachePointers2D(
               dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
               const std::vector<CellID>& cells,
               std::vector<poisson::CellCache3D<cgvar::SIZE> >& cellCache) {
      cellCache.clear();

      for (size_t c=0; c<cells.size(); ++c) {
         // DO_NOT_COMPUTE cells are skipped
         if (mpiGrid[cells[c]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         //if (mpiGrid[cells[c]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;

         // Calculate cell i/j/k indices
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cells[c]);

         CellCache3D<cgvar::SIZE> cache;
         cache.cellID = cells[c];
         cache.cell = mpiGrid[cells[c]];
         cache[0]   = mpiGrid[cells[c]]->parameters;

         #ifdef DEBUG_POISSON_CG
         if (cache.cell == NULL) {
            stringstream s;
            s << "ERROR, NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
            s << "\t Cell ID " << cells[c] << endl;
            cerr << s.str();
            exit(1);
         }
         #endif

         spatial_cell::SpatialCell* dummy = NULL;
         switch (mpiGrid[cells[c]]->sysBoundaryFlag) {
            case sysboundarytype::DO_NOT_COMPUTE:
               break;
            case sysboundarytype::NOT_SYSBOUNDARY:
               // Fetch pointers to this cell's (cell) parameters array, 
               // and pointers to +/- xyz face neighbors' arrays               
               indices[0] -= 1; cache[1] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
               indices[0] += 2; cache[2] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
               indices[0] -= 1;
            
               if (indices[1] == 2) {
                  indices[1] -= 1;
                  dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
                  if (dummy->sysBoundaryFlag == sysboundarytype::ANTISYMMETRIC) {
                     cache[3] = cache[0];
                  } else {
                     cache[3] = dummy->get_cell_parameters();
                  }
               } else {
                  indices[1] -= 1; cache[3] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
               }

               indices[1] += 2; cache[4] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
               indices[1] -= 1;
               break;

            case sysboundarytype::ANTISYMMETRIC:
               continue;
               
               // Get +/- x-neighbor pointers
               indices[0] -= 1;
               dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
               if (dummy == NULL) cache[1] = bndryCellParams;
               else               cache[1] = dummy->parameters;
               indices[0] += 2;
               dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
               if (dummy == NULL) cache[2] = bndryCellParams;
               else               cache[2] = dummy->parameters;
               indices[0] -= 1;

               // Set +/- y-neighbors both point to +y neighbor 
               // if we are at the lower y-boundary, otherwise set both 
               // y-neighbors point to -y neighbor.
               if (indices[1] == 1) {
                  cache[3] = cache[0];
                  indices[1] += 1;
                  dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
                  if (dummy == NULL) {
                     cache[4] = bndryCellParams;
                  } else {
                     cache[4] = dummy->parameters;
                  }
                  indices[1] -= 1;
               } else {
                  cache[4] = cache[0];
                  indices[1] -= 1;
                  dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
                  if (dummy == NULL) {
                     cache[3] = bndryCellParams;
                  } else {
                     cache[3] = dummy->parameters;
                  }
                  indices[1] += 1;
               }
               break;

            default:
               indices[0] -= 1;
               dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
               //if (dummy == NULL) cache[1] = bndryCellParams;
               if (dummy == NULL) continue;
               else               cache[1] = dummy->parameters;
               indices[0] += 2;
               dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
               //if (dummy == NULL) cache[2] = bndryCellParams;
               if (dummy == NULL) continue;
               else               cache[2] = dummy->parameters;
               indices[0] -= 1;

               indices[1] -= 1; 
               dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
               //if (dummy == NULL) cache[3] = bndryCellParams;
               if (dummy == NULL) continue;
               else               cache[3] = dummy->parameters;
               indices[1] += 2;
               dummy = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
               //if (dummy == NULL) cache[4] = bndryCellParams;
               if (dummy == NULL) continue;
               else               cache[4] = dummy->parameters;
               indices[1] -= 1;
               break;
         }

         cellCache.push_back(cache);
      } // for-loop over spatial cells
   }
/*
   void PoissonSolverCG::cachePointers3D(
            dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const std::vector<CellID>& cells,std::vector<poisson::CellCache3D>& redCache,
            std::vector<poisson::CellCache3D>& blackCache) {
      redCache.clear();
      blackCache.clear();
      
      for (size_t c=0; c<cells.size(); ++c) {
         // Calculate cell i/j/k indices
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cells[c]);

         if ((indices[0] + indices[1]%2 + indices[2]%2) % 2 == RED) {
            CellCache3D cache;

            // Cells on domain boundaries are not iterated
            if (mpiGrid[cells[c]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            
            // Fetch pointers to this cell's (cell) parameters array, 
            // and pointers to +/- xyz face neighbors' arrays
            cache.cell = mpiGrid[cells[c]];
            cache[0] = mpiGrid[cells[c]]->parameters;
            
            indices[0] -= 1; cache[1] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] += 2; cache[2] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] -= 1;
            
            indices[1] -= 1; cache[3] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] += 2; cache[4] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] -= 1;
            
            indices[2] -= 1; cache[5] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] += 2; cache[6] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] -= 1;
            
            redCache.push_back(cache);
         } else {
            CellCache3D cache;
            
            // Cells on domain boundaries are not iterated
            if (mpiGrid[cells[c]]->sysBoundaryFlag != 1) continue;
            
            // Fetch pointers to this cell's (cell) parameters array,
            // and pointers to +/- xyz face neighbors' arrays
            cache.cell = mpiGrid[cells[c]];
            cache[0] = mpiGrid[cells[c]]->parameters;
            
            indices[0] -= 1; cache[1] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] += 2; cache[2] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[0] -= 1;
            
            indices[1] -= 1; cache[3] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] += 2; cache[4] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[1] -= 1;
            
            indices[2] -= 1; cache[5] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] += 2; cache[6] = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ]->parameters;
            indices[2] -= 1;
            
            blackCache.push_back(cache);
         }
      }
   }
 */

   void PoissonSolverCG::bvalue(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                std::vector<CellCache3D<cgvar::SIZE> >& cells) {
      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cells[c].cellID);
         CellCache3D<cgvar::SIZE>& cell = cells[c];
         
         if (indices[0] < 2 && (indices[1] > 1 && indices[1] < Parameters::ycells_ini-2)) {
            Real RHS = -3*cell.parameters[0][CellParams::PHI]
                       + cell.parameters[3][CellParams::PHI] + cell.parameters[4][CellParams::PHI];
            cell.variables[cgvar::R] = cell.variables[cgvar::B] - RHS;
            continue;
         }

         if ((indices[1] < 2 || indices[1] > Parameters::ycells_ini-3) 
          && (indices[0] > 1 && indices[0] < Parameters::xcells_ini-2)) {
            Real RHS = -3*cell.parameters[0][CellParams::PHI]
                       + cell.parameters[1][CellParams::PHI] + cell.parameters[2][CellParams::PHI];
            cell.variables[cgvar::R] = cell.variables[cgvar::B] - RHS;
            continue;
         }

         if (indices[0] == 1 && (indices[1] == 1 || indices[1] == Parameters::ycells_ini-2)) {
            cell.variables[cgvar::R] = 0;
            continue;
         }
      }
   }
   
   bool PoissonSolverCG::solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;

      // If mesh partitioning has changed, recalculate pointer caches
      if (Parameters::meshRepartitioned == true) {
         phiprof::start("Pointer Caching");
         if (Poisson::is2D == true) {
            cachePointers2D(mpiGrid,mpiGrid.get_local_cells_on_process_boundary(POISSON_NEIGHBORHOOD_ID),bndryCellPointers);
            cachePointers2D(mpiGrid,mpiGrid.get_local_cells_not_on_process_boundary(POISSON_NEIGHBORHOOD_ID),innerCellPointers);
         } else {
            //cachePointers3D(mpiGrid,mpiGrid.get_local_cells_on_process_boundary(POISSON_NEIGHBORHOOD_ID),bndryCellPointers);
            //cachePointers3D(mpiGrid,mpiGrid.get_local_cells_not_on_process_boundary(POISSON_NEIGHBORHOOD_ID),innerCellPointers);
         }
         phiprof::stop("Pointer Caching");
      }

      // Calculate charge density and sync
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_RHOQ_TOT,false);

      //for (size_t c=0; c<bndryCellPointers.size(); ++c) calculateChargeDensity(bndryCellPointers[c].cell);
      //mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      //for (size_t c=0; c<innerCellPointers.size(); ++c) calculateChargeDensity(innerCellPointers[c].cell);
      
      phiprof::start("Charge Density");
      #pragma omp parallel
      {
         const int tid = omp_get_thread_num();

         #pragma omp for nowait
         for (size_t c=0; c<bndryCellPointers.size(); ++c) {
            calculateChargeDensitySingle(bndryCellPointers[c].cell);
         }
         
         if (tid == 0) {
            mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
         }
         
         #pragma omp for nowait
         for (size_t c=0; c<innerCellPointers.size(); ++c) {
            calculateChargeDensitySingle(innerCellPointers[c].cell);
         }
      }
      phiprof::stop("Charge Density");

      phiprof::start("MPI (RHOQ)");
      mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      phiprof::stop("MPI (RHOQ)");

      // Initialize
      if (startIteration(mpiGrid) == false) {
         cerr << "(POISSON CG) ERROR has occurred in startIteration in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }

      iterations = 0;
      Real relPotentialChange = 0;
      do {
         const int N_iterations = 1;

         SpatialCell::set_mpi_transfer_type(Transfer::CELL_RHOQ_TOT,false);
         mpiGrid.update_copies_of_remote_neighbors(POISSON_NEIGHBORHOOD_ID);
         //bvalue(mpiGrid,innerCellPointers);
         //bvalue(mpiGrid,bndryCellPointers);

         if (calculateAlpha() == false) {
            logFile << "(POISSON SOLVER CG) ERROR: Failed to calculate 'alpha' in ";
            logFile << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }
         if (update_x_r() == false) {
            logFile << "(POISSON SOLVER CG) ERROR: Failed to update x and r vectors in ";
            logFile << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }
         if (update_p(mpiGrid) == false) {
            logFile << "(POISSON SOLVER CG) ERROR: Failed to update p vector in ";
            logFile << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }
         
         iterations += N_iterations;
         
         //cerr << iterations << '\t' << globalVariables[cgglobal::R_MAX] << endl;
         //if (mpiGrid.get_rank() == 0) {
         //   cerr << iterations << "\t" << globalVariables[cgglobal::R_MAX] << endl;
         //}
         
         if (iterations >= Poisson::maxIterations) break;
         if (globalVariables[cgglobal::R_MAX] < Poisson::maxAbsoluteError) break;
      } while (true);      

      if (calculateElectrostaticField(mpiGrid) == false) {
         logFile << "(POISSON SOLVER CG) ERROR: Failed to calculate electrostatic field in ";
         logFile << __FILE__ << ":" << __LINE__ << endl << write;
         success = false;
      }
      
      error<cgvar::SIZE>(innerCellPointers);
      error<cgvar::SIZE>(bndryCellPointers);
      
      return success;
   }

   /**
    * Upon successful return, CellParams::PHI_TMP has the correct value of P0 
    * on all cells (local and buffered).
    * @param mpiGrid Parallel grid library.
    * @return If true, CG solver is ready to iterate.*/
   bool PoissonSolverCG::startIteration(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      bool success = true;
      if (startIteration(bndryCellPointers) == false) success = false;
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
      mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);

      if (startIteration(innerCellPointers) == false) success = false;
      mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);

      return success;
   }

   /** Calculate the value of alpha for this iteration.
    * @return If true, CG solver is ready to iterate.*/
   bool PoissonSolverCG::startIteration(std::vector<CellCache3D<cgvar::SIZE> >& cells) {
      phiprof::start("start iteration");

      #pragma omp parallel for
      for (size_t c=0; c<cells.size(); ++c) {
         // Convert charge density to charge density times cell size squared / epsilon0:
         const Real DX2 = cells[c].parameters[0][CellParams::DX]*cells[c].parameters[0][CellParams::DX];
         cells[c].variables[cgvar::B] = -cells[c].parameters[0][CellParams::RHOQ_TOT]*DX2;

         // Calculate R0:
         Real RHS = -4*cells[c].parameters[0][CellParams::PHI]
                  + cells[c].parameters[1][CellParams::PHI] + cells[c].parameters[2][CellParams::PHI]
                  + cells[c].parameters[3][CellParams::PHI] + cells[c].parameters[4][CellParams::PHI];
         cells[c].variables[cgvar::R] = cells[c].variables[cgvar::B] - RHS;

         // Calculate P0:
         const Real P0 = cells[c].variables[cgvar::R];
         cells[c].parameters[0][CellParams::PHI_TMP] = P0;
      }
      phiprof::stop("start iteration",cells.size(),"Spatial Cells");

      return true;
   }
   
   bool PoissonSolverCG::update_p(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::start("R transpose * R");
      
      Real R_T_R = 0;
      Real t_start = 0;
      Real t_total = 0;
      if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();

      #pragma omp parallel reduction(+:R_T_R)
      {
         #pragma omp for
         for (size_t c=0; c<bndryCellPointers.size(); ++c) {
            CellCache3D<cgvar::SIZE>& cell = bndryCellPointers[c];
            R_T_R += cell.variables[cgvar::R]*cell.variables[cgvar::R];
         }

         #pragma omp for
         for (size_t c=0; c<innerCellPointers.size(); ++c) {
            CellCache3D<cgvar::SIZE>& cell = innerCellPointers[c];
            R_T_R += cell.variables[cgvar::R]*cell.variables[cgvar::R];
         }
      }
      const size_t N_cells = bndryCellPointers.size()+innerCellPointers.size();

      if (Parameters::prepareForRebalance == true) t_total = (MPI_Wtime() - t_start);
      phiprof::stop("R transpose * R",N_cells,"Spatial Cells");

      // Reduce value to master, calculate beta parameter and broadcast it to all processes:
      phiprof::start("MPI (R transp R)");
      Real global_R_T_R;
      MPI_Reduce(&R_T_R,&global_R_T_R,1,MPI_Type<Real>(),MPI_SUM,0,MPI_COMM_WORLD);
      globalVariables[cgglobal::BETA] = global_R_T_R / (globalVariables[cgglobal::R_T_R] + 100*numeric_limits<Real>::min());
      globalVariables[cgglobal::R_T_R] = global_R_T_R;
      MPI_Bcast(globalVariables,cgglobal::SIZE,MPI_Type<Real>(),0,MPI_COMM_WORLD);
      phiprof::stop("MPI (R transp R)");

      phiprof::start("update P");
      t_start = 0;
      if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();

      #pragma omp parallel
      {
         const int tid = omp_get_thread_num();
         #pragma omp for nowait
         for (size_t c=0; c<bndryCellPointers.size(); ++c) {
            CellCache3D<cgvar::SIZE>& cell = bndryCellPointers[c];
            cell.parameters[0][CellParams::PHI_TMP]
              = cell.variables[cgvar::R]
              + globalVariables[cgglobal::BETA]*cell.parameters[0][CellParams::PHI_TMP];
         }

         // Send new P values to neighbor processes:
         if (tid == 0) {
            SpatialCell::set_mpi_transfer_type(Transfer::CELL_PHI,false);
            mpiGrid.start_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
         }

         #pragma omp for nowait
         for (size_t c=0; c<innerCellPointers.size(); ++c) {
            CellCache3D<cgvar::SIZE>& cell = innerCellPointers[c];
            cell.parameters[0][CellParams::PHI_TMP]
              = cell.variables[cgvar::R]
              + globalVariables[cgglobal::BETA]*cell.parameters[0][CellParams::PHI_TMP];
         }
      }
      
      // Measure computation time if needed
      if (Parameters::prepareForRebalance == true) {
         const size_t N = max((size_t)1,bndryCellPointers.size()+innerCellPointers.size());
         t_total += (MPI_Wtime() - t_start);
         const Real t_average = t_total / N;

         #pragma omp parallel
         {
            #pragma omp for nowait
            for (size_t c=0; c<bndryCellPointers.size(); ++c) {
               CellCache3D<cgvar::SIZE>& cell = bndryCellPointers[c];
               cell.parameters[0][CellParams::LBWEIGHTCOUNTER] += t_average;
            }
            #pragma omp for nowait
            for (size_t c=0; c<innerCellPointers.size(); ++c) {
               CellCache3D<cgvar::SIZE>& cell = innerCellPointers[c];
               cell.parameters[0][CellParams::LBWEIGHTCOUNTER] += t_average;
            }            
         }
      }
      phiprof::stop("update P",innerCellPointers.size()+bndryCellPointers.size(),"Spatial Cells");

      // Wait for MPI to complete:
      phiprof::start("MPI (update P)");
      mpiGrid.wait_remote_neighbor_copy_updates(POISSON_NEIGHBORHOOD_ID);
      phiprof::stop("MPI (update P)");

      return true;
   }

   /**
    * 
    * @return */
   bool PoissonSolverCG::update_x_r() {
      phiprof::start("update x and r");
      
      Real t_start = 0;
      if (Parameters::prepareForRebalance == true) t_start = MPI_Wtime();
      Real R_max = -numeric_limits<Real>::max();
      Real* thread_R_max = new Real[omp_get_max_threads()];
      
      #pragma omp parallel
      {         
         Real my_R_max = 0;
         #pragma omp for nowait
         for (size_t c=0; c<bndryCellPointers.size(); ++c) {
            CellCache3D<cgvar::SIZE>& cell = bndryCellPointers[c];
            cell.parameters[0][CellParams::PHI] += globalVariables[cgglobal::ALPHA]*cell.parameters[0][CellParams::PHI_TMP];
            cell.variables[cgvar::R] -= globalVariables[cgglobal::ALPHA]*cell.variables[cgvar::A_TIMES_P];
            if (fabs(cell.variables[cgvar::R]) > my_R_max) my_R_max = fabs(cell.variables[cgvar::R]);
         }
         #pragma omp for nowait
         for (size_t c=0; c<innerCellPointers.size(); ++c) {
            CellCache3D<cgvar::SIZE>& cell = innerCellPointers[c];
            cell.parameters[0][CellParams::PHI] += globalVariables[cgglobal::ALPHA]*cell.parameters[0][CellParams::PHI_TMP];
            cell.variables[cgvar::R] -= globalVariables[cgglobal::ALPHA]*cell.variables[cgvar::A_TIMES_P];
            if (fabs(cell.variables[cgvar::R]) > my_R_max) my_R_max = fabs(cell.variables[cgvar::R]);
         }
         const int tid = omp_get_thread_num();
         thread_R_max[tid] = my_R_max;
      }

      for (int tid=1; tid<omp_get_max_threads(); ++tid) {
         if (thread_R_max[tid] > thread_R_max[0]) thread_R_max[0] = thread_R_max[tid];
      }
      R_max = thread_R_max[0];
      delete [] thread_R_max; thread_R_max = NULL;

      size_t N_cells = bndryCellPointers.size() + innerCellPointers.size();

      // Measure computation time if needed
      if (Parameters::prepareForRebalance == true) {
         const Real t_average = (MPI_Wtime() - t_start) / max((size_t)1,N_cells);

         #pragma omp parallel
         {
            #pragma omp for nowait
            for (size_t c=0; c<bndryCellPointers.size(); ++c) {
               CellCache3D<cgvar::SIZE>& cell = bndryCellPointers[c];
               cell.parameters[0][CellParams::LBWEIGHTCOUNTER] += t_average;
            }
            #pragma omp for nowait
            for (size_t c=0; c<innerCellPointers.size(); ++c) {
               CellCache3D<cgvar::SIZE>& cell = innerCellPointers[c];
               cell.parameters[0][CellParams::LBWEIGHTCOUNTER] += t_average;
            }            
         }
      }
      
      phiprof::stop("update x and r",N_cells,"Spatial Cells");

      phiprof::start("MPI (x and r)");
      MPI_Allreduce(&R_max,&(globalVariables[cgglobal::R_MAX]),1,MPI_Type<Real>(),MPI_MAX,MPI_COMM_WORLD);
      phiprof::stop("MPI (x and r)");

      return true;
   }

} // namespace poisson
