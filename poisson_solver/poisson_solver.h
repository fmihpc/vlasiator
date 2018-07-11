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
 * File:   poisson_solver.h
 * Author: sandroos
 *
 * Created on January 14, 2015, 1:42 PM
 */

#ifndef POISSON_SOLVER_H
#define	POISSON_SOLVER_H

#ifdef _OPENMP
   #include <omp.h>
#endif
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../mpiconversion.h"
#include "../object_factory.h"
#include "../spatial_cell.hpp"

namespace poisson {

   template<unsigned int VARS>
   struct CellCache2D {
      spatial_cell::SpatialCell* cell;
      Real* parameters[5];
      Real*& operator[](const int& i) {return parameters[i];}
      Real variables[VARS];
   };

   template<unsigned int VARS>
   struct CellCache3D {
#warning TEMP remove
       CellID cellID;
       spatial_cell::SpatialCell* cell;
       Real* parameters[7];
       Real*& operator[](const int& i) {return parameters[i];}
       Real* const& operator[](const int& i) const {return parameters[i];}
       Real variables[VARS];
    };

    class PoissonSolver {
    public:
       PoissonSolver();
       virtual ~PoissonSolver();

       // ***** DECLARATIONS OF VIRTUAL MEMBER FUNCTIONS ***** //

       virtual bool calculateBackgroundField(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                             const std::vector<CellID>& cells);
       virtual bool calculateChargeDensity(spatial_cell::SpatialCell* cell);
       virtual void calculateChargeDensitySingle(spatial_cell::SpatialCell* cell);
       template<unsigned int VARS>
       bool calculateElectrostaticField2D(const std::vector<poisson::CellCache3D<VARS> >& cells);
       template<unsigned int VARS>
       bool calculateElectrostaticField3D(const std::vector<poisson::CellCache3D<VARS> >& cells);
       /*template<unsigned int VARS>
       virtual bool checkGaussLaw(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                  const std::vector<poisson::CellCache3D<VARS> >& cells,
                                  Real& efieldFlux,Real& totalCharge);*/
       template<unsigned int VARS>
       Real error(std::vector<poisson::CellCache3D<VARS> >& cells);
       virtual Real maxError2D(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
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
      static Real maxAbsoluteError;                /**< Potential iteration is stopped if maximum potential 
                                                    * error drops below this value.*/
      static uint maxIterations;                   /**< Maximum number of iterations allowed, only 
                                                    * has an effect on iterative solvers.*/
      static Real minRelativePotentialChange;      /**< Iterative solvers keep on iterating the solution 
                                                    * until the change in potential during successive 
                                                    * iterations is less than this value.*/
      static std::vector<Real*> localCellParams;   /**< Pointers to spatial cell parameters, order 
						    * is the same as in getLocalCells() vector.*/
      static bool timeDependentBackground;         /**< If true, the background field / charge density is 
                                                    * time-dependent and must be recalculated each time step.*/
      
      static void cacheCellParameters(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
				      const std::vector<CellID>& cells);
   };

   bool initialize(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
   bool finalize();
   bool solve(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

   
   
   
   template<unsigned int VARS> inline
   bool PoissonSolver::calculateElectrostaticField2D(const std::vector<poisson::CellCache3D<VARS> >& cells) {
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
   }
   
   template<unsigned int VARS> inline
   bool PoissonSolver::calculateElectrostaticField3D(const std::vector<poisson::CellCache3D<VARS> >& cells) {
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
   }

   /** Estimate the error in the numerical solution of the electrostatic potential.
    * The error is calculated as the maximum deviation of nabla^2 (phi) + rho_q/epsilon_0 
    * from zero value. This function only works for a 2D (xy) solver.
    * This function must be called simultaneously by all MPI processes.
    * @param mpiGrid Parallel grid.
    * @return The error in Poisson equation solution. The return value is 
    * the same at all MPI processes.*/
   template<unsigned int VARS> inline
   Real PoissonSolver::error(std::vector<poisson::CellCache3D<VARS> >& cells) {
      phiprof::start("error evaluation");
      Real maxError = -std::numeric_limits<Real>::max();
      Real* threadMaxError = new Real[omp_get_max_threads()];
      
      #pragma omp parallel
      {
         const int tid = omp_get_thread_num();
         Real myError = 0;

         #pragma omp for
         for (size_t c=0; c<cells.size(); ++c) {
            if (cells[c].cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
               cells[c].parameters[0][CellParams::PHI_TMP] = 0;
               continue;
            }
            
            const Real DX2 = cells[c].parameters[0][CellParams::DX]*cells[c].parameters[0][CellParams::DX];
            const Real rho_q   = cells[c].parameters[0][CellParams::RHOQ_TOT];
            const Real phi_111 = cells[c].parameters[0][CellParams::PHI];
            const Real phi_011 = cells[c].parameters[1][CellParams::PHI]; // -x neighbor
            const Real phi_211 = cells[c].parameters[2][CellParams::PHI]; // +x neighbor
            const Real phi_101 = cells[c].parameters[3][CellParams::PHI]; // -y neighbor
            const Real phi_121 = cells[c].parameters[4][CellParams::PHI]; // +y neighbor
         
            const Real RHS = phi_011+phi_211+phi_101+phi_121-4*phi_111;
            const Real cellError = fabs(-rho_q*DX2-RHS);
            cells[c].parameters[0][CellParams::PHI_TMP] = cellError;
            if (cellError > myError) myError = cellError;
         }
         threadMaxError[tid] = myError;
      }

      // Calculate the maximum error over all per-thread values to variable maxError
      for (int i=1; i<omp_get_max_threads(); ++i) {
         if (threadMaxError[i] > threadMaxError[0]) threadMaxError[0] = threadMaxError[i];
      }
      maxError = threadMaxError[0];
      delete [] threadMaxError; threadMaxError = NULL;
      
      // Reduce the maximum error to all processes
      Real globalMaxError;
      MPI_Allreduce(&maxError,&globalMaxError,1,MPI_Type<Real>(),MPI_MAX,MPI_COMM_WORLD);
      phiprof::stop("error evaluation",cells.size(),"Spatial Cells");

      return globalMaxError;
   }

} // namespace poisson

#endif	// POISSON_SOLVER_H

