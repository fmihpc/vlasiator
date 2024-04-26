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
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdint.h>

#ifdef _OPENMP
   #include <omp.h>
#endif
#include <zoltan.h>
#include <phiprof.hpp>
#include <dccrg.hpp>

#include "../vlasovmover.h"
#include "../spatial_cell_wrapper.hpp"
#include "../grid.h"
#include "../definitions.h"
#include "../iowrite.h"
#include "../object_wrapper.h"

#include "../vlasovsolver/cpu_moments.h"
//#include "cpu_acc_semilag.hpp"
//#include "cpu_trans_map.hpp"

using namespace std;
using namespace spatial_cell;

#warning TESTING can be removed later
static void writeVelMesh(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const vector<CellID>& cells = getLocalCells();

   static int counter=-1;
   if (counter < 0) {
      counter = Parameters::systemWrites.size();
      Parameters::systemWriteDistributionWriteStride.push_back(1);
      Parameters::systemWriteName.push_back("velocity-trans");
      Parameters::systemWriteDistributionWriteXlineStride.push_back(0);
      Parameters::systemWriteDistributionWriteYlineStride.push_back(0);
      Parameters::systemWriteDistributionWriteZlineStride.push_back(0);
      Parameters::systemWriteTimeInterval.push_back(-1.0);
      Parameters::systemWrites.push_back(0);
   }
   writeGrid(mpiGrid,NULL,counter,true);
   ++Parameters::systemWrites[counter];
}

Real calculateTotalMass(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID) {
   const vector<CellID>& local_cells = getLocalCells();
   Real sum=0.0;
   for (size_t c=0; c<local_cells.size(); ++c) {
      const CellID cellID = local_cells[c];
      SpatialCell* cell = mpiGrid[cellID];
      
      for (vmesh::LocalID blockLID=0; blockLID<cell->get_number_of_velocity_blocks(popID); ++blockLID) {
         const Real* parameters = cell->get_block_parameters(blockLID);
         const Realf* data = cell->get_data(blockLID);
         Real blockMass = 0.0;
         for (int i=0; i<WID3; ++i) {
            blockMass += data[i];
         }
         const Real DV3 = parameters[BlockParams::DVX]*parameters[BlockParams::DVY]*parameters[BlockParams::DVZ];
         sum += blockMass*DV3;
      }      
   }
   
   Real globalMass=0.0;
   MPI_Allreduce(&sum,&globalMass,1,MPI_Type<Real>(),MPI_SUM,MPI_COMM_WORLD);
   return globalMass;
}

/*!
  
  Propagates the distribution function in spatial space. 
  
  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

 * REQUIREMENTS: Remote neighbor distribution functions must've 
 * been synchronized before calling this function.
*/

void calculateSpatialTranslation(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,Real dt) {
   return;

   /*
   typedef Parameters P;
   int trans_timer;

   // DEBUG remove
   dt = 1.0;
   
   //phiprof::start("semilag-trans");
   phiprof::start("compute_cell_lists");
   const vector<CellID>& local_cells = getLocalCells();
   phiprof::stop("compute_cell_lists");

   // Note: mpiGrid.is_local( cellID ) == true if cell is local

   static int cntr=0;
   if (cntr == 0) {
      writeVelMesh(mpiGrid);
      if (mpiGrid.get_rank() == 0) {
         cout << "Initial mass is " << calculateTotalMass(mpiGrid) << endl;
      }
      cntr=1;
   }
   
   static int dim=2;

   // Generate target mesh
   phiprof::start("target mesh generation");
   for (size_t c=0; c<local_cells.size(); ++c) {
      createTargetMesh(mpiGrid,local_cells[c],dim,false);
   }
   phiprof::stop("target mesh generation");

   #warning DEBUG remove me
   for (size_t c=0; c<local_cells.size(); ++c) {
      SpatialCell* spatial_cell = mpiGrid[local_cells[c]];
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh_temporary();
      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
      spatial_cell->swap(vmesh,blockContainer);
   }
//   writeVelMesh(mpiGrid);

   phiprof::start("mapping");
   //if (P::xcells_ini > 1) {
      for (size_t c=0; c<local_cells.size(); ++c) {
         if (do_translate_cell(mpiGrid[local_cells[c]])) {
            trans_map_1d(mpiGrid,local_cells[c],dim,dt);
         }
      }
   //}
   phiprof::stop("mapping");

   for (size_t c=0; c<local_cells.size(); ++c) {
      SpatialCell* spatial_cell = mpiGrid[local_cells[c]];
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh_temporary();
      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
      //spatial_cell->swap(vmesh,blockContainer);
      vmesh.clear();
      blockContainer.clear();
   }
   writeVelMesh(mpiGrid);
   
   if (mpiGrid.get_rank() == 0) {
      cout << "Total mass (dim=" << dim << ") is " << calculateTotalMass(mpiGrid) << endl;
   }
   
   --dim;
   if (dim < 0) dim = 2;

   // Mapping complete, update moments //
   phiprof::start("compute-moments-n-maxdt");
   // Note: Parallelization over blocks is not thread-safe
   #pragma omp  parallel for
   for (size_t c=0; c<local_cells.size(); ++c) {
      SpatialCell* SC=mpiGrid[local_cells[c]];
      
      const Real dx=SC->parameters[CellParams::DX];
      const Real dy=SC->parameters[CellParams::DY];
      const Real dz=SC->parameters[CellParams::DZ];
      SC->parameters[CellParams::RHO_R  ] = 0.0;
      SC->parameters[CellParams::VX_R] = 0.0;
      SC->parameters[CellParams::VY_R] = 0.0;
      SC->parameters[CellParams::VZ_R] = 0.0;
      SC->parameters[CellParams::P_11_R ] = 0.0;
      SC->parameters[CellParams::P_22_R ] = 0.0;
      SC->parameters[CellParams::P_33_R ] = 0.0;
      
      //Reset spatial max DT
      SC->parameters[CellParams::MAXRDT]=numeric_limits<Real>::max();
      for (vmesh::LocalID block_i=0; block_i<SC->get_number_of_velocity_blocks(); ++block_i) {
         const Real* const blockParams = SC->get_block_parameters(block_i);

         //compute maximum dt. Algorithm has a CFL condition, since it
         //is written only for the case where we have a stencil
         //supporting max translation of one cell
         for (unsigned int i=0; i<WID;i+=WID-1) {
            const Real Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
            const Real Vy = blockParams[BlockParams::VYCRD] + (i+HALF)*blockParams[BlockParams::DVY];
            const Real Vz = blockParams[BlockParams::VZCRD] + (i+HALF)*blockParams[BlockParams::DVZ];
            
            if(fabs(Vx)!=ZERO) SC->parameters[CellParams::MAXRDT]=min(dx/fabs(Vx),SC->parameters[CellParams::MAXRDT]);
            if(fabs(Vy)!=ZERO) SC->parameters[CellParams::MAXRDT]=min(dy/fabs(Vy),SC->parameters[CellParams::MAXRDT]);
            if(fabs(Vz)!=ZERO) SC->parameters[CellParams::MAXRDT]=min(dz/fabs(Vz),SC->parameters[CellParams::MAXRDT]);
         }
         
         //compute first moments for this block
         if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
           cpu_calcVelocityFirstMoments(
                                        SC,
                                        block_i,			 
                                        CellParams::RHO_R,
                                        CellParams::VX_R,
                                        CellParams::VY_R,
                                        CellParams::VZ_R
                                       );   //set first moments after translation
      }
      // Second iteration needed as rho has to be already computed when computing pressure
      for (vmesh::LocalID block_i=0; block_i< SC->get_number_of_velocity_blocks(); ++block_i){
         //compute second moments for this block
         if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
           cpu_calcVelocitySecondMoments(
                                         SC,
                                         block_i,			  
                                         CellParams::RHO_R,
                                         CellParams::VX_R,
                                         CellParams::VY_R,
                                         CellParams::VZ_R,
                                         CellParams::P_11_R,
                                         CellParams::P_22_R,
                                         CellParams::P_33_R
                                        );   //set second moments after translation
      }
   }
   phiprof::stop("compute-moments-n-maxdt");
   //phiprof::stop("semilag-trans");
    */
}

/*
  --------------------------------------------------
  Acceleration (velocity space propagation)
  --------------------------------------------------
*/
void calculateAcceleration(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,Real dt) {
   return;
   /*
   typedef Parameters P;
   const vector<CellID>& cells = getLocalCells();
   vector<CellID> propagatedCells;
   // Iterate through all local cells and propagate distribution functions 
   // in velocity space. Ghost cells (spatial cells at the boundary of the simulation 
   // volume) do not need to be propagated:

   
   //    if(dt > 0) { // FIXME this has to be deactivated to support regular projects but it breaks test_trans support most likely, with this on dt stays 0
   //do not propagate for zero or negative dt. Typically dt==0 when
   //acceleration is turned off. 
   //Aet initial cells to propagate
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* SC = mpiGrid[cells[c]];
      //disregard boundary cells
      //do not integrate cells with no blocks  (well, do not computes in practice)
      if (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY &&
          SC->get_number_of_velocity_blocks() != 0) {
         propagatedCells.push_back(cells[c]);
      }
   }

   //Semilagrangian acceleration
   phiprof::start("semilag-acc");
   //#pragma omp parallel for schedule(dynamic,1)
   for (size_t c=0; c<propagatedCells.size(); ++c) {
      const CellID cellID = propagatedCells[c];
      //generate pseudo-random order which is always the same irrespectiive of parallelization, restarts, etc
      srand(P::tstep + cellID);
      uint map_order=rand()%3;
      phiprof::start("cell-semilag-acc");
      cpu_accelerate_cell(mpiGrid[cellID],map_order,dt);
      phiprof::stop("cell-semilag-acc");
   }
   phiprof::stop("semilag-acc");

   phiprof::start("Compute moments");
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      //compute moments after acceleration
      mpiGrid[cellID]->parameters[CellParams::RHO_V  ] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::VX_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::VY_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::VZ_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::P_11_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::P_22_V] = 0.0;
      mpiGrid[cellID]->parameters[CellParams::P_33_V] = 0.0;

      for (vmesh::LocalID block_i=0; block_i<mpiGrid[cellID]->get_number_of_velocity_blocks(); ++block_i) {
         cpu_calcVelocityFirstMoments(
                                      mpiGrid[cellID],
                                      block_i,
                                      CellParams::RHO_V,
                                      CellParams::VX_V,
                                      CellParams::VY_V,
                                      CellParams::VZ_V
                                     );   //set first moments after acceleration
      }

      // Second iteration needed as rho has to be already computed when computing pressure
      for (vmesh::LocalID block_i=0; block_i<mpiGrid[cellID]->get_number_of_velocity_blocks(); ++block_i) {
         cpu_calcVelocitySecondMoments(
                                       mpiGrid[cellID],
                                       block_i,
                                       CellParams::RHO_V,
                                       CellParams::VX_V,
                                       CellParams::VY_V,
                                       CellParams::VZ_V,
                                       CellParams::P_11_V,
                                       CellParams::P_22_V,
                                       CellParams::P_33_V
                                      );   //set second moments after acceleration
      }
   }
   phiprof::stop("Compute moments");
    */
}


/*--------------------------------------------------
  Functions for computing moments
  --------------------------------------------------*/
void calculateInterpolatedVelocityMoments(
                                          dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                          const int cp_rhom,
                                          const int cp_vx,
                                          const int cp_vy,
                                          const int cp_vz,
                                          const int cp_rhoq,
                                          const int cp_p11,
                                          const int cp_p22,
                                          const int cp_p33
                                         ) {
   const vector<CellID>& cells = getLocalCells();
   
   //Iterate through all local cells (excl. system boundary cells):
#pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      if(SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SC->parameters[cp_rhom  ] = 0.5* ( SC->parameters[CellParams::RHOM_R] + SC->parameters[CellParams::RHOM_V] );
         SC->parameters[cp_vx] = 0.5* ( SC->parameters[CellParams::VX_R] + SC->parameters[CellParams::VX_V] );
         SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R] + SC->parameters[CellParams::VY_V] );
         SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R] + SC->parameters[CellParams::VZ_V] );
         SC->parameters[cp_rhoq  ] = 0.5* ( SC->parameters[CellParams::RHOQ_R] + SC->parameters[CellParams::RHOQ_V] );
         SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R] + SC->parameters[CellParams::P_11_V] );
         SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R] + SC->parameters[CellParams::P_22_V] );
         SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R] + SC->parameters[CellParams::P_33_V] );
      }
   }
}

void calculateCellVelocityMoments(SpatialCell* SC,
                                  bool doNotSkip // default: false
                                 ) {
   /*
   // if doNotSkip == true then the first clause is false and we will never return, i.e. always compute
   // otherwise we skip DO_NOT_COMPUTE cells
   // or boundary cells of layer larger than 1
   if (!doNotSkip &&
       (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
	(SC->sysBoundaryLayer != 1  &&
	 SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
       ) return;

   // Clear old moments
   Real* cellParams = SC->get_cell_parameters();
   cellParams[CellParams::RHO  ] = 0.0;
   cellParams[CellParams::VX] = 0.0;
   cellParams[CellParams::VY] = 0.0;
   cellParams[CellParams::VZ] = 0.0;
   cellParams[CellParams::P_11 ] = 0.0;
   cellParams[CellParams::P_22 ] = 0.0;
   cellParams[CellParams::P_33 ] = 0.0;

   // Calculate first moments
   for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Temporary array for storing this species' contribution
      Real array[4];
      for (int i=0; i<4; ++i) array[i] = 0;
      
      // Pointers to this species' data
      const Realf* data = SC->get_data(popID);
      const Real* blockParams = SC->get_block_parameters(popID);
      
      for (vmesh::LocalID blockLID=0; blockLID<SC->get_number_of_velocity_blocks(popID); ++blockLID) {
         blockVelocityFirstMoments(
                  data,
                  blockParams,
                  array
         );
         data += SIZE_VELBLOCK;
         blockParams += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }
      
      const Real massRatio = getObjectWrapper().particleSpecies[popID].mass / physicalconstants::MASS_PROTON;
      cellParams[CellParams::RHO  ] += array[0]*massRatio;
      cellParams[CellParams::VX] += array[1]*massRatio;
      cellParams[CellParams::VY] += array[2]*massRatio;
      cellParams[CellParams::VZ] += array[3]*massRatio;
   } // for-loop over particle species

   // Second iteration needed as rho has to be already computed when computing pressure
   for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Temporary array for storing this species' contribution
      Real array[3];
      for (int i=0; i<3; ++i) array[i] = 0;
      
      // Pointers to this species' data
      const Realf* data = SC->get_data(popID);
      const Real* blockParams = SC->get_block_parameters(popID);
      
      for (vmesh::LocalID blockLID=0; blockLID<SC->get_number_of_velocity_blocks(popID); ++blockLID) {
         blockVelocitySecondMoments(
                  data,
                  blockParams,
                  cellParams,
                  CellParams::VX,
                  CellParams::VY,
                  CellParams::VZ,
                  array
         );
         data += SIZE_VELBLOCK;
         blockParams += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }
      
      cellParams[CellParams::P_11] += array[0]*getObjectWrapper().particleSpecies[popID].mass;
      cellParams[CellParams::P_22] += array[1]*getObjectWrapper().particleSpecies[popID].mass;
      cellParams[CellParams::P_33] += array[2]*getObjectWrapper().particleSpecies[popID].mass;
   } // for-loop over particle species
    */
}

void calculateInitialVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   /*const vector<CellID>& cells = getLocalCells();
   phiprof::start("Calculate moments"); 
   // Iterate through all local cells (incl. system boundary cells):
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      calculateCellVelocityMoments(SC);
      // WARNING the following is sane as this function is only called by initializeGrid.
      // We need initialized _DT2 values for the dt=0 field propagation done in the beginning.
      // Later these will be set properly.
      SC->parameters[CellParams::RHO_DT2] = SC->parameters[CellParams::RHO];
      SC->parameters[CellParams::VX_DT2] = SC->parameters[CellParams::VX];
      SC->parameters[CellParams::VY_DT2] = SC->parameters[CellParams::VY];
      SC->parameters[CellParams::VZ_DT2] = SC->parameters[CellParams::VZ];
      SC->parameters[CellParams::P_11_DT2] = SC->parameters[CellParams::P_11];
      SC->parameters[CellParams::P_22_DT2] = SC->parameters[CellParams::P_22];
      SC->parameters[CellParams::P_33_DT2] = SC->parameters[CellParams::P_33];

   }
   phiprof::stop("Calculate moments"); 
   */
}
