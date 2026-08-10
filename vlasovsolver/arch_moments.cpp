/*
 * This file is part of Vlasiator.
 * Copyright 2024-2025 University of Helsinki, CSC
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

#include <phiprof.hpp>
#include "arch_moments.h"
#include "vlasovmover.h"
#include "../object_wrapper.h"
#include "../fieldsolver/fs_common.h" // divideIfNonZero()

#ifdef USE_GPU
#include "gpu_moments.h"
#endif

using namespace std;

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include contributions from
 * all existing particle populations.
 * @param cell Spatial cell.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param doNotSkip If false, DO_NOT_COMPUTE cells are skipped.*/
void calculateCellMoments(spatial_cell::SpatialCell* cell,
                          const bool& computeSecond,
                          const bool& computePopulationMomentsOnly,
                          const bool& doNotSkip) {

   // Called once per cell. If doNotSkip == true, then DO_NOT_COMPUTE cells aren't skipped.
   if (!doNotSkip && cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
       return;
   }

   // Clear old moments to zero value
   if (computePopulationMomentsOnly == false) {
      cell->parameters[CellParams::RHOM  ] = 0.0;
      cell->parameters[CellParams::VX] = 0.0;
      cell->parameters[CellParams::VY] = 0.0;
      cell->parameters[CellParams::VZ] = 0.0;
      cell->parameters[CellParams::RHOQ  ] = 0.0;
      cell->parameters[CellParams::P_11] = 0.0;
      cell->parameters[CellParams::P_22] = 0.0;
      cell->parameters[CellParams::P_33] = 0.0;
      cell->parameters[CellParams::P_23] = 0.0;
      cell->parameters[CellParams::P_13] = 0.0;
      cell->parameters[CellParams::P_12] = 0.0;
   }

   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
      #endif
      const uint nBlocks = cell->get_velocity_mesh(popID)->size();
      Population &pop = cell->get_population(popID);
      if (nBlocks == 0) {
         pop.RHO = 0;
         for (int i=0; i<3; ++i) {
            pop.V[i]=0;
         }
         for (int i=0; i<nMom2; ++i) {
            pop.P[i]=0;
         }
         continue;
      }

      vmesh::MeshParameters& vMeshprint=vmesh::getMeshWrapper()->velocityMeshesCreation->at(popID);
      species::Species& species=getObjectWrapper().particleSpecies[popID];  

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real charge = getObjectWrapper().particleSpecies[popID].charge;

      // Temporary array where the moments for this species are accumulated
      Real array[nMom1] = {0};
	 
      int Ref=0;
      int MaxRef=0;
      if(P::activateVamr) {
	    Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	    MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	    if(Ref < MaxRef){
		  // Update the Refined parameter to know the cells that need to be integrated
	      changeRefined(cell,popID);
		}
	    // Calculate species' contribution to first velocity moments with vAMR
	    blockVelocityFirstMomentsVamr(blockContainer,
				      array,
				      nBlocks);
      }else {
	    // Calculate species' contribution to first velocity moments
	    blockVelocityFirstMoments(blockContainer,
				  array,
				  nBlocks);
      }
  

      pop.RHO = array[0];
      pop.V[0] = divideIfNonZero(array[1], array[0]);
      pop.V[1] = divideIfNonZero(array[2], array[0]);
      pop.V[2] = divideIfNonZero(array[3], array[0]);

      if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
	  //Sum of all the partial integrals and sharing of the final result between the different vAMR grids of the same species
	    for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	      Population &pop2 = cell->get_population(popID2);
	      pop.RHO += pop2.RHO;
	      pop.V[0] += pop2.V[0];
	      pop.V[1] += pop2.V[1];
	      pop.V[2] += pop2.V[2];
	    };
	    for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	      Population &pop2 = cell->get_population(popID2);
	      pop2.RHO =  pop.RHO;
	      pop2.V[0] = pop.V[0];
	      pop2.V[1] = pop.V[1];
	      pop2.V[2] = pop.V[2];
	    };
      }; 

      if (!computePopulationMomentsOnly) {
         // Store species' contribution to bulk velocity moments
         cell->parameters[CellParams::RHOM  ] += array[0]*mass;
         cell->parameters[CellParams::VX] += array[1]*mass;
         cell->parameters[CellParams::VY] += array[2]*mass;
         cell->parameters[CellParams::VZ] += array[3]*mass;
         cell->parameters[CellParams::RHOQ  ] += array[0]*charge;
      }
   } // for-loop over particle species

   if(!computePopulationMomentsOnly) {
      cell->parameters[CellParams::VX] = divideIfNonZero(cell->parameters[CellParams::VX], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VY] = divideIfNonZero(cell->parameters[CellParams::VY], cell->parameters[CellParams::RHOM]);
      cell->parameters[CellParams::VZ] = divideIfNonZero(cell->parameters[CellParams::VZ], cell->parameters[CellParams::RHOM]);
   }

   // Compute second moments only if requested
   if (computeSecond == false) {
      return;
   }

   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
      #endif
      const uint nBlocks = cell->get_velocity_mesh(popID)->size();
      if (nBlocks == 0) {
         continue;
      }

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;

      // Temporary array for storing moments
      Real array[nMom2] = {0};

      int Ref=0;
      int MaxRef=0;
      if(P::activateVamr) {
	    int Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	    int MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	    // Calculate species' contribution to second velocity moments with vAMR
	    blockVelocitySecondMomentsVamr(blockContainer,
				       cell->parameters[CellParams::VX],
				       cell->parameters[CellParams::VY],
				       cell->parameters[CellParams::VZ],
				       array,
				       nBlocks);
      }else{
	    // Calculate species' contribution to second velocity moments
	    blockVelocitySecondMoments(blockContainer,
				   cell->parameters[CellParams::VX],
				   cell->parameters[CellParams::VY],
				   cell->parameters[CellParams::VZ],
				   array,
				   nBlocks);
      }
      // Store species' contribution to bulk velocity moments
      Population &pop = cell->get_population(popID);
      for (size_t i=0; i<nMom2; ++i) {
         pop.P[i] = mass * array[i];
      }

      if (!computePopulationMomentsOnly) {
         cell->parameters[CellParams::P_11] += pop.P[0];
         cell->parameters[CellParams::P_22] += pop.P[1];
         cell->parameters[CellParams::P_33] += pop.P[2];
         cell->parameters[CellParams::P_23] += pop.P[3];
         cell->parameters[CellParams::P_13] += pop.P[4];
         cell->parameters[CellParams::P_12] += pop.P[5];
      }

      if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
	  //Sum of all the partial integrals and sharing of the final result between the different vAMR grids of the same species
        for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	      Population &pop2 = cell->get_population(popID2);
	      for (size_t i=0; i<nMom2; ++i) {
	        pop.P[i] +=   pop2.P[i];
	      };
	    };
	    for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	      Population &pop2 = cell->get_population(popID2);
	      for (size_t i=0; i<nMom2; ++i) {
	        pop2.P[i] = pop.P[i];
	      };
	    };
      };
      
   } // for-loop over particle species

}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _R variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param initialCompute If true, force re-calculation of outflow L1 sysboundary cell moments.
  (otherwise skipped as their VDF contents are not kept up to date)
*/
void calculateMoments_R(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond,
   const bool initialCompute) {

   // override with optimized GPU version to launch
   // single kernel accessing all cells at once (10x faster)
   #ifdef USE_GPU
   gpu_calculateMoments_R(mpiGrid,cells,computeSecond,initialCompute);
   return;
   #endif

   phiprof::Timer computeMomentsTimer {"Compute _R moments"};
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }
         // Clear old moments to zero value
         if (popID == 0) {
            cell->parameters[CellParams::RHOM_R  ] = 0.0;
            cell->parameters[CellParams::VX_R] = 0.0;
            cell->parameters[CellParams::VY_R] = 0.0;
            cell->parameters[CellParams::VZ_R] = 0.0;
            cell->parameters[CellParams::RHOQ_R  ] = 0.0;
            cell->parameters[CellParams::P_11_R] = 0.0;
            cell->parameters[CellParams::P_22_R] = 0.0;
            cell->parameters[CellParams::P_33_R] = 0.0;
            cell->parameters[CellParams::P_23_R] = 0.0;
            cell->parameters[CellParams::P_13_R] = 0.0;
            cell->parameters[CellParams::P_12_R] = 0.0;
         }

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         Population &pop = cell->get_population(popID);
         if (nBlocks == 0) {
            pop.RHO_R = 0;
            for (int i=0; i<3; ++i) {
               pop.V_R[i]=0;
            }
            for (int i=0; i<nMom2; ++i) {
               pop.P_R[i]=0;
            }
            continue;
         }
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array where the moments for this species are accumulated
         Real array[nMom1] = {0};

         int Ref=0;
	     int MaxRef=0;
	     if(P::activateVamr) {
	       Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	       MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	       if(Ref < MaxRef){
   			 // Update the Refined parameter to know the cells that need to be integrated
	         changeRefined(cell,popID);
		   }
	       // Calculate species' contribution to first velocity moments with vAMR
	       blockVelocityFirstMomentsVamr(blockContainer,
					 array,
					 nBlocks);
	     }else {
	       // Calculate species' contribution to first velocity moments
	       blockVelocityFirstMoments(blockContainer,
				     array,
				     nBlocks);
	     }
	
         // Store species' contribution to bulk velocity moments
         pop.RHO_R = array[0];
         pop.V_R[0] = divideIfNonZero(array[1], array[0]);
         pop.V_R[1] = divideIfNonZero(array[2], array[0]);
         pop.V_R[2] = divideIfNonZero(array[3], array[0]);

	     if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
		 //Sum of all the partial integrals and sharing of the final result between the different vAMR grids of the same species
	       for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	         Population &pop2 = cell->get_population(popID2);
	         pop.RHO_R += pop2.RHO_R;
	         pop.V_R[0] += pop2.V_R[0];
	         pop.V_R[1] += pop2.V_R[1];
	         pop.V_R[2] += pop2.V_R[2];
	       };
	       for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	         Population &pop2 = cell->get_population(popID2);
	         pop2.RHO_R =  pop.RHO_R;
	         pop2.V_R[0] = pop.V_R[0];
	         pop2.V_R[1] = pop.V_R[1];
	         pop2.V_R[2] = pop.V_R[2];
	       };
	     };

         cell->parameters[CellParams::RHOM_R  ] += array[0]*mass;
         cell->parameters[CellParams::VX_R] += array[1]*mass;
         cell->parameters[CellParams::VY_R] += array[2]*mass;
         cell->parameters[CellParams::VZ_R] += array[3]*mass;
         cell->parameters[CellParams::RHOQ_R  ] += array[0]*charge;
      } // for-loop over spatial cells
   } // for-loop over particle species

   #pragma omp parallel for schedule(static)
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* cell = mpiGrid[cells[c]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
         continue;
      }
      cell->parameters[CellParams::VX_R] = divideIfNonZero(cell->parameters[CellParams::VX_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VY_R] = divideIfNonZero(cell->parameters[CellParams::VY_R], cell->parameters[CellParams::RHOM_R]);
      cell->parameters[CellParams::VZ_R] = divideIfNonZero(cell->parameters[CellParams::VZ_R], cell->parameters[CellParams::RHOM_R]);
   }

   // Compute second moments only if requested.
   if (computeSecond == false) {
      return;
   }

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         if (nBlocks == 0) {
            continue;
         }
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;

         // Temporary array where species' contribution to 2nd moments is accumulated
         Real array[nMom2] = {0};

         int Ref=0;
	 	 int MaxRef=0;
	 	 if(P::activateVamr) {
	   	    int Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	   	    int MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	   	    // Calculate species' contribution to second velocity moments with Vam
	   	    blockVelocitySecondMomentsVamr(blockContainer,
			 		     cell->parameters[CellParams::VX],
					     cell->parameters[CellParams::VY],
					     cell->parameters[CellParams::VZ],
					     array,
					     nBlocks);
	 	 }else{
	   	    // Calculate species' contribution to second velocity moments
	   	    blockVelocitySecondMoments(blockContainer,
				     cell->parameters[CellParams::VX],
				     cell->parameters[CellParams::VY],
				     cell->parameters[CellParams::VZ],
				     array,
				     nBlocks);
	 	 }
	 
         // Store species' contribution to 2nd bulk velocity moments
         Population &pop = cell->get_population(popID);
         for (size_t i = 0; i < nMom2; ++i) {
            pop.P_R[i] = mass * array[i];
         }

         cell->parameters[CellParams::P_11_R] += pop.P_R[0];
         cell->parameters[CellParams::P_22_R] += pop.P_R[1];
         cell->parameters[CellParams::P_33_R] += pop.P_R[2];
         cell->parameters[CellParams::P_23_R] += pop.P_R[3];
         cell->parameters[CellParams::P_13_R] += pop.P_R[4];
         cell->parameters[CellParams::P_12_R] += pop.P_R[5];

	 	 if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
	 	 //Sum of all the partial integrals and sharing of the final result between the different vAMR grids of the same species
	   	    for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	     	   Population &pop2 = cell->get_population(popID2);
	     	   for (size_t i=0; i<nMom2; ++i) {
	       	      pop.P_R[i] +=   pop2.P_R[i];
	     	   };
	   		};
	   		for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	     	   Population &pop2 = cell->get_population(popID2);
	     	   for (size_t i=0; i<nMom2; ++i) {
	       	      pop2.P_R[i] = pop.P_R[i];
	     	   };
	   		};
	 	 };
	 
      } // for-loop over spatial cells
   } // for-loop over particle species
}

/** Calculate zeroth, first, and (possibly) second bulk velocity moments for the
 * given spatial cell. Additionally, for each species, calculate the maximum
 * spatial time step so that CFL(spatial)=1. The calculated moments include
 * contributions from all existing particle populations. The calculated moments
 * are stored to SpatialCell::parameters in _V variables.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing the spatial cells to be calculated.
 * @param computeSecond If true, second velocity moments are calculated.
 * @param initialCompute If true, force re-calculation of outflow L1 sysboundary cell moments.
  (otherwise skipped as their VDF contents are not kept up to date)
*/
void calculateMoments_V(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond,
   const bool initialCompute) {

   // override with optimized GPU version to launch
   // single kernel accessing all cells at once (10x faster)
   #ifdef USE_GPU
   gpu_calculateMoments_V(mpiGrid,cells,computeSecond,initialCompute);
   return;
   #endif

   phiprof::Timer computeMomentsTimer {"Compute _V moments"};
   // Loop over all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }

         // Clear old moments to zero value
         if (popID == 0) {
            cell->parameters[CellParams::RHOM_V  ] = 0.0;
            cell->parameters[CellParams::VX_V] = 0.0;
            cell->parameters[CellParams::VY_V] = 0.0;
            cell->parameters[CellParams::VZ_V] = 0.0;
            cell->parameters[CellParams::RHOQ_V  ] = 0.0;
            cell->parameters[CellParams::P_11_V] = 0.0;
            cell->parameters[CellParams::P_22_V] = 0.0;
            cell->parameters[CellParams::P_33_V] = 0.0;
            cell->parameters[CellParams::P_23_V] = 0.0;
            cell->parameters[CellParams::P_13_V] = 0.0;
            cell->parameters[CellParams::P_12_V] = 0.0;
         }

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         Population &pop = cell->get_population(popID);
         if (nBlocks == 0) {
            pop.RHO_V = 0;
            for (int i=0; i<3; ++i) {
               pop.V_V[i]=0;
            }
            for (int i=0; i<nMom2; ++i) {
               pop.P_V[i]=0;
            }
            continue;
         }

         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array where the moments for this species are accumulated
         Real array[nMom1] = {0};

         int Ref=0;
	     int MaxRef=0;
	     if(P::activateVamr) {
	       Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	       MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	       if(Ref < MaxRef){
			 // Update the Refined parameter to know the cells that need to be integrated
	         changeRefined(cell,popID);
		   }
	       // Calculate species' contribution to first velocity moments with vAMR
	       blockVelocityFirstMomentsVamr(blockContainer,
					 array,
					 nBlocks);
	     }else {
	       // Calculate species' contribution to first velocity moments
	       blockVelocityFirstMoments(blockContainer,
				     array,
				     nBlocks);
	     }
	 
         // Store species' contribution to bulk velocity moments
         pop.RHO_V = array[0];
         pop.V_V[0] = divideIfNonZero(array[1], array[0]);
         pop.V_V[1] = divideIfNonZero(array[2], array[0]);
         pop.V_V[2] = divideIfNonZero(array[3], array[0]);
      
	     if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
		 //Sum of all the partial integrals and sharing of the final result between the different vAMR grids of the same species
	       for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	         Population &pop2 = cell->get_population(popID2);
	         pop.RHO_V += pop2.RHO_V;
	         pop.V_V[0] += pop2.V_V[0];
	         pop.V_V[1] += pop2.V_V[1];
	         pop.V_V[2] += pop2.V_V[2];
	       };
	       for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	         Population &pop2 = cell->get_population(popID2);
	         pop2.RHO_V =  pop.RHO_V;
	         pop2.V_V[0] = pop.V_V[0];
	         pop2.V_V[1] = pop.V_V[1];
	         pop2.V_V[2] = pop.V_V[2];
	       };
	     };

         cell->parameters[CellParams::RHOM_V  ] += array[0]*mass;
         cell->parameters[CellParams::VX_V] += array[1]*mass;
         cell->parameters[CellParams::VY_V] += array[2]*mass;
         cell->parameters[CellParams::VZ_V] += array[3]*mass;
         cell->parameters[CellParams::RHOQ_V  ] += array[0]*charge;

      } // for-loop over spatial cells
   } // for-loop over particle species

   #pragma omp parallel for schedule(static)
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* cell = mpiGrid[cells[c]];
      if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
         continue;
      }
      cell->parameters[CellParams::VX_V] = divideIfNonZero(cell->parameters[CellParams::VX_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VY_V] = divideIfNonZero(cell->parameters[CellParams::VY_V], cell->parameters[CellParams::RHOM_V]);
      cell->parameters[CellParams::VZ_V] = divideIfNonZero(cell->parameters[CellParams::VZ_V], cell->parameters[CellParams::RHOM_V]);
   }

   // Compute second moments only if requested
   if (computeSecond == false) {
      return;
   }

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      #pragma omp parallel for schedule(dynamic,1)
      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];

         if (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         if (cell->sysBoundaryFlag == sysboundarytype::OUTFLOW && cell->sysBoundaryLayer != 1 && !initialCompute) { // these should have been handled by the boundary code
            continue;
         }

         #ifdef USE_GPU
         vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
         #endif
         const uint nBlocks = cell->get_velocity_mesh(popID)->size();
         if (nBlocks == 0) {
            continue;
         }

         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         const Real charge = getObjectWrapper().particleSpecies[popID].charge;

         // Temporary array where moments are stored
         Real array[nMom2] = {0};

         int Ref=0;
	     int MaxRef=0;
	     if(P::activateVamr) {
	       int Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	       int MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	       // Calculate species' contribution to second velocity moments with vAMR
	       blockVelocitySecondMomentsVamr(blockContainer,
					  cell->parameters[CellParams::VX],
					  cell->parameters[CellParams::VY],
					  cell->parameters[CellParams::VZ],
					  array,
					  nBlocks);
	    }else{
	       // Calculate species' contribution to second velocity moments
	       blockVelocitySecondMoments(blockContainer,
				      cell->parameters[CellParams::VX],
				      cell->parameters[CellParams::VY],
				      cell->parameters[CellParams::VZ],
				      array,
				      nBlocks);
	    }
      
         // Store species' contribution to 2nd bulk velocity moments
         Population &pop = cell->get_population(popID);
         for (size_t i = 0; i < nMom2; ++i) {
            pop.P_V[i] = mass * array[i];
         }

         cell->parameters[CellParams::P_11_V] += pop.P_V[0];
         cell->parameters[CellParams::P_22_V] += pop.P_V[1];
         cell->parameters[CellParams::P_33_V] += pop.P_V[2];
         cell->parameters[CellParams::P_23_V] += pop.P_V[3];
         cell->parameters[CellParams::P_13_V] += pop.P_V[4];
         cell->parameters[CellParams::P_12_V] += pop.P_V[5];
	 
	    if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
		//Sum of all the partial integrals and sharing of the final result between the different vAMR grids of the same species
	      for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	        Population &pop2 = cell->get_population(popID2);
	        for (size_t i=0; i<nMom2; ++i) {
	          pop.P_V[i] +=   pop2.P_V[i];
	        };
	      };
	      for (uint popID2=(popID-MaxRef); popID2<popID; ++popID2) {
	        Population &pop2 = cell->get_population(popID2);
	        for (size_t i=0; i<nMom2; ++i) {
	          pop2.P_V[i] = pop.P_V[i];
	        };
	      };
	    };
	 
      } // for-loop over spatial cells
   } // for-loop over particle species
}

// Update communication from level R+1 to level R
void vamr_transfer_values(
			  dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			  const std::vector<CellID>& cells,
			  const uint popID) {

  // Loop over all particle species
#pragma omp parallel for schedule(dynamic,1)
  for (size_t c=0; c<cells.size(); ++c) {
    SpatialCell* cell = mpiGrid[cells[c]];
	Realf minValue = 0; //cell->getVelocityBlockMinValue(popID);
	  
#ifdef USE_GPU
      vmesh::VelocityMesh* vmesh    = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf    = cell->get_velocity_mesh(popID+1);
      vmesh::VelocityBlockContainer* blockContainer = cell->dev_get_velocity_blocks(popID);
      vmesh::VelocityBlockContainer* blockContainerraf = cell->dev_get_velocity_blocks(popID+1);
#else
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf    = cell->get_velocity_mesh(popID+1);
      vmesh::VelocityBlockContainer* blockContainer = cell->get_velocity_blocks(popID);
      vmesh::VelocityBlockContainer* blockContainerraf = cell->get_velocity_blocks(popID+1);
#endif
    const uint nBlocks = cell->get_velocity_mesh(popID)->size();
    Population &pop = cell->get_population(popID);
    if (nBlocks == 0) {
	  continue;
    }

    Realf *data = blockContainer->getData();
    Realf *dataraf = blockContainerraf->getData();
	uint8_t *ghost = blockContainerraf->getGhost();
  
    for (vmesh::LocalID localID=0; localID<vmesh->size(); ++localID) {
	  const vmesh::GlobalID globalID = vmesh->getGlobalID(localID);	   
	   
	  vmesh::LocalID Indices[3];
	  vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);

	  for (int i=0; i<2; ++i) {
	    for (int j=0; j<2; ++j) {
	      for (int k=0; k<2; ++k) {
	        //Indices of the refined blocks
	        //Each coarse block contains 8 refined blocks
	        vmesh::LocalID Indicesraf[3];
	        Indicesraf[0] = 2*Indices[0]+i ;
	        Indicesraf[1] = 2*Indices[1]+j ;
	        Indicesraf[2] = 2*Indices[2]+k ;

	        const vmesh::GlobalID globalIDraf=vmeshraf->getGlobalID(Indicesraf);
			//Check if the refined cell exists
	        if (globalIDraf ==  vmeshraf->invalidGlobalID()) {
		      continue;
	        }else{  
		      const vmesh::LocalID  localIDraf=vmeshraf->getLocalID(globalIDraf);
		      if (localIDraf == vmeshraf->invalidLocalID()) {
		        continue;
		      }else{
		        for (int i2=0; i2<2; ++i2) {
		          for (int j2=0; j2<2; ++j2) {
		            for (int k2=0; k2<2; ++k2) {
					//Loop over the 8 coarsed cells located on the WID3 refined cells (the refined block)
	        	      Realf summ=0;
					  //  if( data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]>cell->getVelocityBlockMinValue(popID+1)){ //old criteria used to improve communication
			          Realf datasave= data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
			          data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]=0;
			          for (int i3=0; i3<2; ++i3) {
			            for (int j3=0; j3<2; ++j3) {
			              for (int k3=0; k3<2; ++k3) {
						  //Loop over the 8 refined cells contained in the 1 coarsed cell
							 if(dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]>minValue){
							 // if(ghost[localIDraf]==1 && dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]>minValue){ //Another criteria
				              data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]+= dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]/8.0;
				              summ+=1.0;
			                }else{
							  dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]=datasave;
							}
			              }
			            }
			          }
			         if (summ!=8.0){
			            data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]=datasave;
			          }
					/*}else{
				      for (int i3=0; i3<2; ++i3) {
				        for (int j3=0; j3<2; ++j3) {
				      	  for (int k3=0; k3<2; ++k3) {
						    dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]=data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
				      	  }
				    	}
				  	  }
					  }*/ //Removed because the boundaries were looking too coarse 
						
		            };
		          };
		        };//loop over refined cells
		      };
	        };
	      };//loop over refined blocks
	    };
	  };
    };
  }
}

// Update the Refined parameter to know the cells that need to be integrated
 void changeRefined(spatial_cell::SpatialCell* cell,
 const uint popID){

     vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
     vmesh::VelocityMesh* vmeshraf = cell->get_velocity_mesh(popID+1);
     
     uint8_t *refined =cell->get_velocity_blocks(popID)->getRefined();
//Parallelisation is just in case; this function may already be part of a parallelised region     
#pragma omp parallel for schedule(dynamic,1)       
     for (vmesh::LocalID localID=0; localID<vmesh->size(); ++localID) {
       const vmesh::GlobalID globalID = vmesh->getGlobalID(localID);	   
	   
       vmesh::LocalID Indices[3];
       vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);

       for (int i=0; i<2; ++i) {
	     for (int j=0; j<2; ++j) {
	       for (int k=0; k<2; ++k) {
	       //Indices of the refined blocks
	       //Each coarse block contains 8 refined blocks
	       vmesh::LocalID Indicesraf[3];
	       Indicesraf[0] = 2*Indices[0]+i ;
	       Indicesraf[1] = 2*Indices[1]+j ;
	       Indicesraf[2] = 2*Indices[2]+k ;

	       const vmesh::GlobalID globalIDraf=vmeshraf->getGlobalID(Indicesraf);
		   //We will check if the refined block exists
	       if (globalIDraf ==  vmeshraf->invalidGlobalID()) {
	         //each refined block represent (WID/2)^3 coarse cells (8 in our case) 
	         refined[localID*WID3+cellIndex(2*i,2*j,2*k)]=false;
	         refined[localID*WID3+cellIndex(2*i+1,2*j,2*k)]=false;
	         refined[localID*WID3+cellIndex(2*i,2*j+1,2*k)]=false;
	         refined[localID*WID3+cellIndex(2*i,2*j,2*k+1)]=false;
	         refined[localID*WID3+cellIndex(2*i+1,2*j+1,2*k)]=false;
	         refined[localID*WID3+cellIndex(2*i,2*j+1,2*k+1)]=false;
	         refined[localID*WID3+cellIndex(2*i+1,2*j,2*k+1)]=false;
	         refined[localID*WID3+cellIndex(2*i+1,2*j+1,2*k+1)]=false; 
	       }else{  
	         const vmesh::LocalID  localIDraf=vmeshraf->getLocalID(globalIDraf);
	         if (localIDraf == vmeshraf->invalidLocalID()) {
		       refined[localID*WID3+cellIndex(2*i,2*j,2*k)]=false;
		       refined[localID*WID3+cellIndex(2*i+1,2*j,2*k)]=false;
		       refined[localID*WID3+cellIndex(2*i,2*j+1,2*k)]=false;
		       refined[localID*WID3+cellIndex(2*i,2*j,2*k+1)]=false;
		       refined[localID*WID3+cellIndex(2*i+1,2*j+1,2*k)]=false;
		       refined[localID*WID3+cellIndex(2*i,2*j+1,2*k+1)]=false;
		       refined[localID*WID3+cellIndex(2*i+1,2*j,2*k+1)]=false;
		       refined[localID*WID3+cellIndex(2*i+1,2*j+1,2*k+1)]=false; 
	         }else{
			   //The refine block cell
		       refined[localID*WID3+cellIndex(2*i,2*j,2*k)]=true;
		       refined[localID*WID3+cellIndex(2*i+1,2*j,2*k)]=true;
		       refined[localID*WID3+cellIndex(2*i,2*j+1,2*k)]=true;
		       refined[localID*WID3+cellIndex(2*i,2*j,2*k+1)]=true;
		       refined[localID*WID3+cellIndex(2*i+1,2*j+1,2*k)]=true;
		       refined[localID*WID3+cellIndex(2*i,2*j+1,2*k+1)]=true;
		       refined[localID*WID3+cellIndex(2*i+1,2*j,2*k+1)]=true;
		       refined[localID*WID3+cellIndex(2*i+1,2*j+1,2*k+1)]=true; 
	         };
	       };
	     };
	   };
     };
   };
  }

//Every vamr_refinedStep we check all the velocity cells with the vAMR criterion with the 1st order
void RefinedOrder1(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells){

#pragma omp parallel for schedule(dynamic,1)
for (size_t c=0; c<cells.size(); ++c) {
SpatialCell* cell = mpiGrid[cells[c]];
std::unordered_set<vmesh::GlobalID> ListBlockExist[getObjectWrapper().particleSpecies.size()];
  
  for (int popID=(getObjectWrapper().particleSpecies.size()-2); popID>-1; --popID) {

    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 && getObjectWrapper().particleSpecies[popID].RefinementLevel<getObjectWrapper().particleSpecies[popID].MaxRefinementLevel){
       
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf    = cell->get_velocity_mesh(popID+1);
      Realf *data = cell->get_velocity_blocks(popID)->getData();
      
      for (vmesh::LocalID localID=0; localID<vmesh->size(); ++localID) { 
		vmesh::GlobalID globalID = vmesh->getGlobalID(localID);
		vmesh::LocalID Indices[3];
		vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);

		if (getObjectWrapper().particleSpecies[popID].RefinementLevel==0) {
	  		ListBlockExist[popID].insert(globalID);
		}
      
		for (int i=0; i<2; ++i) {
	  	  for (int j=0; j<2; ++j) {
	    	for (int k=0; k<2; ++k) {
			//Loop over the future refined blocks R+1
	      	  Realf Datagros = 0; 
	   
	          for (int i2=0; i2<2; ++i2) {
				for (int j2=0; j2<2; ++j2) {
		       	  for (int k2=0; k2<2; ++k2) {
				  //Sum over the 8 cells of level R in order to reproduce the coarse cell R-1
		    	    Datagros += data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  }
				}
	      	  }
	      	  Datagros/=8.0;
	    	      	     		    
	      	  Realf D = abs( data[localID*WID3+cellIndex(1+i,1+j,1+k)] - Datagros ); 
			  // The idea is to always compare the 8 central cells of level R with the reproduce R-1 cells  

	      	  vmesh::LocalID Indicesraf[3];
	      	  Indicesraf[0] = 2*Indices[0]+i;
	      	  Indicesraf[1] = 2*Indices[1]+j;
	      	  Indicesraf[2] = 2*Indices[2]+k;

              Realf Dcomp = cell->getVelocityBlockMinValue(popID);

              if(getObjectWrapper().particleSpecies[popID].CriteriaMethod==1){ //epsilon
                Dcomp= getObjectWrapper().particleSpecies[popID].CriteriaValue;
              }else if(getObjectWrapper().particleSpecies[popID].CriteriaMethod==2){ //epsilon*2^-R
 	            Dcomp= getObjectWrapper().particleSpecies[popID].CriteriaValue/(1u << getObjectWrapper().particleSpecies[popID].RefinementLevel);
              }
              if (D > Dcomp){
				// We should create a new block for R+1
				int addWidthV = 1; //getObjectWrapper().particleSpecies[popID+1].sparseBlockAddWidthV; 
				for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
		  		  for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
		    		for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
		      		  const vmesh::GlobalID globalIDraf = vmeshraf->getGlobalID(Indicesraf[0]+offset_vx,Indicesraf[1]+offset_vy,Indicesraf[2]+offset_vz);
		      		  if (globalIDraf==  vmeshraf->invalidGlobalID()) {
						// std::cout<< " GlobalID bug not normal" << "Indices[0]+offset_vx" << Indices[0]+offset_vx << "Indices[1]+offset_vy " << Indices[1]+offset_vy << "Indices[1]+offset_vy" << Indices[2]+offset_vz <<std::endl;
		      		  }else{
		      			ListBlockExist[popID+1].insert(globalIDraf);
		    		  }
						
		  			}
				  } 
	     		}
				  
			  }
	      	
	    	}
	  	 }
	  } 

    }
    
   }
    
  }

  //We will now ensure that each level R+1 exist on the level R with the proper ghost cells between both
  
  for (int popID=(getObjectWrapper().particleSpecies.size()-2); popID>-1; --popID) {

    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 && getObjectWrapper().particleSpecies[popID].RefinementLevel<getObjectWrapper().particleSpecies[popID].MaxRefinementLevel){

      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf = cell->get_velocity_mesh(popID+1);
      
      for (vmesh::GlobalID globalIDraf : ListBlockExist[popID+1]) {
	
	    vmesh::LocalID Indicesraf[3];
		vmeshraf->getIndices(globalIDraf, Indicesraf[0], Indicesraf[1], Indicesraf[2]);

		vmesh::LocalID Indices[3];
		Indices[0] = Indicesraf[0]/2;
		Indices[1] = Indicesraf[1]/2;
		Indices[2] = Indicesraf[2]/2;
	
		int addWidthV = 0; //getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV; //Replace species.sparseBlockAddWidthV
		for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
	  	  for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
	    	for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
	      	  const vmesh::GlobalID globalID = vmesh->getGlobalID(Indices[0]+offset_vx,Indices[1]+offset_vy,Indices[2]+offset_vz);
	      	  ListBlockExist[popID].insert(globalID);
	    	}
	  	  }
		}
      }

    }
  }

  //create the cells if needed
  for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {

    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 ){

      vmesh::VelocityMesh* vmesh  = cell->get_velocity_mesh(popID);
      vmesh::LocalID Localsize= vmesh->size();
      
      for (vmesh::LocalID localID=0; localID<Localsize; ++localID) { //If the blocks don't need to exist anymore, they are removed
		vmesh::VelocityMesh* vmesh  = cell->get_velocity_mesh(popID); 
		vmesh::GlobalID globalID = vmesh->getGlobalID(localID);
		if (ListBlockExist[popID].find(globalID) == ListBlockExist[popID].end()) {
	  	  cell->remove_velocity_block(globalID,popID);
	      Localsize-=1;
	      localID-=1;	    
	    }
      }

      
      if(getObjectWrapper().particleSpecies[popID].RefinementLevel!=0){
	      
		for (vmesh::GlobalID globalID : ListBlockExist[popID]) {
	  
	  	  if(cell->add_velocity_block(globalID,popID)){ //True if it's a new block
	        //need to adapt the creation to the level of refinement needed
   	        vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
	 	    vmesh::VelocityMesh* vmeshgros= cell->get_velocity_mesh(popID-1);
		    Realf *datagros = cell->get_velocity_blocks(popID-1)->getData();
		    Realf *data     = cell->get_velocity_blocks(popID)->getData(); 
			  
	        vmesh::LocalID Indices[3];
	        vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);
			// i, j and k indicate the separation of the coarse block into eight refined blocks
	        int i = Indices[0]%2;
	        int j = Indices[1]%2;
	        int k = Indices[2]%2;
	    
	    	vmesh::LocalID Indicesgros[3];
	    	Indicesgros[0] = (Indices[0]-i)/2;
	    	Indicesgros[1] = (Indices[1]-j)/2;
	    	Indicesgros[2] = (Indices[2]-k)/2;

		    vmesh::GlobalID globalIDgros=vmeshgros->getGlobalID(Indicesgros);
		    vmesh::LocalID  localIDgros=vmeshgros->getLocalID(globalIDgros);
	    	vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
		    uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
		    ghost[0]=1;
	  
	    	for (int i2=0; i2<2; ++i2) {
	      	  for (int j2=0; j2<2; ++j2) {
				for (int k2=0; k2<2; ++k2) {
				  //Loop over the 8 coarsed cells located on the WID3 refined cells (the refined block)
		  		  data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		 		  data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2+1)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2+1)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2+1)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		  		  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2+1)]=datagros[localIDgros*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
				}
	      	  }
		    }
		
	      }else{
		    vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
		    vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
		    uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
		    ghost[0]=1;
		  }
		}
      }else{
	    for (vmesh::GlobalID globalID : ListBlockExist[popID]) {
	      cell->add_velocity_block(globalID,popID);
	      //ghost not important for this one
	      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
	      vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
	      uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
	      ghost[0]=1;
	    }
      }
      
    }
  }
}     
}

//Every vamr_refinedStep we check all the velocity cells with the vAMR criterion with the 3rd order
void RefinedOrder3(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells){

  Realf M[3][3][3];
  Realf A[3] = {1.0/8.0, 1.0, -1.0/8.0};

  for (int i0 = 0; i0 < 3; i0++) {
    for (int j0 = 0; j0 < 3; j0++) {
      for (int k0 = 0; k0 < 3; k0++) {
		M[i0][j0][k0] = A[i0] * A[j0] * A[k0];
      }
    }
  }

#pragma omp parallel for schedule(dynamic,1)
for (size_t c=0; c<cells.size(); ++c) {
  SpatialCell* cell = mpiGrid[cells[c]];
  std::unordered_set<vmesh::GlobalID> ListBlockExist[getObjectWrapper().particleSpecies.size()];
  
  for (int popID=(getObjectWrapper().particleSpecies.size()-2); popID>-1; --popID) {

    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 && getObjectWrapper().particleSpecies[popID].RefinementLevel<getObjectWrapper().particleSpecies[popID].MaxRefinementLevel){
       
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf    = cell->get_velocity_mesh(popID+1);
      Realf *data = cell->get_velocity_blocks(popID)->getData();
      
      for (vmesh::LocalID localID=0; localID<vmesh->size(); ++localID) { 
		vmesh::GlobalID globalID = vmesh->getGlobalID(localID);
		vmesh::LocalID Indices[3];
		vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);

		if (getObjectWrapper().particleSpecies[popID].RefinementLevel==0) {
	  	  ListBlockExist[popID].insert(globalID);
		}

       	bool Verif[8];
		for (int fill=0; fill<8; ++fill){
		  Verif[fill]=true;
		}
	  	//We will need 1 ghost cell of level R-1 in every direction, meaning that we need a 4x4x4 grid of level R-1 cells for our block of level R
		Realf Datagros[4][4][4];
	  
		for (int neighbour_vx=0; neighbour_vx<4; ++neighbour_vx) {
		  int n_vx=(neighbour_vx+1)/2 -1; //gives the shift in block
		  int i=abs((neighbour_vx+1)%2);
		  for (int neighbour_vy=0; neighbour_vy<4; ++neighbour_vy) {
		    int n_vy=(neighbour_vy+1)/2 -1;
		    int j=abs((neighbour_vy+1)%2);
		    for (int neighbour_vz=0; neighbour_vz<4; ++neighbour_vz) {
		      int n_vz=(neighbour_vz+1)/2 -1;
		      int k=abs((neighbour_vz+1)%2);
	      
		      vmesh::LocalID localIDneigh = vmesh->invalidLocalID();
	      
		      if(n_vx==0 && n_vy==0 && n_vz==0){
			  //The cell of level R-1 will be defined on the actual block of level R
			    localIDneigh = localID;
		      }else{
			  //The cell of level R-1 will be defined on a neighbouring block of level R
				const vmesh::GlobalID globalIDneigh = vmesh->getGlobalID(Indices[0]+n_vx,Indices[1]+n_vy,Indices[2]+n_vz);
			    if (globalIDneigh !=  vmesh->invalidGlobalID()) {
			      localIDneigh=vmesh->getLocalID(globalIDneigh);
			    }
		      }

		      if (localIDneigh !=  vmesh->invalidLocalID()) {
			  //The block exists, meaning that we can compute the cell value for level R-1
			    Realf Datagrossum=0;
			    for (int i2=0; i2<2; ++i2) {
			      for (int j2=0; j2<2; ++j2) {
			        for (int k2=0; k2<2; ++k2) {
			          Datagrossum += data[localIDneigh*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
			        }
			      }
			    }
			    //Computation of the grid R-1
			    Datagros[neighbour_vx][neighbour_vy][neighbour_vz] = Datagrossum/8.0;
		      }else{
			    for(int iloop = 0; iloop< (2-abs(n_vx)); ++iloop){
			      for(int jloop = 0; jloop< (2-abs(n_vy)); ++jloop){
			        for(int kloop = 0; kloop< (2-abs(n_vz)); ++kloop){
					//Loop over all the future refined blocks that will be impacted by the non-existing block
				    //It is technically impossible for the block in which we have the local ID to return invalidLocalID(), but this has been taken into account
			          Verif[(std::max(n_vx,0)+iloop)+2*(std::max(n_vy,0)+jloop)+4*(std::max(n_vz,0)+kloop)]=false;
			        }
			      }
			    }
		  
		      }
	      
		    }
		  }
		}
	
		
	    //Loop over the future refined blocks R+1
		for (int i=0; i<2; ++i) {
	  	  for (int j=0; j<2; ++j) {
		    for (int k=0; k<2; ++k) {
		      if (Verif[i+2*j+4*k]){
			    Realf Datagrossum = 0;
	
		        for (int i2=-1; i2<2; ++i2) {
			      for (int j2=-1; j2<2; ++j2) {
			        for (int k2=-1; k2<2; ++k2) {
					  //Sum over the 8 cells of level R in order to reproduce the coarse cell R-1
			          Datagrossum += M[1-i2*(1-2*i)][1-j2*(1-2*j)][1-k2*(1-2*k)]*Datagros[1+i+i2][1+j+j2][1+k+k2];
			        }
			      }
		        }
	    	      	     		    
		        Realf D = abs( data[localID*WID3+cellIndex(1+i,1+j,1+k)] - Datagrossum );
				// The idea is to always compare the 8 central cells of level R with the reproduce R-1 cells	      

	      	    vmesh::LocalID Indicesraf[3];
	      	    Indicesraf[0] = 2*Indices[0]+i ;
	      	    Indicesraf[1] = 2*Indices[1]+j ;
	      	    Indicesraf[2] = 2*Indices[2]+k ;

                  Realf Dcomp = cell->getVelocityBlockMinValue(popID);

                  if(getObjectWrapper().particleSpecies[popID].CriteriaMethod==1){ //epsilon
                    Dcomp= getObjectWrapper().particleSpecies[popID].CriteriaValue;
                  }else if(getObjectWrapper().particleSpecies[popID].CriteriaMethod==2){ //epsilon*2^-R
 	                Dcomp= getObjectWrapper().particleSpecies[popID].CriteriaValue/(1u << getObjectWrapper().particleSpecies[popID].RefinementLevel);
                  }
                  if (D > Dcomp){
				  // We should create a new block for R+1
			      int addWidthV = 1;//getObjectWrapper().particleSpecies[popID+1].sparseBlockAddWidthV; 
				  for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
		  		    for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
		    		  for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
		      		    const vmesh::GlobalID globalIDraf = vmeshraf->getGlobalID(Indicesraf[0]+offset_vx,Indicesraf[1]+offset_vy,Indicesraf[2]+offset_vz);
		      		    if (globalIDraf==  vmeshraf->invalidGlobalID()) {
						  // std::cout<< " GlobalID bug not normal" << "Indices[0]+offset_vx" << Indices[0]+offset_vx << "Indices[1]+offset_vy " << Indices[1]+offset_vy << "Indices[1]+offset_vy" << Indices[2]+offset_vz <<std::endl;
		      		    }else{
		      			  ListBlockExist[popID+1].insert(globalIDraf);
		    		    }
					  }
				   }
	     		 }
				  
			   }
				  
		      }
	       	}
	  	  }
		} 

      }
    }
  }

  //We will now ensure that each level R+1 exist on the level R with the proper ghost cells between both
  for (int popID=(getObjectWrapper().particleSpecies.size()-2); popID>-1; --popID) {
    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 && getObjectWrapper().particleSpecies[popID].RefinementLevel<getObjectWrapper().particleSpecies[popID].MaxRefinementLevel){

      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf  = cell->get_velocity_mesh(popID+1);
      
      for (vmesh::GlobalID globalIDraf : ListBlockExist[popID+1]) {
	
		vmesh::LocalID Indicesraf[3];
		vmeshraf->getIndices(globalIDraf, Indicesraf[0], Indicesraf[1], Indicesraf[2]);

		vmesh::LocalID Indices[3];
		Indices[0] = Indicesraf[0]/2;
		Indices[1] = Indicesraf[1]/2;
		Indices[2] = Indicesraf[2]/2;
	
		int addWidthV = 1;//getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
		for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
	 	  for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
	    	for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
	      	  const vmesh::GlobalID globalID = vmesh->getGlobalID(Indices[0]+offset_vx,Indices[1]+offset_vy,Indices[2]+offset_vz);
	      	  if (globalID==  vmesh->invalidGlobalID()) {
				std::cout<< " GlobalID bug not normal" << "Indices[0]+offset_vx" << Indices[0]+offset_vx << "Indices[1]+offset_vy " << Indices[1]+offset_vy << "Indices[1]+offset_vy" << Indices[2]+offset_vz <<std::endl;
	      	  }
	      	  ListBlockExist[popID].insert(globalID);
	    	}
	  	  }
		}
		  
      }
    }
  }

  //create the cells if needed
  for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0){

      vmesh::VelocityMesh* vmeshinit  = cell->get_velocity_mesh(popID);
      vmesh::LocalID Localsize= vmeshinit->size();
      
      for (vmesh::LocalID localID=0; localID<Localsize; ++localID) { //If the blocks don't need to exist anymore, they are removed
		vmesh::VelocityMesh* vmeshinit2  = cell->get_velocity_mesh(popID);
		vmesh::GlobalID globalID = vmeshinit2->getGlobalID(localID);
		if (ListBlockExist[popID].find(globalID) == ListBlockExist[popID].end()) {
	  	  cell->remove_velocity_block(globalID,popID);
	  	  Localsize-=1;
	  	  localID-=1;	    
		}
      }

      
      if(getObjectWrapper().particleSpecies[popID].RefinementLevel!=0){
		for (vmesh::GlobalID globalID : ListBlockExist[popID]) {
	 	  if(cell->add_velocity_block(globalID,popID)){
	    	//True if it's a new block
	    	//need to adapt the creation to the level of refinement needed
		    vmesh::VelocityMesh* vmesh  = cell->get_velocity_mesh(popID);
	    	vmesh::VelocityMesh* vmeshgros  = cell->get_velocity_mesh(popID-1);
	    	Realf *datagros = cell->get_velocity_blocks(popID-1)->getData();
	    	Realf *data = cell->get_velocity_blocks(popID)->getData(); 

	    	vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
	    	if (localIDcreated==  vmesh->invalidLocalID()) {
	      	  //std::cout<< " localIDcreated bug not normal"  <<std::endl;
	    	}
	    
	    	vmesh::LocalID Indices[3];
	    	vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);
			// i, j and k indicate the separation of the coarse block into eight refined blocks
	    	int i = Indices[0]%2;
	    	int j = Indices[1]%2;
	    	int k = Indices[2]%2;

	    	vmesh::LocalID Indicesgrosinit[3];
	    	Indicesgrosinit[0] = (Indices[0]-i)/2;
	    	Indicesgrosinit[1] = (Indices[1]-j)/2;
	    	Indicesgrosinit[2] = (Indices[2]-k)/2;

	    	for (int i2=0; i2<2; ++i2) {
	      	  for (int j2=0; j2<2; ++j2) {
				for (int k2=0; k2<2; ++k2) {
                //Loop over the 8 coarsed cells located on the WID3 refined cells (the refined block)
				// We ensure that the cells are well initialised to 0
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2+1)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2+1)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2+1)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2+1)]=0.0;
		  		  for (int i4=-1; i4<2; ++i4){
		    		for (int j4=-1; j4<2; ++j4){
		      		  for (int k4=-1; k4<2; ++k4){
					  // Loop over the coarse neighbour cell and the actual cell
				      // i5, j5 and k5 are the indices of the coarse cell that will be used for interpolation
				      //It can be negative or >WID-1 if we need to go to another block			      
						int i5=2*i+i2+i4;
						int j5=2*j+j2+j4;
						int k5=2*k+k2+k4;

						vmesh::LocalID Indicesgros[3];
						Indicesgros[0] = Indicesgrosinit[0];
						Indicesgros[1] = Indicesgrosinit[1];
						Indicesgros[2] = Indicesgrosinit[2];
					 	//We check whether we need to move to a different block.
						if(i5<0){
			  			  Indicesgros[0]-= 1;
						}else if(i5> WID-1){
			  			  Indicesgros[0]+= 1;
						}
			
						if(j5<0){
			  			  Indicesgros[1]-= 1;
						}else if(j5> WID-1){
			  			  Indicesgros[1]+= 1;
						}
			
						if(k5<0){
			  			  Indicesgros[2]-= 1;
						}else if(k5> WID-1){
			  			  Indicesgros[2]+= 1;
						}

						vmesh::GlobalID globalIDgros=vmeshgros->getGlobalID(Indicesgros);
						if (globalIDgros==  vmeshgros->invalidGlobalID()) {
			 			  //std::cout<< " GlobalID bug not normal"<< "Indices[0] " << Indices[0] << "Indicesgros[0] "<< Indicesgros[0] << " GlobalID " << globalID << " PopID " << popID <<std::endl;
						}
	     
						vmesh::LocalID  localIDgros=vmeshgros->getLocalID(globalIDgros);
						if (localIDgros==  vmeshgros->invalidLocalID()) {
			  			  //std::cout<< " localIDgros bug not normal"<< "Indices[0] " << Indices[0] << "Indicesgros[0] "<< Indicesgros[0] << " GlobalIDgros " << globalIDgros << " i5, j5, k5 " << i5 << " " << j5 << " " << k5 << " "  <<std::endl;
						}
						// We change the cell coordinates if we move to a different block
						i5 = (i5 + WID)%WID;
						j5 = (j5 + WID)%WID;
						k5 = (k5 + WID)%WID;

						if (i5<0 || i5>(WID-1) || j5<0 || j5>(WID-1) ||  k5<0 || k5>(WID-1) ) {
			  			  //std::cout<< " localidgros bug not normal"<< "Indices[0] " << Indices[0] << "Indicesgros[0] "<< Indicesgros[0] << " GlobalIDgros " << globalIDgros << " i5, j5, k5 " << i5 << " " << j5 << " " << k5 << " "  <<std::endl;
						}
			
						Realf contribution=datagros[localIDgros*WID3+cellIndex(i5,j5,k5)];

						data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2)]+=M[i4+1][j4+1][k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2)]+=M[-i4+1][j4+1][k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2)]+=M[i4+1][-j4+1][k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2+1)]+=M[i4+1][j4+1][-k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2)]+=M[-i4+1][-j4+1][k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2+1)]+=M[-i4+1][j4+1][-k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2+1)]+=M[i4+1][-j4+1][-k4+1]*contribution;
						data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2+1)]+=M[-i4+1][-j4+1][-k4+1]*contribution;
				       
		      		  }		
		    		}
				  }
				}
	      	  }
	    	} 
				  
		  }else{
		    vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
		    vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
		    uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
		    ghost[0]=1;
		  }
		}
      }else{
	    for (vmesh::GlobalID globalID : ListBlockExist[popID]) {
	      cell->add_velocity_block(globalID,popID);
	      //ghost not important for this one
	      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
	      vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
	      uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
	      ghost[0]=1;
	    }
      }
	          
    }
  }
}     
}

//Every vamr_refinedStep we check all the velocity cells with the vAMR criterion with the 5th order
void RefinedOrder5(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells){
  
  
  Realf M[5][5][5];
  Realf A[5] = {-3.0/128.0, 11.0/64.0, 1.0, -11.0/64.0, 3.0/128.0};
  
  for (int i0 = 0; i0 < 5; i0++) {
    for (int j0 = 0; j0 < 5; j0++) {
      for (int k0 = 0; k0 < 5; k0++) {
	    M[i0][j0][k0] = A[i0] * A[j0] * A[k0];
      }
    }
  }

  #pragma omp parallel for schedule(dynamic,1)
  for (size_t c=0; c<cells.size(); ++c) {
    SpatialCell* cell = mpiGrid[cells[c]];
    std::unordered_set<vmesh::GlobalID> ListBlockExist[getObjectWrapper().particleSpecies.size()];
  
    for (int popID=(getObjectWrapper().particleSpecies.size()-2); popID>-1; --popID) {
      if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 && getObjectWrapper().particleSpecies[popID].RefinementLevel<getObjectWrapper().particleSpecies[popID].MaxRefinementLevel){
       
        vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
        vmesh::VelocityMesh* vmeshraf    = cell->get_velocity_mesh(popID+1);
        Realf *data = cell->get_velocity_blocks(popID)->getData();
      
        for (vmesh::LocalID localID=0; localID<vmesh->size(); ++localID) { 
	      vmesh::GlobalID globalID = vmesh->getGlobalID(localID);
	      vmesh::LocalID Indices[3];
	      vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);

	      if (getObjectWrapper().particleSpecies[popID].RefinementLevel==0) {
	        ListBlockExist[popID].insert(globalID);
	      }

	      bool Verif=true;
		  //We will need 2 ghost cells of level R-1 in every direction, meaning that we need a 6x6x6 grid of level R-1 cells for our block of level R
	      Realf Datagros[6][6][6];
	  
	      for (int neighbour_vx=0; neighbour_vx<6 && Verif; ++neighbour_vx) {
	        int n_vx=neighbour_vx/2 -1; //gives the shift in block
	        int i=abs(neighbour_vx%2);
	        for (int neighbour_vy=0; neighbour_vy<6 && Verif; ++neighbour_vy) {
	          int n_vy=neighbour_vy/2 -1;
	          int j=abs(neighbour_vy%2);
	          for (int neighbour_vz=0; neighbour_vz<6 && Verif; ++neighbour_vz) {
	            int n_vz=neighbour_vz/2 -1;
	            int k=abs(neighbour_vz%2);
	      
	            vmesh::LocalID localIDneigh = vmesh->invalidLocalID();
	      
	            if(n_vx==0 && n_vy==0 && n_vz==0){
				//The cell of level R-1 will be defined on the actual block of level R
		          localIDneigh = localID;
	            }else{
				//The cell of level R-1 will be defined on a neighbouring block of level R
		          const vmesh::GlobalID globalIDneigh = vmesh->getGlobalID(Indices[0]+n_vx,Indices[1]+n_vy,Indices[2]+n_vz);
		          if (globalIDneigh !=  vmesh->invalidGlobalID()) {
		            localIDneigh=vmesh->getLocalID(globalIDneigh);
		          }
	            }

	            if (localIDneigh !=  vmesh->invalidLocalID()) {
				//The block exists, meaning that we can compute the cell value for level R-1
		          Realf Datagrossum=0;
		          for (int i2=0; i2<2; ++i2) {
		            for (int j2=0; j2<2; ++j2) {
		              for (int k2=0; k2<2; ++k2) {
		                Datagrossum += data[localIDneigh*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
		              }
		            }
		          }
				  //Computation of the grid R-1
		          Datagros[neighbour_vx][neighbour_vy][neighbour_vz] = Datagrossum/8.0;
	            }else{
		        //If one of the neighbouring blocks does not exist, it will directly impact all future R-1 blocks
		          Verif=false;
	            }
	      
	          }
	        }
	      }
	  
	      
	      if (Verif){
	        for (int i=0; i<2; ++i) {
	          for (int j=0; j<2; ++j) {
	            for (int k=0; k<2; ++k) {
				//Loop over the future refined blocks R+1
	              Realf Datagrossum = 0;
	
		          for (int i2=-2; i2<3; ++i2) {
		            for (int j2=-2; j2<3; ++j2) {
		              for (int k2=-2; k2<3; ++k2) {
					  //Sum over the 8 cells of level R in order to reproduce the coarse cell R-1
		                Datagrossum += M[2-i2*(1-2*i)][2-j2*(1-2*j)][2-k2*(1-2*k)]*Datagros[2+i+i2][2+j+j2][2+k+k2];
		              }
		            }
		          }

		          Realf D = abs( data[localID*WID3+cellIndex(1+i,1+j,1+k)] - Datagrossum );
				  // The idea is to always compare the 8 central cells of level R with the reproduce R-1 cells
	      	      vmesh::LocalID Indicesraf[3];
	      	      Indicesraf[0] = 2*Indices[0]+i ;
	      	      Indicesraf[1] = 2*Indices[1]+j ;
	      	      Indicesraf[2] = 2*Indices[2]+k ;

                  Realf Dcomp = cell->getVelocityBlockMinValue(popID);

                  if(getObjectWrapper().particleSpecies[popID].CriteriaMethod==1){ //epsilon
                    Dcomp= getObjectWrapper().particleSpecies[popID].CriteriaValue;
                  }else if(getObjectWrapper().particleSpecies[popID].CriteriaMethod==2){ //epsilon*2^-R
 	                Dcomp= getObjectWrapper().particleSpecies[popID].CriteriaValue/(1u << getObjectWrapper().particleSpecies[popID].RefinementLevel);
                  }
                  if (D > Dcomp){
				  // We should create a new block for R+1
				  int addWidthV = 1;//getObjectWrapper().particleSpecies[popID+1].sparseBlockAddWidthV; 
				  for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
		  		    for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
		    		  for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
		      		    const vmesh::GlobalID globalIDraf = vmeshraf->getGlobalID(Indicesraf[0]+offset_vx,Indicesraf[1]+offset_vy,Indicesraf[2]+offset_vz);
		      		    if (globalIDraf==  vmeshraf->invalidGlobalID()) {
						  // std::cout<< " GlobalID bug not normal" << "Indices[0]+offset_vx" << Indices[0]+offset_vx << "Indices[1]+offset_vy " << Indices[1]+offset_vy << "Indices[1]+offset_vy" << Indices[2]+offset_vz <<std::endl;
		      		    }else{
		      			  ListBlockExist[popID+1].insert(globalIDraf);
		    		    }
						
		  			  }
				    } 
	     		  }
				  
			    }
					
	          }
	    	}
	  	  }
		} 

      }
    }
  }

  //We will now ensure that each level R+1 exist on the level R with the proper ghost cells between both
  for (int popID=(getObjectWrapper().particleSpecies.size()-2); popID>-1; --popID) {
    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0 && getObjectWrapper().particleSpecies[popID].RefinementLevel<getObjectWrapper().particleSpecies[popID].MaxRefinementLevel){

      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      vmesh::VelocityMesh* vmeshraf    = cell->get_velocity_mesh(popID+1);
      
      for (vmesh::GlobalID globalIDraf : ListBlockExist[popID+1]) {
	
		vmesh::LocalID Indicesraf[3];
		vmeshraf->getIndices(globalIDraf, Indicesraf[0], Indicesraf[1], Indicesraf[2]);

		vmesh::LocalID Indices[3];
		Indices[0] = Indicesraf[0]/2;
		Indices[1] = Indicesraf[1]/2;
		Indices[2] = Indicesraf[2]/2;
	
		int addWidthV = 1;//getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV; //Replace species.sparseBlockAddWidthV
		  for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
	  		for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
	    	  for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
	      		const vmesh::GlobalID globalID = vmesh->getGlobalID(Indices[0]+offset_vx,Indices[1]+offset_vy,Indices[2]+offset_vz);
	      		if (globalID==  vmesh->invalidGlobalID()) {
				  std::cout<< " GlobalID bug not normal" << "Indices[0]+offset_vx" << Indices[0]+offset_vx << "Indices[1]+offset_vy " << Indices[1]+offset_vy << "Indices[1]+offset_vy" << Indices[2]+offset_vz <<std::endl;
	      		}
	      		ListBlockExist[popID].insert(globalID);
	    	  }
	  		}
		  }
        }

    }
  }

  //create the cells if needed
  for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
    if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0){

      vmesh::VelocityMesh* vmeshinit  = cell->get_velocity_mesh(popID);
      vmesh::LocalID Localsize= vmeshinit->size();
      
      for (vmesh::LocalID localID=0; localID<Localsize; ++localID) { //If the blocks don't need to exist anymore, they are removed
		vmesh::VelocityMesh* vmeshinit2  = cell->get_velocity_mesh(popID);
		vmesh::GlobalID globalID = vmeshinit2->getGlobalID(localID);
		if (ListBlockExist[popID].find(globalID) == ListBlockExist[popID].end()) {
	  	  cell->remove_velocity_block(globalID,popID);
	  	  Localsize-=1;
	  	  localID-=1;	    
		}
      }
      
      if(getObjectWrapper().particleSpecies[popID].RefinementLevel!=0){
		for (vmesh::GlobalID globalID : ListBlockExist[popID]) {
	  	  if(cell->add_velocity_block(globalID,popID)){
	        //True if it's a new block
	        //need to adapt the creation to the level of refinement needed
	        vmesh::VelocityMesh* vmesh  = cell->get_velocity_mesh(popID);
	        vmesh::VelocityMesh* vmeshgros  = cell->get_velocity_mesh(popID-1);
	        Realf *datagros = cell->get_velocity_blocks(popID-1)->getData();
	        Realf *data = cell->get_velocity_blocks(popID)->getData(); 

	        vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
	        if (localIDcreated==  vmesh->invalidLocalID()) {
	    	  std::cout<< " localIDcreated bug not normal"  <<std::endl;
	        }
	    
	        vmesh::LocalID Indices[3];
	        vmesh->getIndices(globalID, Indices[0], Indices[1], Indices[2]);
			// i, j and k indicate the separation of the coarse block into eight refined blocks
	        int i = Indices[0]%2;
	        int j = Indices[1]%2;
	        int k = Indices[2]%2;

	        vmesh::LocalID Indicesgrosinit[3];
	        Indicesgrosinit[0] = (Indices[0]-i)/2;
	        Indicesgrosinit[1] = (Indices[1]-j)/2;
	        Indicesgrosinit[2] = (Indices[2]-k)/2;

	        for (int i2=0; i2<2; ++i2) {
	          for (int j2=0; j2<2; ++j2) {
			    for (int k2=0; k2<2; ++k2) {
                //Loop over the 8 coarsed cells located on the WID3 refined cells (the refined block)
				// We ensure that the cells are well initialised to 0
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2+1)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2+1)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2+1)]=0.0;
				  data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2+1)]=0.0;
		  		  for (int i4=-2; i4<3; ++i4){
		    	    for (int j4=-2; j4<3; ++j4){
		      		  for (int k4=-2; k4<3; ++k4){
			          // Loop over the coarse neighbour cell and the actual cell
				      // i5, j5 and k5 are the indices of the coarse cell that will be used for interpolation
				      //It can be negative or >WID-1 if we need to go to another block
					    int i5=2*i+i2+i4;
					    int j5=2*j+j2+j4;
					    int k5=2*k+k2+k4;

					    vmesh::LocalID Indicesgros[3];
					    Indicesgros[0] = Indicesgrosinit[0];
					    Indicesgros[1] = Indicesgrosinit[1];
					    Indicesgros[2] = Indicesgrosinit[2];
			      		//We check whether we need to move to a different block
					    if(i5<0){
			  			  Indicesgros[0]-= 1;
					    }else if(i5> WID-1){
			  			  Indicesgros[0]+= 1;
					    }
			
					    if(j5<0){
			  			  Indicesgros[1]-= 1;
					    }else if(j5> WID-1){
			  			  Indicesgros[1]+= 1;
					    }
			
					    if(k5<0){
			  			  Indicesgros[2]-= 1;
					    }else if(k5> WID-1){
			  			  Indicesgros[2]+= 1;
					    }

					    vmesh::GlobalID globalIDgros=vmeshgros->getGlobalID(Indicesgros);
					    if (globalIDgros==  vmeshgros->invalidGlobalID()) {
			  			  std::cout<< " GlobalID bug not normal"<< "Indices[0] " << Indices[0] << "Indicesgros[0] "<< Indicesgros[0] << " GlobalID " << globalID << " PopID " << popID <<std::endl;
					    }
	     
					    vmesh::LocalID  localIDgros=vmeshgros->getLocalID(globalIDgros);
					    if (localIDgros==  vmeshgros->invalidLocalID()) {
			  			  std::cout<< " localIDgros bug not normal"<< "Indices[0] " << Indices[0] << "Indicesgros[0] "<< Indicesgros[0] << " GlobalIDgros " << globalIDgros << " i5, j5, k5 " << i5 << " " << j5 << " " << k5 << " "  <<std::endl;
					    }
						// We change the cell coordinates if we move to a different block
					    i5 = (i5 + WID)%WID;
					    j5 = (j5 + WID)%WID;
		 			    k5 = (k5 + WID)%WID;
			
					    Realf contribution=datagros[localIDgros*WID3+cellIndex(i5,j5,k5)];

						data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2)]+=M[i4+2][j4+2][k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2)]+=M[-i4+2][j4+2][k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2)]+=M[i4+2][-j4+2][k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2,2*j2,2*k2+1)]+=M[i4+2][j4+2][-k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2)]+=M[-i4+2][-j4+2][k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2,2*k2+1)]+=M[-i4+2][j4+2][-k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2,2*j2+1,2*k2+1)]+=M[i4+2][-j4+2][-k4+2]*contribution;
					    data[localIDcreated*WID3+cellIndex(2*i2+1,2*j2+1,2*k2+1)]+=M[-i4+2][-j4+2][-k4+2]*contribution;
				       
		      		  }		
		    	    }
		  		  }
		  
			    }
	  		  }
		    } 

  	      }else{
		    vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
		    vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
		    uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
		    ghost[0]=1;
		  }
		}
      }else{
	    for (vmesh::GlobalID globalID : ListBlockExist[popID]) {
	      cell->add_velocity_block(globalID,popID);
	      //ghost not important for this one
	      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
	      vmesh::LocalID  localIDcreated=vmesh->getLocalID(globalID);
	      uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost(localIDcreated);
	      ghost[0]=1;
	    }
      }
	             
    }
  }
}     
}

//Update the ghost parameter and remove the newly created cells that don't respect the criterion
void SmallRefinedOrder1(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells){

#pragma omp parallel for schedule(dynamic,1)
for (size_t c=0; c<cells.size(); ++c) {
SpatialCell* cell = mpiGrid[cells[c]];

 for (uint popID=1; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
   if(getObjectWrapper().particleSpecies[popID].MaxRefinementLevel>0){
     vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
     Realf *data = cell->get_velocity_blocks(popID)->getData();
     uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost();

     vmesh::LocalID Localsize= vmesh->size();
     for (vmesh::LocalID localID=0; localID<Localsize; ++localID) {
       if(ghost[localID]==0){ 
	     // Here the idea will be to check if we keep the newly created blocks of level R+1, i.e. with ghost=0
	     Realf Datagrosgros = 0; //R-1 cell of the same size of the block of level R+1
	     Realf Datagros = 0; //R cell in the middle of the block of level R , this cell is not consistent with a real cell of level R but here it's only used for estimation
	   
	     for (int i=0; i<WID3; ++i) {
	       Datagrosgros += data[localID*WID3+i];
	     }		    	  	  
	     Datagrosgros/=WID3;
		 // i2, j2, k2 are not the same as always because we here estimate a non-existing central cell with level R+1
	     for (int i2=0; i2<2; ++i2) {
	       for (int j2=0; j2<2; ++j2) {
	         for (int k2=0; k2<2; ++k2) {
	           Datagros += data[localID*WID3+cellIndex(1+i2,1+j2,1+k2)];
	         }
	       }
	     }   	  
	     Datagros/=8.0;	      	     		    
	     Realf D = abs( Datagrosgros - Datagros );
		 // The idea is to always compare the central cell of level R with the reproduce R-1 cell	  
         Realf Dcomp = cell->getVelocityBlockMinValue(popID-1);

         if(getObjectWrapper().particleSpecies[popID-1].CriteriaMethod==1){ //epsilon
     	   Dcomp= getObjectWrapper().particleSpecies[popID-1].CriteriaValue;
         }else if(getObjectWrapper().particleSpecies[popID-1].CriteriaMethod==2){ //epsilon*2^-R
           Dcomp= getObjectWrapper().particleSpecies[popID-1].CriteriaValue/(1u << getObjectWrapper().particleSpecies[popID-1].RefinementLevel);
         }
         if (D > Dcomp){
	       // We keep the block
	       ghost[localID]=1; // This could be changed to 2 to distinguish the newly created cells from those that already exist
	     }else{
	       //If the block don't need to exist anymore, it is removed
	       vmesh::VelocityMesh* vmesh  = cell->get_velocity_mesh(popID); 
	       vmesh::GlobalID globalID = vmesh->getGlobalID(localID);
	       cell->remove_velocity_block(globalID,popID);
	       Localsize-=1;
	       localID-=1;
	     }
       }	      	
     }
	   
   }
  }
}
}

//Fixed ghost to 1 for all the initial velocity cells to avoid destruction
void GhostFixation(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells){

#pragma omp parallel for schedule(dynamic,1)
for (size_t c=0; c<cells.size(); ++c) {
  SpatialCell* cell = mpiGrid[cells[c]];
  for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
      uint8_t *ghost = cell->get_velocity_blocks(popID)->getGhost();     
      for (vmesh::LocalID localID=0; localID<vmesh->size(); ++localID) {
	  ghost[localID]=1;
    }	  
  }
}    
}
