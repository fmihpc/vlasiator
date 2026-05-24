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
	      changeRefined(cell,popID);
	    }
	    // Calculate species' contribution to first velocity moments with Vamr
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

      //      std::cout<< " popID  "<< popID << " a pour RHO  "<< pop.RHO <<std::endl;
      if(P::activateVamr && Ref==MaxRef && MaxRef>0 ){
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
	         changeRefined(cell,popID);
	       }
	       // Calculate species' contribution to first velocity moments with Vamr
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

	 // Temporary array for storing moments
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
	         changeRefined(cell,popID);
	       }
	       // Calculate species' contribution to first velocity moments with Vamr
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

	     // Temporary array for storing moments
	     Real array[nMom2] = {0};
      
	     int Ref=0;
	     int MaxRef=0;
	     if(P::activateVamr) {
	       int Ref=getObjectWrapper().particleSpecies[popID].RefinementLevel;
	       int MaxRef=getObjectWrapper().particleSpecies[popID].MaxRefinementLevel;
	       // Calculate species' contribution to second velocity moments with Vamr
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


void vamr_transfer_values(
			  dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			  const std::vector<CellID>& cells,
			  const uint popID) {

  phiprof::Timer computeMomentsTimer {"vamr_transfer_value"};
  // Loop over all particle species
#pragma omp parallel for schedule(dynamic,1)
  for (size_t c=0; c<cells.size(); ++c) {
    SpatialCell* cell = mpiGrid[cells[c]];

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
	        	      Realf summ=0;
			          Realf datasave= data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)];
			          data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]=0;
			          for (int i3=0; i3<2; ++i3) {
			            for (int j3=0; j3<2; ++j3) {
			              for (int k3=0; k3<2; ++k3) {
			                if(dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]>0){
				              data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]+= dataraf[localIDraf*WID3+cellIndex(2*i2+i3,2*j2+j3,2*k2+k3)]/8.0;
				              summ+=1.0;
			                }
			              }
			            }
			          }
			          if (summ!=8.0){
			            data[localID*WID3+cellIndex(2*i+i2,2*j+j2,2*k+k2)]=datasave;
			          }
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


 void changeRefined(spatial_cell::SpatialCell* cell,
 const uint popID){

     vmesh::VelocityMesh* vmesh    = cell->get_velocity_mesh(popID);
     vmesh::VelocityMesh* vmeshraf = cell->get_velocity_mesh(popID+1);
     //
     uint8_t *refined =cell->get_velocity_blocks(popID)->getRefined();

     int dedans=0;
     int dehors=0;
     
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
