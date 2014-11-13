/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute
*/

#ifndef CPU_ACC_MAP_H
#define CPU_ACC_MAP_H

#include  "vec4.h"
#include "algorithm"
#include "cmath"
#include "utility"
#include "common.h"
#include "spatial_cell.hpp"
#include "cpu_acc_sort_blocks.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#define MAX_BLOCKS_PER_DIM 100

//index in the temporary and padded column data values array. Each
//column has an empty block in ether end.
#define i_pcolumnv(num_k_blocks, k_block, j, k) ( (j) * WID * (num_k_blocks + 2) +  (k) + ( k_block + 1 ) * WID )

using namespace spatial_cell;

struct PropagParams {
   Real intersection;
   Real intersection_di;
   Real intersection_dj;
   Real intersection_dk;
   
   Real dv;
   Real v_min;
   Real inv_dv;

   int dimension;
   int i_mapped;
   int j_mapped;
   int k_mapped;
   int refMul;
   uint8_t refLevel;
   uint8_t maxRefLevel;
   vmesh::GlobalID Nx;
   vmesh::GlobalID Ny;
   vmesh::GlobalID k_cell_global_target_max;
};

void map_1d(SpatialCell* spatial_cell,PropagParams& params,
	    vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
	    vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer);

void generateTargetMesh(SpatialCell* spatial_cell,const std::vector<vmesh::LocalID>& blocks,PropagParams& params,
			const uint8_t& targetRefLevel,const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh);

/* 
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt)

   TODO: parallelize with openMP over block-columns. If one also
   pre-creates new blocks in a separate loop first (serial operation),
   then the openmp parallization would scale well (better than over
   spatial cells), and would not need synchronization.
   
*/

bool map_1d(SpatialCell* spatial_cell,Transform<Real,3,Affine>& fwd_transform,Transform<Real,3,Affine>& bwd_transform,int dimension,int propag) {
   // Move the old velocity mesh and data to the variables below (very fast)
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> vmesh;
   vmesh::VelocityBlockContainer<vmesh::LocalID> blockContainer;
   spatial_cell->swap(vmesh,blockContainer);

   // Sort the blocks according to their refinement levels (very fast)
   const uint8_t maxRefLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMaxAllowedRefinementLevel();
   std::vector<std::vector<vmesh::LocalID> > blocks(maxRefLevel+1);
   for (vmesh::LocalID block=0; block<vmesh.size(); ++block) {
      uint8_t refLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getRefinementLevel(vmesh.getGlobalID(block));
      blocks[refLevel].push_back(block);
   }

   // Computer intersections etc.
   PropagParams propagParams;
   propagParams.maxRefLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMaxAllowedRefinementLevel();
   propagParams.Nx = SpatialCell::get_velocity_grid_length(propagParams.maxRefLevel)[0];
   propagParams.Ny = SpatialCell::get_velocity_grid_length(propagParams.maxRefLevel)[1];
   propagParams.dimension = dimension;
   switch (propag) {
    case 0:
      compute_intersections_1st(spatial_cell, bwd_transform, fwd_transform, dimension, propagParams.maxRefLevel,
				propagParams.intersection,propagParams.intersection_di,propagParams.intersection_dj,propagParams.intersection_dk);
      break;
    case 1:
      compute_intersections_2nd(spatial_cell, bwd_transform, fwd_transform, dimension, propagParams.maxRefLevel,
				propagParams.intersection,propagParams.intersection_di,propagParams.intersection_dj,propagParams.intersection_dk);
      break;
    case 2:
      compute_intersections_3rd(spatial_cell, bwd_transform, fwd_transform, dimension, propagParams.maxRefLevel,
				propagParams.intersection,propagParams.intersection_di,propagParams.intersection_dj,propagParams.intersection_dk);
      break;
    default:
      std::cerr << "error in map1d" << std::endl; exit(1);
      break;
   }

   propagParams.dimension = dimension;
   propagParams.dv    = SpatialCell::get_velocity_grid_cell_size(propagParams.maxRefLevel)[dimension];
   propagParams.v_min = SpatialCell::get_velocity_grid_min_limits()[dimension];
   propagParams.inv_dv = 1.0/propagParams.dv;
   propagParams.k_cell_global_target_max = SpatialCell::get_velocity_grid_length(propagParams.maxRefLevel)[dimension]*WID;
   switch (dimension) {
    case 0: {
      propagParams.i_mapped = 2;
      propagParams.j_mapped = 1;
      propagParams.k_mapped = 0;
      
      // Swap intersection di and dk
      const Real tmp=propagParams.intersection_di;
      propagParams.intersection_di=propagParams.intersection_dk;
      propagParams.intersection_dk=tmp;}
      break;
    case 1: {
      propagParams.i_mapped = 0;
      propagParams.j_mapped = 2;
      propagParams.k_mapped = 1;
      
      // Swap intersection dj and dk
      const Real tmp=propagParams.intersection_dj;
      propagParams.intersection_dj=propagParams.intersection_dk;
      propagParams.intersection_dk=tmp;}
      break;
    case 2:
      propagParams.i_mapped = 0;
      propagParams.j_mapped = 1;
      propagParams.k_mapped = 2;
      break;
    default:
      std::cerr << "error in map1d" << std::endl; exit(1);
      break;
   }
      
   bool rvalue = true;

   phiprof::start("mesh generation");
   for (uint8_t r=0; r<blocks.size(); ++r) {
      for (uint8_t rr=r; rr<blocks.size(); ++rr) {
	 propagParams.refLevel = rr;
	 generateTargetMesh(spatial_cell,blocks[rr],propagParams,r,vmesh);
      }
   }
   phiprof::stop("mesh generation");

   phiprof::start("mapping");
   map_1d(spatial_cell,propagParams,vmesh,blockContainer);
   phiprof::stop("mapping");

   // Merge values from coarse blocks to refined blocks wherever the same domain 
   // is covered by overlapping blocks (at different refinement levels)
   //spatial_cell->merge_values();
   
   //if (spatial_cell->checkMesh() == false) {
      //std::cerr << "error(s) in mesh, exiting" << std::endl; exit(1);
   //}

   rvalue = false;
}

void generateTargetMesh(SpatialCell* spatial_cell,const std::vector<vmesh::LocalID>& blocks,PropagParams& params,
			const uint8_t& targetRefLevel,const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh) {
   params.refMul = pow(2,(params.maxRefLevel-params.refLevel));
   const int baseMul = pow(2,params.maxRefLevel-targetRefLevel);

   for (size_t b=0; b<blocks.size(); ++b) {
      const vmesh::LocalID blockLID = blocks[b];
      const vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);
      vmesh::LocalID sourceIndex[3];
      vmesh.getIndices(blockGID,params.refLevel,sourceIndex[0],sourceIndex[1],sourceIndex[2]);
      
      switch (params.dimension) {
       case 0: {
	  const vmesh::LocalID temp=sourceIndex[2];
	  sourceIndex[2]=sourceIndex[0];
	  sourceIndex[0]=temp;}
	 break;
       case 1: {
	  vmesh::LocalID temp=sourceIndex[2];
	  sourceIndex[2]=sourceIndex[1];
	  sourceIndex[1]=temp;}
	 break;
       case 2:
	 break;
      }
      sourceIndex[0] *= WID*params.refMul;
      sourceIndex[1] *= WID*params.refMul;
      sourceIndex[2] *= WID*params.refMul;

      vmesh::GlobalID k_trgt_min = std::numeric_limits<vmesh::GlobalID>::max();
      vmesh::GlobalID k_trgt_max = 0;
      for (int k=0; k<2; ++k) for (int j=0; j<2; ++j) for (int i=0; i<2; ++i) {
	 const int i_cell = i*WID*params.refMul;
	 const int j_cell = j*WID*params.refMul;
	 const int k_cell = k*WID*params.refMul;

	 Realf intersection_min_base = params.intersection +
	   sourceIndex[0]*params.intersection_di +
	   (sourceIndex[1]+j_cell)*params.intersection_dj;

	 Realf intersection_min = intersection_min_base + i_cell*params.intersection_di;
	 Realf v_top = params.v_min + (sourceIndex[2]+k_cell)*params.dv;
	 vmesh::GlobalID lagrangian_gk_top = static_cast<vmesh::GlobalID>((v_top-intersection_min)/params.intersection_dk);

	 if (lagrangian_gk_top < k_trgt_min) k_trgt_min = lagrangian_gk_top;
	 if (lagrangian_gk_top > k_trgt_max) k_trgt_max = lagrangian_gk_top;
      }

      vmesh::GlobalID targetBlockIndex[3];
      //targetBlockIndex[params.i_mapped] = sourceIndex[params.i_mapped] / (WID*baseMul);
      //targetBlockIndex[params.j_mapped] = sourceIndex[params.j_mapped] / (WID*baseMul);
      targetBlockIndex[params.i_mapped] = sourceIndex[0] / (WID*baseMul);
      targetBlockIndex[params.j_mapped] = sourceIndex[1] / (WID*baseMul);
      k_trgt_min /= (WID*baseMul);
      k_trgt_max /= (WID*baseMul);
      for (vmesh::GlobalID k=k_trgt_min; k<=k_trgt_max; ++k) {
	 targetBlockIndex[params.k_mapped] = k;
	 vmesh::GlobalID targetBlock = vmesh.getGlobalID(targetRefLevel,targetBlockIndex);

	 if (targetRefLevel == 0) {
	    spatial_cell->add_velocity_block(targetBlock);
	 } else {
	    targetBlock = vmesh.getParent(targetBlock);
	    std::map<vmesh::GlobalID,vmesh::LocalID> insertedBlocks;
	    spatial_cell->refine_block(targetBlock,insertedBlocks);
	 }
      }
   }
}

/**
 * @param spatial_cell Propagated spatial cell, contains a valid (unsorted) target mesh.
 * @param params Parameters needed in the propagation.
 * @param vmesh The source mesh.
 * @param blockContainer Source mesh data.*/
void map_1d(SpatialCell* spatial_cell,PropagParams& params,
	   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
	   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer) {

   std::vector<vmesh::GlobalID> removeList;
   
   for (vmesh::LocalID targetLID=0; targetLID<spatial_cell->get_number_of_velocity_blocks(); ++targetLID) {
      vmesh::GlobalID targetGID = spatial_cell->get_velocity_block_global_id(targetLID);
      params.refLevel = vmesh.getRefinementLevel(targetGID);
      params.refMul = std::pow(2,params.maxRefLevel-params.refLevel);
      vmesh::LocalID targetIndex[3];
      vmesh.getIndices(targetGID,params.refLevel,targetIndex[0],targetIndex[1],targetIndex[2]);

      // Pointer to target data
      Realf* data_trgt = spatial_cell->get_data(targetLID);
      
      switch (params.dimension) {
       case 0: {
	  const vmesh::LocalID temp=targetIndex[2];
	  targetIndex[2]=targetIndex[0];
	  targetIndex[0]=temp;
       }
	 break;
       case 1: {
	  vmesh::LocalID temp=targetIndex[2];
	  targetIndex[2]=targetIndex[1];
	  targetIndex[1]=temp;
       }
	 break;
       case 2:
	 break;
      }
      targetIndex[0] *= WID*params.refMul;
      targetIndex[1] *= WID*params.refMul;
      targetIndex[2] *= WID*params.refMul;

      Real accum=0;

      // Note: the i,j,k indices below are transposed indices:
      for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
	 // Target cell global i/j/k indices (transposed):
	 vmesh::GlobalID i_cell_trgt = targetIndex[0] + i*params.refMul;
	 vmesh::GlobalID j_cell_trgt = targetIndex[1] + j*params.refMul;
	 vmesh::GlobalID k_cell_trgt = targetIndex[2] + k*params.refMul;

	 Realf intersection_min = params.intersection + i_cell_trgt*params.intersection_di + j_cell_trgt*params.intersection_dj;

	 // Velocities where the bottom and top z-faces of the target cell maps to:
	 Realf v_src_bot = k_cell_trgt*params.intersection_dk + intersection_min;
	 Realf v_src_top = (k_cell_trgt+params.refMul)*params.intersection_dk + intersection_min;

	 Realf v_src_bot_l = (v_src_bot-params.v_min)/params.dv;
	 Realf v_src_top_l = (v_src_top-params.v_min)/params.dv;

	 // Source cell global k-index:
	 vmesh::GlobalID k_cell_src_bot = static_cast<vmesh::GlobalID>((v_src_bot-params.v_min)/params.dv);
	 vmesh::GlobalID k_cell_src_top = static_cast<vmesh::GlobalID>((v_src_top-params.v_min)/params.dv);
	 if (k_cell_src_top < k_cell_src_bot) {
	    vmesh::GlobalID tmp=k_cell_src_bot;
	    k_cell_src_bot = k_cell_src_top;
	    k_cell_src_top = tmp;
	 }

	 // Source cell global i/j/k indices (untransposed):
	 vmesh::GlobalID srcIndex[3];
	 srcIndex[params.i_mapped] = i_cell_trgt;
	 srcIndex[params.j_mapped] = j_cell_trgt;

	 vmesh::GlobalID k_cell_src = k_cell_src_bot;
	 while (k_cell_src <= k_cell_src_top) {
	    // Find the source block
	    srcIndex[params.k_mapped] = k_cell_src;
	    uint8_t srcRefLevel = params.refLevel;
	    vmesh::GlobalID sourceGID = vmesh.findBlockDown(srcRefLevel,srcIndex);
	    int srcRefMul = std::pow(2,params.maxRefLevel-srcRefLevel);
	    
	    if (sourceGID == vmesh.invalidGlobalID()) {
	       ++k_cell_src; 
	       continue;
	    }
	    
	    Realf* data_src = blockContainer.getData(vmesh.getLocalID(sourceGID));

	    vmesh::GlobalID i_block_src = srcIndex[params.i_mapped] / (WID*srcRefMul);
	    vmesh::GlobalID j_block_src = srcIndex[params.j_mapped] / (WID*srcRefMul);
	    vmesh::GlobalID k_block_src = k_cell_src/(WID*srcRefMul);
	    int k_cell_bot = (k_cell_src - k_block_src*(WID*srcRefMul)) / srcRefMul;

	    int srcCellIndex[3];
	    srcCellIndex[params.i_mapped] = (i_cell_trgt - i_block_src*(WID*srcRefMul)) / srcRefMul;
	    srcCellIndex[params.j_mapped] = (j_cell_trgt - j_block_src*(WID*srcRefMul)) / srcRefMul;
	    srcCellIndex[params.k_mapped] = k_cell_bot;

	    Realf v_bot = std::max(v_src_bot_l,(Realf)(k_block_src*(WID*srcRefMul)+(k_cell_bot  )*srcRefMul));
	    Realf v_top = std::min(v_src_top_l,(Realf)(k_block_src*(WID*srcRefMul)+(k_cell_bot+1)*srcRefMul));

	    Realf v_int_min = (v_bot - (k_block_src*(WID*srcRefMul)+k_cell_bot*srcRefMul)) / srcRefMul;
	    Realf v_int_max = (v_top - (k_block_src*(WID*srcRefMul)+k_cell_bot*srcRefMul)) / srcRefMul;

	    Realf factor = std::pow(2,params.refLevel-srcRefLevel);
	    
	    int trgtCellIndex[3];
	    trgtCellIndex[params.i_mapped] = i;
	    trgtCellIndex[params.j_mapped] = j;
	    trgtCellIndex[params.k_mapped] = k;
	    const int trgtCell = vblock::index(trgtCellIndex[0],trgtCellIndex[1],trgtCellIndex[2]);
	    
	    data_trgt[trgtCell] += data_src[vblock::index(srcCellIndex[0],srcCellIndex[1],srcCellIndex[2])] * (v_int_max-v_int_min) * factor;
	    accum += data_src[vblock::index(srcCellIndex[0],srcCellIndex[1],srcCellIndex[2])] * (v_int_max-v_int_min) * factor;
	    k_cell_src = k_block_src*(WID*srcRefMul)+(k_cell_bot+1)*srcRefMul;
	 }
      }

      if (accum < 1e-30) {
	 removeList.push_back(targetGID);
      }
   }

   for (size_t b=0; b<removeList.size(); ++b) spatial_cell->remove_velocity_block(removeList[b]);
}

#endif   
