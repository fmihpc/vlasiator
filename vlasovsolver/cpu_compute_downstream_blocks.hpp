#ifndef CPU_DOWNSTREAM_BLOCKS_H
#define CPU_DOWNSTREAM_BLOCKS_H

#include "spatial_cell.hpp"
#include <Eigen/Geometry>
#include "../parameters.h"

using namespace Eigen;

inline void increment_cell_data(SpatialCell* spatial_cell,
                           const unsigned int fcell_i, const unsigned int fcell_j,const unsigned int fcell_k,
                           const Real value){
   

   const unsigned int block_i=fcell_i/WID;
   const unsigned int block_j=fcell_j/WID;
   const unsigned int block_k=fcell_k/WID;
   const unsigned int cell_i=fcell_i-block_i*WID;
   const unsigned int cell_j=fcell_j-block_j*WID;
   const unsigned int cell_k=fcell_k-block_k*WID;
   const unsigned int block = block_i + block_j * SpatialCell::vx_length + block_k * SpatialCell::vx_length * SpatialCell::vy_length;
   const unsigned int cell  = cell_i + cell_j * WID + cell_k * WID2;
   
   if(block_i>= SpatialCell::vx_length ||
      block_j>= SpatialCell::vy_length ||
      block_k>= SpatialCell::vz_length ||
      cell>64 ) {
     cout<< "warning, outside bounds"<<endl;
   //outside outer bounds -> cannot mark
      return;
   }        

   //add block if it does not exist yet
   
   if(spatial_cell->add_velocity_block(block)) {
     //block existed, or was created successfully
      Velocity_Block* block_ptr = spatial_cell->at(block);
      block_ptr->data[cell]+=value;
   }
   
}



/*
  Mark downstream cells for a particular velocity-cell
  
  \par spatial_cell Spatial cell that is integrated
  \par v velocity lower left corner of velocity-cell after rotation
*/

inline void mark_downstream_cells(SpatialCell* spatial_cell,const Array3d v){

   //note that fcell_ijk can also have a value of -1 right at the edge(!)
   //When we compute fcell_i we compare inside which shifted-by-halfvxyz
   //cube the velocity is. This gives us information about the eight
   //potential targets to which we will add density. fcell_ijk is the
   //index of the lower left one
   const int fcell_i=(v[0] - SpatialCell::vx_min-0.5*SpatialCell::cell_dvx) / SpatialCell::cell_dvx;
   const int fcell_j=(v[1] - SpatialCell::vy_min-0.5*SpatialCell::cell_dvy) / SpatialCell::cell_dvy;
   const int fcell_k=(v[2] - SpatialCell::vz_min-0.5*SpatialCell::cell_dvz) / SpatialCell::cell_dvz;
   const unsigned int fcell_p1_i=fcell_i+1;
   const unsigned int fcell_p1_j=fcell_j+1;
   const unsigned int fcell_p1_k=fcell_k+1;

   if (fcell_i < 0 || fcell_j < 0 || fcell_k < 0 ) {
      //outside of out velocity space edge (!)
      return;
   }
   
   //mark all target cells with 1.0 in data
   increment_cell_data(spatial_cell, fcell_i   , fcell_j   , fcell_k   , 1.0);
   increment_cell_data(spatial_cell, fcell_p1_i, fcell_j   , fcell_k   , 1.0);
   increment_cell_data(spatial_cell, fcell_i   , fcell_p1_j, fcell_k   , 1.0);
   increment_cell_data(spatial_cell, fcell_i   , fcell_j   , fcell_p1_k, 1.0);
   increment_cell_data(spatial_cell, fcell_i   , fcell_p1_j, fcell_p1_k, 1.0);
   increment_cell_data(spatial_cell, fcell_p1_i, fcell_j   , fcell_p1_k, 1.0);
   increment_cell_data(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_k   , 1.0);
   increment_cell_data(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_p1_k, 1.0);

}




/* Computes which blocks the distribution function will flow. These
 * will be used to compute the upstream Lagrangian cells. The upstream
 * blocks are also created here if they did not exist. data values for
 * cells are assumed to be zero(!), and will still be zero after this function
  \par spatial_cell Spatial cell that is integrated
  \par transform The euclidian acceleration transform (forward int time)
  \par downstream_blocks List of blocks which will have content 
*/

void compute_downstream_blocks(SpatialCell *spatial_cell,const Transform<Real,3,Affine>& transform,std::vector<unsigned int>& downstream_blocks) {
   
   // Make a copy of the blocklist as we don't want to iterate over blocks added by this function
   std::vector<unsigned int> blocks;
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      blocks.push_back(spatial_cell->velocity_block_list[block_i]);
   }

   
   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
     Velocity_Block* block_ptr = spatial_cell->at(blocks[block_i]);
      const Real dvx=block_ptr->parameters[BlockParams::DVX];
      const Real dvy=block_ptr->parameters[BlockParams::DVY];
      const Real dvz=block_ptr->parameters[BlockParams::DVZ];
      /* shifted to start in middle of cells*/
      const Real block_start_vx=block_ptr->parameters[BlockParams::VXCRD] + 0.5*dvx;
      const Real block_start_vy=block_ptr->parameters[BlockParams::VYCRD] + 0.5*dvy;
      const Real block_start_vz=block_ptr->parameters[BlockParams::VZCRD] + 0.5*dvz;
      
      //loop over cells 
      for (unsigned int cell_xi = 0; cell_xi < WID; cell_xi++) {
         for (unsigned int cell_yi = 0; cell_yi < WID; cell_yi++) {
            for (unsigned int cell_zi = 0; cell_zi < WID; cell_zi++) {
               /*this is the center of the cell*/
               const Eigen::Matrix<Real,3,1> s_node_position(block_start_vx + cell_xi*dvx,
                                                             block_start_vy + cell_yi*dvy,
                                                             block_start_vz + cell_zi*dvz);
               /* and where it was rotated*/
               const Eigen::Matrix<Real,3,1> s_node_position_tf=transform*s_node_position;
               mark_downstream_cells(spatial_cell,s_node_position_tf.matrix());
                        
            }
         }
      }
   }
   
   //compute downstream block list based on marked cell datas. 
   downstream_blocks.clear();
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
     unsigned int block = spatial_cell->velocity_block_list[block_i];
     Velocity_Block* block_ptr = spatial_cell->at(block);
     bool is_downstream_block=false;
     for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
       if(block_ptr->data[cell] >0.0 )
	 is_downstream_block=true;
       block_ptr->data[cell] = 0.0;
     }
     if(is_downstream_block)
       downstream_blocks.push_back(block);
   }

}

#endif
