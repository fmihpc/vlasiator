#ifndef CPU_CIC_H
#define CPU_CIC_H

#include "spatial_cell.hpp"
#include <Eigen/Geometry>
using namespace Eigen;

inline void cic_increment_cell_value(SpatialCell* spatial_cell,
                                     const unsigned int fcell_i, const unsigned int fcell_j,const unsigned int fcell_k,
                                     const unsigned int n_subcells, const double value){
   

   const unsigned int block_i=fcell_i/WID;
   const unsigned int block_j=fcell_j/WID;
   const unsigned int block_k=fcell_k/WID;
   const unsigned int cell_i=fcell_i-block_i*WID;
   const unsigned int cell_j=fcell_j-block_j*WID;
   const unsigned int cell_k=fcell_k-block_k*WID;
   
   if(block_i>= SpatialCell::vx_length ||
      block_j>= SpatialCell::vy_length ||
      block_k>= SpatialCell::vz_length)
      return;
         
   const unsigned int block = block_i + block_j * SpatialCell::vx_length + block_k * SpatialCell::vx_length * SpatialCell::vy_length;
   const unsigned int cell  = cell_i + cell_j * WID + cell_k * WID2;
   spatial_cell->increment_value(block,cell,value);
}

/*cloud in cell interpolation*/
//TODO, what about negative indices p_ijk, reformulate
inline void cic_interpolation(SpatialCell* spatial_cell,const Array3d v,const unsigned int n_subcells,const double value) {
   static int count=0;
   const double particle_dvx=SpatialCell::cell_dvx/n_subcells;
   const double particle_dvy=SpatialCell::cell_dvy/n_subcells;
   const double particle_dvz=SpatialCell::cell_dvz/n_subcells;
   const int p_i=(v[0] - SpatialCell::vx_min-0.5*particle_dvx) / particle_dvx;
   const int p_j=(v[1] - SpatialCell::vy_min-0.5*particle_dvy) / particle_dvy;
   const int p_k=(v[2] - SpatialCell::vz_min-0.5*particle_dvz) / particle_dvz;
   const unsigned int fcell_i=p_i/n_subcells;
   const unsigned int fcell_j=p_j/n_subcells;
   const unsigned int fcell_k=p_k/n_subcells;
   const unsigned int fcell_p1_i=(p_i+1)/n_subcells;
   const unsigned int fcell_p1_j=(p_j+1)/n_subcells;
   const unsigned int fcell_p1_k=(p_k+1)/n_subcells;


   const double wx=(fcell_i!=fcell_p1_i)?((v[0]-p_i*particle_dvx - SpatialCell::vx_min-0.5*particle_dvx)/particle_dvx):0.0;
   const double wy=(fcell_j!=fcell_p1_j)?((v[1]-p_j*particle_dvy - SpatialCell::vy_min-0.5*particle_dvy)/particle_dvy):0.0;
   const double wz=(fcell_k!=fcell_p1_k)?((v[2]-p_k*particle_dvz - SpatialCell::vz_min-0.5*particle_dvz)/particle_dvz):0.0;


//   if(p_i<0 || p_j<0 || p_k<0){
      //not goog enough, we need cheap and good test for out of bounds
   //    cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_p1_k, n_subcells,    wx *   wy *   wz *value);   
   //   return;
   // }

   if(fcell_i==fcell_p1_i && fcell_j==fcell_p1_j && fcell_k==fcell_p1_k){
      cic_increment_cell_value(spatial_cell, fcell_i  , fcell_j  , fcell_k  , n_subcells, value);
   }
   else if (fcell_i==fcell_p1_i && fcell_j==fcell_p1_j)  {
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_p1_k, n_subcells,  wz *value);
   }
   else if (fcell_j==fcell_p1_j && fcell_k==fcell_p1_k)  {
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wx)*value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_k   , n_subcells, wx*value);
   }
   else if (fcell_i==fcell_p1_i && fcell_k==fcell_p1_k)  {
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wy)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_k   , n_subcells,  wy *value);
   }
   else if (fcell_i==fcell_p1_i)  {
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wy)*(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_k   , n_subcells,    wy *(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_p1_k, n_subcells,    wy *   wz *value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_p1_k, n_subcells, (1-wy)*   wz *value);
   }
   else if (fcell_j==fcell_p1_j)  {
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wx)*(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_k   , n_subcells, wx*(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_p1_k, n_subcells, (1-wx)* wz *value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_p1_k, n_subcells,    wx * wz *value);
   }
   else if (fcell_k==fcell_p1_k)  {
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wx)*(1-wy)*value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_k   , n_subcells,     wx*(1-wy)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_k   , n_subcells, (1-wx)*   wy *value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_k   , n_subcells,    wx *   wy *value);
   }
   else{
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells, (1-wx)*(1-wy)*(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_k   , n_subcells,     wx*(1-wy)*(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_k   , n_subcells, (1-wx)*   wy *(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_p1_k, n_subcells, (1-wx)*(1-wy)*   wz *value);
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_p1_k, n_subcells, (1-wx)*   wy *   wz *value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_p1_k, n_subcells,    wx *(1-wy)*   wz *value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_k   , n_subcells,    wx *   wy *(1-wz)*value);
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_p1_k, n_subcells,    wx *   wy *   wz *value);   
   }
}





void cic(SpatialCell *spatial_cell,Transform<double,3,Affine>& transform) {
   
   // Make a copy of the blocklist as we don't want to iterate over blocks added by this function
   std::vector<unsigned int> blocks;
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      blocks.push_back(spatial_cell->velocity_block_list[block_i]);
   }

   
   /*copy distribution function values into the flux table, and zero the existing distribution function*/
   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      const unsigned int block = blocks[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         block_ptr->fx[cell] = block_ptr->data[cell];
         block_ptr->data[cell] = 0.0;
      }
   }
   
   const unsigned int n_subcells=Parameters::semiLagSubcellsPerDim;
   interpolated_block iblock(HINGED_HYPERPLANE);
   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      const unsigned int block = blocks[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      iblock.load_block(block_ptr);

      const double dvx=block_ptr->parameters[BlockParams::DVX]/n_subcells;
      const double dvy=block_ptr->parameters[BlockParams::DVY]/n_subcells;
      const double dvz=block_ptr->parameters[BlockParams::DVZ]/n_subcells;
      const double block_start_vx=block_ptr->parameters[BlockParams::VXCRD] + 0.5*dvx;
      const double block_start_vy=block_ptr->parameters[BlockParams::VYCRD] + 0.5*dvy;
      const double block_start_vz=block_ptr->parameters[BlockParams::VZCRD] + 0.5*dvz;
      
      //loop over internal points in block
      for (unsigned int cell_xi = 0; cell_xi < WID*n_subcells; cell_xi++) {
         for (unsigned int cell_yi = 0; cell_yi < WID*n_subcells; cell_yi++) {
            for (unsigned int cell_zi = 0; cell_zi < WID*n_subcells; cell_zi++) {
               const Vector3d s_node_position(block_start_vx + cell_xi*dvx,
                                              block_start_vy + cell_yi*dvy,
                                              block_start_vz + cell_zi*dvz);

               const Vector3d s_node_position_tf=transform*s_node_position;
               double value=iblock.get_value(s_node_position[0],s_node_position[1],s_node_position[2])/(n_subcells*n_subcells*n_subcells);
               cic_interpolation(spatial_cell,s_node_position_tf.matrix(),n_subcells,value);
               
               //scaling, just to test things...
               /*
               double value=iblock.get_value(s_node_position[0],s_node_position[1],s_node_position[2]);
               const Vector3d s_node_position_tf(n_subcells*s_node_position[0],
                                                 n_subcells*s_node_position[1],
                                                 n_subcells*s_node_position[2]);
               ngp_interpolation(spatial_cell,s_node_position_tf.matrix(),n_subcells,value);
               */
            }
         }
      }
   }
}

#endif
