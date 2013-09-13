#ifndef CPU_CIC_H
#define CPU_CIC_H

#include "spatial_cell.hpp"
#include <Eigen/Geometry>
#include "../parameters.h"

using namespace Eigen;

inline void cic_increment_cell_value(SpatialCell* spatial_cell,
                                     const unsigned int fcell_i, const unsigned int fcell_j,const unsigned int fcell_k,
                                     const unsigned int n_subcells, const Real value){
   

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
      //outside outer bounds -> density is lost
      const Real DV3=SpatialCell::cell_dvx*SpatialCell::cell_dvy*SpatialCell::cell_dvz;
      spatial_cell->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*value;
      return;
   }         

   spatial_cell->increment_value(block,cell,value);
}



/*volme integrated version of cic*/
/*
  \par spatial_cell Spatial cell that is integrated
  \par v velocity lower left corner of subcell after rotation
  \par iblock interpolated block describing the block to be rotated
  \par v_source Original position of v prior to rotation
  \par ib_cellid  cellid in interpolated block corresponding to v_source
  \par n_subcells Number of subcells per dimension
*/

inline void volintegrated_interpolation(SpatialCell* spatial_cell,const Array3d v,interpolated_block& iblock,const Array3d v_source,unsigned int ib_cellid,const unsigned int n_subcells){
   static int count=0;
   const Real particle_dvx=SpatialCell::cell_dvx/n_subcells;
   const Real particle_dvy=SpatialCell::cell_dvy/n_subcells;
   const Real particle_dvz=SpatialCell::cell_dvz/n_subcells;

   //note that p_ijk can also have a value of -1 right at the edge(!)
   //When we compute p_i we compare inside which shifted-by-halfvxyz
   //cube the velocity is. This gives us information about the eight
   //potential targets to which we will add density. p_ijk is the
   //index of the lower left one
   const int p_i=(v[0] - SpatialCell::vx_min-0.5*particle_dvx) / particle_dvx;
   const int p_j=(v[1] - SpatialCell::vy_min-0.5*particle_dvy) / particle_dvy;
   const int p_k=(v[2] - SpatialCell::vz_min-0.5*particle_dvz) / particle_dvz;
   const unsigned int fcell_i=p_i/n_subcells;
   const unsigned int fcell_j=p_j/n_subcells;
   const unsigned int fcell_k=p_k/n_subcells;
   const unsigned int fcell_p1_i=(p_i+1)/n_subcells;
   const unsigned int fcell_p1_j=(p_j+1)/n_subcells;
   const unsigned int fcell_p1_k=(p_k+1)/n_subcells;
   const Real subcell_vol_frac=1.0/(n_subcells*n_subcells*n_subcells); //relative part of the total volume for the one subcell represented by v

   //midpoint of the lower left cell
   const Real midpoint_x=p_i*particle_dvx + SpatialCell::vx_min+0.5*particle_dvx;
   const Real midpoint_y=p_j*particle_dvy + SpatialCell::vy_min+0.5*particle_dvy;
   const Real midpoint_z=p_k*particle_dvz + SpatialCell::vz_min+0.5*particle_dvz;

   //weights for CIC
   const Real wx=(v[0]-midpoint_x)/particle_dvx;
   const Real wy=(v[1]-midpoint_y)/particle_dvy;
   const Real wz=(v[2]-midpoint_z)/particle_dvz;
//includes no rotation
   const Real to_source_frame_x=v_source[0]-v[0];
   const Real to_source_frame_y=v_source[1]-v[1];
   const Real to_source_frame_z=v_source[2]-v[2];

   
   const Real ic_vx    = 0.5*(midpoint_x+v[0])+to_source_frame_x;
   const Real ic_vy    = 0.5*(midpoint_y+v[1])+to_source_frame_y;
   const Real ic_vz    = 0.5*(midpoint_z+v[2])+to_source_frame_z;
   const Real ic_p1_vx = ic_vx+0.5*particle_dvx;
   const Real ic_p1_vy = ic_vy+0.5*particle_dvy;
   const Real ic_p1_vz = ic_vz+0.5*particle_dvz;


   if (p_i < 0 || p_j < 0 || p_k < 0 ) {
      //we do not try to handle the edge of (total) velocity space in any smart way, this cell is essentiall lost
      const Real DV3=SpatialCell::cell_dvx*SpatialCell::cell_dvy*SpatialCell::cell_dvz;
      spatial_cell->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*iblock.get_value(ib_cellid,v_source[0],v_source[1],v_source[2]);
      return;
   }
   
   
   if(fcell_i==fcell_p1_i && fcell_j==fcell_p1_j && fcell_k==fcell_p1_k){
      cic_increment_cell_value(spatial_cell, fcell_i  , fcell_j  , fcell_k  , n_subcells,
                               iblock.get_value(ib_cellid,v_source[0],v_source[1],v_source[2])*subcell_vol_frac);
   }
   //   else if (fcell_i==fcell_p1_i && fcell_j==fcell_p1_j)  {
   //   else if (fcell_j==fcell_p1_j && fcell_k==fcell_p1_k)  {
   //   else if (fcell_i==fcell_p1_i && fcell_k==fcell_p1_k)  {
   //   else if (fcell_i==fcell_p1_i)  {
   //   else if (fcell_j==fcell_p1_j)  {
   //   else if (fcell_k==fcell_p1_k)  {
   else{
     //these are the midpoints of the intersecting cubes formed by the
     //particle cubic clound with the target grid. It is in source
     //frame, as that is what we use to read in from iblock the value
     //of the middle of the cube. For a hinged hyperplane
     //interpolation this will * volume give us and exact integration

   
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_k   , n_subcells,(1-wx)*(1-wy)*(1-wz) * subcell_vol_frac * iblock.get_value(ib_cellid,ic_vx   ,ic_vy   ,ic_vz   ));
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_k   , n_subcells,    wx*(1-wy)*(1-wz) * subcell_vol_frac * iblock.get_value(ib_cellid,ic_p1_vx,ic_vy   ,ic_vz   ));
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_k   , n_subcells,(1-wx)*   wy *(1-wz) * subcell_vol_frac * iblock.get_value(ib_cellid,ic_vx   ,ic_p1_vy,ic_vz   ));
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_j   , fcell_p1_k, n_subcells,(1-wx)*(1-wy)*   wz  * subcell_vol_frac * iblock.get_value(ib_cellid,ic_vx   ,ic_vy   ,ic_p1_vz));
      cic_increment_cell_value(spatial_cell, fcell_i   , fcell_p1_j, fcell_p1_k, n_subcells,(1-wx)*   wy *   wz  * subcell_vol_frac * iblock.get_value(ib_cellid,ic_vx   ,ic_p1_vy,ic_p1_vz));
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_j   , fcell_p1_k, n_subcells,   wx *(1-wy)*   wz  * subcell_vol_frac * iblock.get_value(ib_cellid,ic_p1_vx,ic_vy   ,ic_p1_vz));
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_k   , n_subcells,   wx *   wy *(1-wz) * subcell_vol_frac * iblock.get_value(ib_cellid,ic_p1_vx,ic_p1_vy,ic_vz   ));
      cic_increment_cell_value(spatial_cell, fcell_p1_i, fcell_p1_j, fcell_p1_k, n_subcells,   wx *   wy *   wz  * subcell_vol_frac * iblock.get_value(ib_cellid,ic_p1_vx,ic_p1_vy,ic_p1_vz));
   }
}





void cic(SpatialCell *spatial_cell,Transform<Real,3,Affine>& transform) {
   
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

   //n_subcells should be known at compile time, otherwise mapping is ~ 25% slower (!)
   const unsigned int n_subcells=1;//Parameters::semiLagSubcellsPerDim;
   BlockInterpolationType interpolation;
   switch(Parameters::vlasovSemiLagOrder){
       case 0:
          interpolation=CONSTANT;
          break;
       case 1:
          interpolation=LINEAR;
       case 2:
          interpolation=PARABOLIC;
       default:
          interpolation=LINEAR;
   }
   
   interpolated_block iblock(interpolation);
   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      const unsigned int block = blocks[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      iblock.load_block(block_ptr);

      const Real dvx=block_ptr->parameters[BlockParams::DVX]/n_subcells;
      const Real dvy=block_ptr->parameters[BlockParams::DVY]/n_subcells;
      const Real dvz=block_ptr->parameters[BlockParams::DVZ]/n_subcells;
      const Real block_start_vx=block_ptr->parameters[BlockParams::VXCRD] + 0.5*dvx;
      const Real block_start_vy=block_ptr->parameters[BlockParams::VYCRD] + 0.5*dvy;
      const Real block_start_vz=block_ptr->parameters[BlockParams::VZCRD] + 0.5*dvz;
      
      //loop over internal points in block
      //todo, split into siz loops, outer over cells, inner over subcells. Better memory usage, and we need cellidfor get value (fast version)
      for (unsigned int cell_xi = 0; cell_xi < WID; cell_xi++) {
         for (unsigned int cell_yi = 0; cell_yi < WID; cell_yi++) {
            for (unsigned int cell_zi = 0; cell_zi < WID; cell_zi++) {
               unsigned int ib_cellid=iblock.get_cell_id(cell_xi,cell_yi,cell_zi);
               for (unsigned int subcell_xi = 0; subcell_xi < n_subcells; subcell_xi++) {
                  for (unsigned int subcell_yi = 0; subcell_yi < n_subcells; subcell_yi++) {
                     for (unsigned int subcell_zi = 0; subcell_zi < n_subcells; subcell_zi++) {
                        
                        const Eigen::Matrix<Real,3,1> s_node_position(block_start_vx + (cell_xi*n_subcells+subcell_xi)*dvx,
                                                                        block_start_vy + (cell_yi*n_subcells+subcell_yi)*dvy,
                                                                        block_start_vz + (cell_zi*n_subcells+subcell_zi)*dvz);

                        const Eigen::Matrix<Real,3,1> s_node_position_tf=transform*s_node_position;
                        volintegrated_interpolation(spatial_cell,s_node_position_tf.matrix(),iblock,s_node_position.matrix(),ib_cellid,n_subcells);
                        //scaling, just to test things...
                        /*
                          Real value=iblock.get_value(s_node_position[0],s_node_position[1],s_node_position[2]);
                          const Matrix<Real,3> s_node_position_tf(n_subcells*s_node_position[0],
                          n_subcells*s_node_position[1],
                          n_subcells*s_node_position[2]);
                          ngp_interpolation(spatial_cell,s_node_position_tf.matrix(),n_subcells,value);
                        */
                     }
                  }
               }
            }
         }
      }
   }
}

#endif
