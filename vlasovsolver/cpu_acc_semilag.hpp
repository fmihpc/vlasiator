/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CPU_ACC_SEMILAG_H
#define CPU_ACC_SEMILAG_H

#include "algorithm"
#include "cmath"
#include "utility"

#include "common.h"
#include "spatial_cell.hpp"

#include <Eigen/Geometry>
#include "vlasovsolver/cpu_interpolated_block.hpp"
using namespace std;
using namespace spatial_cell;
using namespace Eigen;



/*
template<typename T> inline T cell_full_id(const T& cell_x,const T& cell_y,const T& cell_z,
                                         const T& block_x,const T& block_y,const T& block_z) {
   return WID*block_x+cell_x +
      WID*SpatialCell::vx_length*(WID*block_y+cell_y) +
      WID2*SpatialCell::vx_length*SpatialCell::vy_length*(WID*block_z+cell_z);
}


template<typename T> inline T block_id(const T& block_x,const T& block_y,const T& block_z){
   return block_x + SpatialCell::vx_length*block_y + SpatialCell::vx_length*SpatialCell::vy_length*block_z;
}
*/



/*Compute transform during on timestep, and update the bulk velocity of the cell*/

Transform<double,3,Affine> compute_acceleration_transformation( SpatialCell* spatial_cell, const double dt) {
   /*total field*/
   const double Bx = spatial_cell->parameters[CellParams::BGBXVOL]+spatial_cell->parameters[CellParams::PERBXVOL];
   const double By = spatial_cell->parameters[CellParams::BGBYVOL]+spatial_cell->parameters[CellParams::PERBYVOL];
   const double Bz = spatial_cell->parameters[CellParams::BGBZVOL]+spatial_cell->parameters[CellParams::PERBZVOL];
   /*perturbed field*/
   const double perBx = spatial_cell->parameters[CellParams::PERBXVOL];
   const double perBy = spatial_cell->parameters[CellParams::PERBYVOL];
   const double perBz = spatial_cell->parameters[CellParams::PERBZVOL];   
   //read in derivatives need for curl of B (only pertrubed, curl of background field is always 0!)
   const double dBXdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy]/spatial_cell->parameters[CellParams::DY];
   const double dBXdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz]/spatial_cell->parameters[CellParams::DZ];
   const double dBYdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx]/spatial_cell->parameters[CellParams::DX];

   const double dBYdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz]/spatial_cell->parameters[CellParams::DZ];
   const double dBZdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx]/spatial_cell->parameters[CellParams::DX];
   const double dBZdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy]/spatial_cell->parameters[CellParams::DY];

   
   const Vector3d B(Bx,By,Bz);
   const Vector3d unit_B(B.normalized());
   const double gyro_period = 2 * M_PI * Parameters::m  / (fabs(Parameters::q) * B.norm());
   //Set maximum timestep limit for this cell, based on a  maximum allowed rotation angle
   //TODO, max angle could be read in from cfg
   spatial_cell->parameters[CellParams::MAXVDT]=gyro_period*(10.0/360.0);
   
  //compute initial moments, based on actual distribution function
   spatial_cell->parameters[CellParams::RHO_V  ] = 0.0;
   spatial_cell->parameters[CellParams::RHOVX_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOVY_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOVZ_V] = 0.0;
   
   for(unsigned int block_i=0; block_i< spatial_cell->number_of_blocks;block_i++){
      unsigned int block = spatial_cell->velocity_block_list[block_i];         
      cpu_calcVelocityMoments(spatial_cell,block,CellParams::RHO_V,CellParams::RHOVX_V,CellParams::RHOVY_V,CellParams::RHOVZ_V);   
   }

   
   
   const double rho=spatial_cell->parameters[CellParams::RHO_V];
   //scale rho for hall term, if user requests
   const double hallRho =  (rho <= Parameters::lorentzHallMinimumRho ) ? Parameters::lorentzHallMinimumRho : rho ;
   const double hallPrefactor = 1.0 / (physicalconstants::MU_0 * hallRho * Parameters::q );

   Vector3d bulk_velocity(spatial_cell->parameters[CellParams::RHOVX_V]/rho,
                                 spatial_cell->parameters[CellParams::RHOVY_V]/rho,
                                 spatial_cell->parameters[CellParams::RHOVZ_V]/rho);   

   /*compute total transformation*/
   Transform<double,3,Affine> total_transform(Matrix4d::Identity());
      
   unsigned int bulk_velocity_substeps; /*!<in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep*/

   if(Parameters::lorentzHallTerm)
      bulk_velocity_substeps=dt/(gyro_period*(0.1/360.0)); 
   else
      bulk_velocity_substeps=1;
   
   /*note, we assume q is positive (pretty good assumption though)*/
   const double substeps_radians=-(2.0*M_PI*dt/gyro_period)/bulk_velocity_substeps; /*!< how many radians each substep is*/
   for(uint i=0;i<bulk_velocity_substeps;i++){
   
      /*rotation origin is the point through which we place our rotation axis (direction of which is unitB)*/
      /*first add bulk velocity (using the total transform computed this far*/
      Vector3d rotation_pivot(total_transform*bulk_velocity);
      
      if(Parameters::lorentzHallTerm) {
         //inlude lorentzHallTerm (we should include, always)      
         rotation_pivot[0]-=hallPrefactor*(dBZdy - dBYdz);
         rotation_pivot[1]-=hallPrefactor*(dBXdz - dBZdx);
         rotation_pivot[2]-=hallPrefactor*(dBYdx - dBXdy);
      }
      
      /*add to transform matrix the small rotation around  pivot
        when added like thism, and not using *= operator, the transformations
        are in the correct order
       */
      total_transform=Translation<double,3>(-rotation_pivot)*total_transform;
      total_transform=AngleAxis<double>(substeps_radians,unit_B)*total_transform;
      total_transform=Translation<double,3>(rotation_pivot)*total_transform;
      //TODO: In which order are these operations done on a point!!!
   }

   /*update bulk velocity, have not yet rotated the dist function*/
   bulk_velocity=total_transform*bulk_velocity;
   spatial_cell->parameters[CellParams::RHOVX_V] = rho*bulk_velocity[0];
   spatial_cell->parameters[CellParams::RHOVY_V] = rho*bulk_velocity[1];
   spatial_cell->parameters[CellParams::RHOVZ_V] = rho*bulk_velocity[2];

   return total_transform;
}

/*!
Propagates the distribution function in velocity space of given real space cell.

TODO:
  now this is all double: enable Real

*/

void cic_increment_cell_value(SpatialCell* spatial_cell,const int p_i,const int p_j,const int p_k,
                         const unsigned int n_subcells, const double value){
   if(p_i<0 ||p_j<0||p_k<0)
      return;

   const unsigned int block_i=p_i/(WID*n_subcells);
   const unsigned int block_j=p_j/(WID*n_subcells);
   const unsigned int block_k=p_k/(WID*n_subcells);
   const unsigned int cell_i=(p_i/n_subcells)%WID;
   const unsigned int cell_j=(p_j/n_subcells)%WID;
   const unsigned int cell_k=(p_k/n_subcells)%WID;  

   /*
     const unsigned int fcell_i=p_i/n_subcells;
     const unsigned int fcell_j=p_j/n_subcells;
     const unsigned int fcell_k=p_k/n_subcells;
     const unsigned int block_i=fcell_i/WID;
     const unsigned int block_j=fcell_j/WID;
     const unsigned int block_k=fcell_k/WID;
     const unsigned int cell_i=fcell_i-block_i*WID;
     const unsigned int cell_j=fcell_j-block_j*WID;
     const unsigned int cell_k=fcell_k-block_k*WID;
   */
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
void cic_interpolation(SpatialCell* spatial_cell,const Array3d v,const unsigned int n_subcells,const double value) {
   static int count=0;
   const double particle_dvx=SpatialCell::cell_dvx/n_subcells;
   const double particle_dvy=SpatialCell::cell_dvy/n_subcells;
   const double particle_dvz=SpatialCell::cell_dvz/n_subcells;
   const int p_i=(v[0] - SpatialCell::vx_min-0.5*particle_dvx) / particle_dvx;
   const int p_j=(v[1] - SpatialCell::vy_min-0.5*particle_dvy) / particle_dvy;
   const int p_k=(v[2] - SpatialCell::vz_min-0.5*particle_dvz) / particle_dvz;
   const double wx=(v[0]-p_i*particle_dvx - SpatialCell::vx_min-0.5*particle_dvx)/particle_dvx;
   const double wy=(v[1]-p_j*particle_dvy - SpatialCell::vy_min-0.5*particle_dvy)/particle_dvy;
   const double wz=(v[2]-p_k*particle_dvz - SpatialCell::vz_min-0.5*particle_dvz)/particle_dvz;
   if(count<5) {
      cout << p_i <<" " << p_j <<" " << p_k <<" "<< wx  <<" "<< wy <<" "<< wz <<endl;
      count++;
   }
   
   cic_increment_cell_value(spatial_cell, p_i  , p_j  , p_k  , n_subcells, (1-wx)*(1-wy)*(1-wz)*value);
   cic_increment_cell_value(spatial_cell, p_i+1, p_j  , p_k  , n_subcells,     wx*(1-wy)*(1-wz)*value);
   cic_increment_cell_value(spatial_cell, p_i  , p_j+1, p_k  , n_subcells, (1-wx)*   wy *(1-wz)*value);
   cic_increment_cell_value(spatial_cell, p_i  , p_j  , p_k+1, n_subcells, (1-wx)*(1-wy)*   wz *value);
   cic_increment_cell_value(spatial_cell, p_i  , p_j+1, p_k+1, n_subcells, (1-wx)*   wy *   wz *value);
   cic_increment_cell_value(spatial_cell, p_i+1, p_j  , p_k+1, n_subcells,    wx *(1-wy)*  -wz *value);
   cic_increment_cell_value(spatial_cell, p_i+1, p_j+1, p_k  , n_subcells,    wx *   wy *(1-wz)*value);
   cic_increment_cell_value(spatial_cell, p_i+1, p_j+1, p_k+1, n_subcells,    wx *   wy *   wz *value);   
}

/*nearest grid point*/
void ngp_interpolation(SpatialCell* spatial_cell,Array3d v,const unsigned int n_subcells,const double value) {
   spatial_cell->increment_value(v[0],v[1],v[2],value);
}


void cpu_accelerate_cell(
   SpatialCell* spatial_cell,
   const double dt) {

   phiprof::start("compute-transform");
   //compute the transform performed in this acceleration
   Transform<double,3,Affine> total_transform= compute_acceleration_transformation(spatial_cell,dt);
   phiprof::stop("compute-transform");

   
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


   /*do not change, current simple linear approximation is based on this.*/
   const Array3d  grid_min(SpatialCell::vx_min,SpatialCell::vy_min,SpatialCell::vz_min);
   const Array3d  block_dv(SpatialCell::block_dvx,SpatialCell::block_dvy,SpatialCell::block_dvz);
   const Array3d  cell_dv(SpatialCell::cell_dvx,SpatialCell::cell_dvy,SpatialCell::cell_dvz);
   
   const unsigned int n_subcells=3;
   interpolated_block iblock;

   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      const unsigned int block = blocks[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      iblock.set_block(block_ptr);

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
               
               const Vector3d s_node_position_tf=total_transform*s_node_position;
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

