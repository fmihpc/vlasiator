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
#include "vlasovsolver/cpu_cic.hpp"
using namespace std;
using namespace spatial_cell;
using namespace Eigen;

/*

TODO

Handle v space boundaries properly and quickly
Use real not double (or own float datatype?)
remove einspline (useless)
Use agner's vectorclass (perhaps the 3 vectors with position?)
test nsubcells and make to parameter


*/

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

void cpu_accelerate_cell(SpatialCell* spatial_cell,const double dt) {

   phiprof::start("compute-transform");
   //compute the transform performed in this acceleration
   Transform<double,3,Affine> total_transform= compute_acceleration_transformation(spatial_cell,dt);
   phiprof::stop("compute-transform");
   cic(spatial_cell,total_transform);
}
   

#endif

