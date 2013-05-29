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
using namespace std;
using namespace spatial_cell;
using namespace Eigen;





bool get_intersecting_cube(Vector3d & intersect_center,
                           Vector3d & intersect_dv,
                           const Vector3d & source_center,
                           const Vector3d & source_dv,
                           const Vector3d & target_center,
                           const Vector3d & target_dv
                           )  {
   for(int i=0;i<3;i++){
      double lower=max(source_center[i]-0.5*source_dv[i],target_center[i]-0.5*target_dv[i]);
      double upper=min(source_center[i]+0.5*source_dv[i],target_center[i]+0.5*target_dv[i]);
      if(lower>=upper)
         return false;
      intersect_center[i]=0.5*(upper+lower);
      intersect_dv[i]=upper-lower;
   }
   return true;
}

                         
                         
   

/*!
Propagates the distribution function in velocity space of given real space cell.

TODO:
  now this is all double: enable Real

*/
void cpu_accelerate_cell(
   SpatialCell* spatial_cell,
   const double dt) {


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
   //Compute maximum timestep limit for this cell, based ona  maximum allowed rotation angle
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

   
   /*transform that scales normal velocity units (m/s) to block indices*/
   // Transform<double,3,Affine> velocity_to_indices(Matrix4d::Identity()); 
   //   velocity_to_indices*=Translation<double,3>(Vector3d(
   
   
   /*compute total transformation*/
   Transform<double,3,Affine> total_transform(Matrix4d::Identity());
      
   const unsigned int bulk_velocity_substeps=dt/(gyro_period*(0.1/360.0)); /*!<in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep*/
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
   

   const int subcells=1;
   const Vector3d cell_dv(SpatialCell::cell_dvx,
                             SpatialCell::cell_dvy,
                             SpatialCell::cell_dvz);
   const Vector3d subcell_dv(cell_dv/subcells);
   

   /*do the actual accerelation operation*/
   /*PERF TODO, 
     Instead of doing these transformations in velocity units, we could do them in index units (assuming dvx=dvy=dvz =>
     would get rid of a lot of computing back and forth of velocities and indices.
     Use Array3d to get rid of index based computation
   */
   /*QUALITY TODO, interpolations, better integration (take into account overlap, hwo to combine with interpolation?*/
   
   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      const unsigned int block = blocks[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         const double cell_vx_min = SpatialCell::get_velocity_cell_vx_min(block, cell);
         const double cell_vy_min = SpatialCell::get_velocity_cell_vy_min(block, cell);
         const double cell_vz_min = SpatialCell::get_velocity_cell_vz_min(block, cell);
         
         const double subcell_rho = block_ptr->fx[cell] /(subcells*subcells*subcells);
	 

         for (int x_i = 0; x_i < subcells; x_i++)
            for (int y_i = 0; y_i < subcells; y_i++)
               for (int z_i = 0; z_i < subcells; z_i++) {
		 
                  Vector3d subcell_center(cell_vx_min + (x_i + 0.5) * subcell_dv[0],
                                          cell_vy_min + (y_i + 0.5) * subcell_dv[1],
                                          cell_vz_min + (z_i + 0.5) * subcell_dv[2]);
                  /*rotate subcell to new position*/
                  subcell_center=total_transform*subcell_center;
		 
                  /*map it back,*/

                  /*go through all 8 corners to find all potential overlapping cubes*/
                  /*hew we disregard to rotation of the subcell cube, ok for small angles*/
                  set< pair<unsigned int,unsigned int>> completed_targets;
                  
                  for(int cornerx=-1;cornerx<=1;cornerx+=2) 
                     for(int cornery=-1;cornery<=1;cornery+=2) 
                        for(int cornerz=-1;cornerz<=1;cornerz+=2) {
                           unsigned int new_block=SpatialCell::get_velocity_block(subcell_center[0]+cornerx*0.5*subcell_dv[0],
                                                                                  subcell_center[1]+cornery*0.5*subcell_dv[1],
                                                                                  subcell_center[2]+cornerz*0.5*subcell_dv[2]);
                           unsigned int new_cell=SpatialCell::get_velocity_cell(new_block,
                                                                                subcell_center[0]+cornerx*0.5*subcell_dv[0],
                                                                                subcell_center[1]+cornery*0.5*subcell_dv[1],
                                                                                subcell_center[2]+cornerz*0.5*subcell_dv[2]);

                           pair<unsigned int,unsigned int> target_cell(new_block,new_cell);
                           /*Only add contributions to cells we have not
                             computed yet!*/
                           if(completed_targets.find(target_cell) == completed_targets.end() ){
                              completed_targets.insert(target_cell);
                              const double new_cell_vx_min = SpatialCell::get_velocity_cell_vx_min(new_block, new_cell);
                              const double new_cell_vy_min = SpatialCell::get_velocity_cell_vy_min(new_block, new_cell);
                              const double new_cell_vz_min = SpatialCell::get_velocity_cell_vz_min(new_block, new_cell);
                              Vector3d new_cell_center(new_cell_vx_min + 0.5 * cell_dv[0],
                                                       new_cell_vy_min + 0.5 * cell_dv[1],
                                                       new_cell_vz_min + 0.5 * cell_dv[2]);
//                              cout<< "target cell is block "<<new_block << " cell "<<new_cell << " at " << new_cell_center[0] <<" "<< new_cell_center[1] <<" "<< new_cell_center[2] <<endl;
                           
                              Vector3d intersecting_center,intersecting_dv;
                              if(get_intersecting_cube(intersecting_center,intersecting_dv,
                                                       subcell_center,subcell_dv,
                                                       new_cell_center,cell_dv)){
                                 /*add if they are intersecting*/
                                 // cout <<"adding to " << new_block << " " << new_cell <<" rho of " << subcell_rho*(intersecting_dv[0]*intersecting_dv[1]*intersecting_dv[2])/( subcell_dv[0]*subcell_dv[1]*subcell_dv[2])<<endl;
                                 spatial_cell->increment_value(new_block,new_cell,
                                                               subcell_rho*
                                                               (intersecting_dv[0]*intersecting_dv[1]*intersecting_dv[2])/
                                                               ( subcell_dv[0]*subcell_dv[1]*subcell_dv[2]));
                                 
                              }
                           }
                        }
	       }
      }
   }
}
   
   

#endif

