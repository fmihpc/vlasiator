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

/*!
Return the relative shared volume of cubes specified by given arguments.

Returned value is relative to length1^3 and is between 0 and 1.
*/
double get_relative_shared_volume(
	const boost::array<double, 3>& center1,
	const double length1,
	const boost::array<double, 3>& center2,
	const double length2
) {
	const double overlapping_x = std::min(length1, std::min(length2, (length1 + length2) / 2 - fabs(center1[0] - center2[0]))),
		overlapping_y = std::min(length1, std::min(length2, (length1 + length2) / 2 - fabs(center1[1] - center2[1]))),
		overlapping_z = std::min(length1, std::min(length2, (length1 + length2) / 2 - fabs(center1[2] - center2[2])));
	return std::max(0.0, overlapping_x)
		* std::max(0.0, overlapping_y)
		* std::max(0.0, overlapping_z)
		/ std::min(length1 * length1 * length1, length2 * length2 * length2);
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
   

   const int subcells=3;
   const double sub_dvx = SpatialCell::cell_dvx / subcells,
      sub_dvy = SpatialCell::cell_dvy / subcells,
      sub_dvz = SpatialCell::cell_dvz / subcells;        

   /*do the actual accerelation operation*/
   /*PERF TODO, instead of doing these transformations in velocity units, we could do them in index units (assuming dvx=dvy=dvz =>
     would get rid of a lot of computing back and forth of velocities and indices. */
   /*QUALITY TODO, interpolations, better integration (take into account overlap, hwo to combine with interpolation?*/
   
   for (unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      const unsigned int block = blocks[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         const double cell_vx_min = SpatialCell::get_velocity_cell_vx_min(block, cell);
         const double cell_vy_min = SpatialCell::get_velocity_cell_vy_min(block, cell);
         const double cell_vz_min = SpatialCell::get_velocity_cell_vz_min(block, cell);
         
         const double subcell_rho = block_ptr->fx[cell] / (subcells*subcells*subcells);

         for (int x_i = 0; x_i < subcells; x_i++)
            for (int y_i = 0; y_i < subcells; y_i++)
               for (int z_i = 0; z_i < subcells; z_i++) {
                  Vector3d subcell_position(cell_vx_min + (x_i + 0.5) * sub_dvx,
                                            cell_vy_min + (y_i + 0.5) * sub_dvy,
                                            cell_vz_min + (z_i + 0.5) * sub_dvz);
                  Vector3d subcell_new_position=total_transform*subcell_position;
                  spatial_cell->increment_value(subcell_new_position[0],
                                                subcell_new_position[1],
                                                subcell_new_position[2],subcell_rho);
               }
      }
   }
}


#endif

