
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


   
   /*copy distribution function values into the flux table, and zero the existing distribution function*/
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      const unsigned int block = spatial_cell->velocity_block_list[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         block_ptr->fx[cell] = block_ptr->data[cell];
         block_ptr->data[cell] = 0.0;
      }
   }
   
   // Make a copy of the blocklist as we don't want to iterate over blocks added by this function
   std::vector<unsigned int> blocks;
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      blocks.push_back(spatial_cell->velocity_block_list[block_i]);
   }

   /*compute total transformation*/
   Transform<double,3,Affine> total_transform;
   const unsigned int bulk_velocity_substeps=dt/(gyro_period*(0.1/360.0)); /*!<in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep*/
   const double substeps_radians=(2.0*M_PI*dt/gyro_period)/bulk_velocity_substeps; /*!< how many radians each substep is*/
   for(uint i=0;i<bulk_velocity_substeps;i++){
   
      /*rotation origin is the point through which we place our rotation axis (direction of which is unitB)*/
      /*first add bulk velocity (using the total transform computed this far*/
      Vector3d rotation_origin(total_transform*bulk_velocity);

      if(Parameters::lorentzHallTerm) {
         //inlude lorentzHallTerm (we should include, always)      
         rotation_origin[0]-=hallPrefactor*(dBZdy - dBYdz);
         rotation_origin[1]-=hallPrefactor*(dBXdz - dBZdx);
         rotation_origin[2]-=hallPrefactor*(dBYdx - dBXdy);
      }
      
      /*add to transform matrix*/
      total_transform*=Translation<double,3>(-rotation_origin);
      total_transform*=AngleAxis<double>(substeps_radians,unit_B);
      total_transform*=Translation<double,3>(rotation_origin);
   }

   cout<< total_transform.matrix();
   
}


#endif

