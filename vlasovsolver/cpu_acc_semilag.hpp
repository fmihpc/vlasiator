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

using namespace std;
using namespace spatial_cell;


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



   
   // don't iterate over blocks added by this function
   std::vector<unsigned int> blocks;
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      blocks.push_back(spatial_cell->velocity_block_list[block_i]);
   }

   
//   orbit_time = 2 * M_PI * Parameters::m  / (fabs(Parameters::q) * B_abs);
        
   /*loop over velocity cells */
   for(unsigned int block_i = 0; block_i < blocks.size(); block_i++) {
      Velocity_Block* block=spatial_cell->at(blocks[block_i]);
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         const unsigned int block=blocks[block_i];
         const double cell_vx_min = SpatialCell::get_velocity_cell_vx_min(block, cell);
         const double cell_vy_min = SpatialCell::get_velocity_cell_vy_min(block, cell);
         const double cell_vz_min = SpatialCell::get_velocity_cell_vz_min(block, cell);
         const double distribution_function = spatial_cell->at(block)->data[cell];

         /*
         if (!spatial_cell->add_velocity_block(target_block)) {
            std::cerr << __FILE__ << ":" << __LINE__
                      << " Couldn't add target velocity block " << target_block
                      << std::endl;
            abort();
         }
         */

      }
   }


/*!
Applies fluxes of the distribution function in given spatial cell.

Overwrites current cell data with fluxes where we stored rotated data.
*/

   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      const unsigned int block = spatial_cell->velocity_block_list[block_i];
      
      Velocity_Block* block_ptr = spatial_cell->at(block);

      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         block_ptr->data[cell] = block_ptr->fx[cell];
      }
   }
}


#endif

