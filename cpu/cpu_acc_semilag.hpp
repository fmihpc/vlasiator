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
#include "boost/array.hpp"
#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "utility"

#include "common.h"
#include "spatial_cell.hpp"

using namespace std;

typedef boost::array<double, 3> coordinate_t;

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
	SpatialCell& spatial_cell,
	const double dt,
	const double steps_per_orbit,
	const double charge,
	const double mass
) {
	const double Bx = spatial_cell.parameters[CellParams::BXVOL],
		By = spatial_cell.parameters[CellParams::BYVOL],
		Bz = spatial_cell.parameters[CellParams::BZVOL],
		Ex = spatial_cell.parameters[CellParams::EXVOL],
		Ey = spatial_cell.parameters[CellParams::EYVOL],
		Ez = spatial_cell.parameters[CellParams::EZVOL];

	// don't iterate over blocks added by this function
	std::vector<unsigned int> blocks;
	for (
		unsigned int block = spatial_cell.velocity_block_list[0], block_i = 0;
		block_i < max_velocity_blocks
			&& spatial_cell.velocity_block_list[block_i] != error_velocity_block;
		block = spatial_cell.velocity_block_list[++block_i]
	) {
		blocks.push_back(block);
	}

	// calculate how many steps to take when tracing particle trajectory
	const double orbit_time = 2 * M_PI * mass / (charge * sqrt(Bx * Bx + By * By + Bz * Bz));
	unsigned int substeps = dt / orbit_time * steps_per_orbit;
	if (substeps < 1) {
		substeps = 1;
	}

	// break every velocity cell into this many subcells when accelerating
	const int x_subcells = 1,
		y_subcells = 1,
		z_subcells = 1,
		subcells = x_subcells * y_subcells * z_subcells;

	const double sub_dvx = cell_dvx / x_subcells,
		sub_dvy = cell_dvy / y_subcells,
		sub_dvz = cell_dvz / z_subcells;

	BOOST_FOREACH(unsigned int block, blocks) {

		for (unsigned int cell = 0; cell < velocity_block_len; cell++) {

			const double cell_vx_min = get_velocity_cell_vx_min(block, cell),
				cell_vy_min = get_velocity_cell_vy_min(block, cell),
				cell_vz_min = get_velocity_cell_vz_min(block, cell);

			// TODO: higher order reconstruction of distribution function in subcells
			const double distribution_function = spatial_cell.at(block)->data[cell] / subcells;

			for (int x_i = 0; x_i < x_subcells; x_i++)
			for (int y_i = 0; y_i < y_subcells; y_i++)
			for (int z_i = 0; z_i < z_subcells; z_i++) {

				// accelerate a particle with given charge and mass
				boost::array<double, 3> current_v = {{
					cell_vx_min + (x_i + 0.5) * sub_dvx,
					cell_vy_min + (y_i + 0.5) * sub_dvy,
					cell_vz_min + (z_i + 0.5) * sub_dvz
				}};

				for (unsigned int substep = 0; substep < substeps; substep++) {
					const double ax = charge * (Ex + current_v[1] * Bz - current_v[2] * By) / mass,
						ay = charge * (Ey + current_v[2] * Bx - current_v[0] * Bz) / mass,
						az = charge * (Ez + current_v[0] * By - current_v[1] * Bx) / mass;

						current_v[0] += ax * dt / substeps;
						current_v[1] += ay * dt / substeps;
						current_v[2] += az * dt / substeps;
				}

				// get all velocity block & cell pairs which overlap with current subcell
				boost::unordered_set<std::pair<unsigned int, unsigned int> > targets;

				const double dvx = sub_dvx / 2,
					dvy = sub_dvy / 2,
					dvz = sub_dvz / 2;

				const boost::array<boost::array<double, 3>, 8> target_coordinates = boost::assign::list_of
					(boost::assign::list_of(current_v[0] - dvx)(current_v[1] - dvy)(current_v[2] - dvz))
					(boost::assign::list_of(current_v[0] - dvx)(current_v[1] - dvy)(current_v[2] + dvz))
					(boost::assign::list_of(current_v[0] - dvx)(current_v[1] + dvy)(current_v[2] - dvz))
					(boost::assign::list_of(current_v[0] - dvx)(current_v[1] + dvy)(current_v[2] + dvz))
					(boost::assign::list_of(current_v[0] + dvx)(current_v[1] - dvy)(current_v[2] - dvz))
					(boost::assign::list_of(current_v[0] + dvx)(current_v[1] - dvy)(current_v[2] + dvz))
					(boost::assign::list_of(current_v[0] + dvx)(current_v[1] + dvy)(current_v[2] - dvz))
					(boost::assign::list_of(current_v[0] + dvx)(current_v[1] + dvy)(current_v[2] + dvz));

				BOOST_FOREACH(coordinate_t coordinate, target_coordinates) {
					const unsigned int target_block = get_velocity_block(
						coordinate[0],
						coordinate[1],
						coordinate[2]
					);

					if (target_block == error_velocity_block) {
						continue;
					}

					if (!spatial_cell.add_velocity_block(target_block)) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Couldn't add target velocity block " << target_block
							<< std::endl;
						abort();
					}

					const unsigned int target_cell = get_velocity_cell(
						target_block,
						coordinate[0],
						coordinate[1],
						coordinate[2]
					);

					if (target_cell == error_velocity_cell) {
						continue;
					}

					targets.insert(std::make_pair(target_block, target_cell));
				}

				// assign fluxes to target cells from current subcell
				for (boost::unordered_set<std::pair<unsigned int, unsigned int> >::const_iterator
					item = targets.begin();
					item != targets.end();
					item++
				) {
					const unsigned int target_block = item->first,
						target_cell = item->second;

					const boost::array<double, 3> target_v = {
						(get_velocity_cell_vx_min(target_block, target_cell)
							+ get_velocity_cell_vx_max(target_block, target_cell)) / 2,
						(get_velocity_cell_vy_min(target_block, target_cell)
							+ get_velocity_cell_vy_max(target_block, target_cell)) / 2,
						(get_velocity_cell_vz_min(target_block, target_cell)
							+ get_velocity_cell_vz_max(target_block, target_cell)) / 2
					};

					// assign value based on shared volume between target and source cells
					const double shared_V = get_relative_shared_volume(
						current_v,
						sub_dvx,
						target_v,
						cell_dvx
					);

					Velocity_Block* block_ptr = spatial_cell.at(target_block);
					block_ptr->fx[target_cell] += distribution_function * shared_V;
				}
			}
		}
	}
}

/*!
Applies fluxes of the distribution function in given spatial cell.

Overwrites current cell data with fluxes and zeroes fluxes.
*/
void apply_fluxes(Spatial_Cell& spatial_cell)
{
	for (
		unsigned int block = spatial_cell.velocity_block_list[0], block_i = 0;
		block_i < max_velocity_blocks
			&& spatial_cell.velocity_block_list[block_i] != error_velocity_block;
		block = spatial_cell.velocity_block_list[++block_i]
	) {
		Velocity_Block* block_ptr = spatial_cell.at(block);

		for (unsigned int cell = 0; cell < velocity_block_len; cell++) {
			if (block_ptr->fx[cell] > 1e6) {
				cout << block_ptr->fx[cell] << " in block " << block << " cell " << cell << endl;
			}
			block_ptr->data[cell] = block_ptr->fx[cell];
			block_ptr->fx[cell] = 0;
		}
	}
}

#endif

