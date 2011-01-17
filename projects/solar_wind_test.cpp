/*
Test particle test in GUMICS E and B fields.

Copyright 2010, 2011 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "boost/array.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "string"

#include "cell_spatial.h"
#include "common.h"
#include "logger.h"
#include "parameters.h"
#include "project.h"
#include "solar_wind.hpp"

using namespace std;
using namespace boost::posix_time;

extern Logger logger;

solar_wind::Solar_wind sw;

// tell to main.cpp that cell parameters should always be recalculated, figure out in calcSimParameters whether to actually
bool cellParametersChanged(creal& t) {return true;}

// in SI units
#define kB 1.38e-23
#define Q_P 1.602176487e-19
#define M_P 1.672621637e-27
#define R_E 6.3712e6

// cells closer than this to origin are ignored by the test particle simulation
#define INNER_BOUNDARY (7 * R_E)

Real get_initial_velocity_cell_value(const Real vx, const Real vy, const Real vz, const Real sw_vx, const Real sw_vy, const Real sw_vz, const Real T)
{
	Real f_x, f_y, f_z;

	f_x = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vx - sw_vx) * (vx - sw_vx) / 2 / kB / T);
	f_y = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vy - sw_vy) * (vy - sw_vy) / 2 / kB / T);
	f_z = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vz - sw_vz) * (vz - sw_vz) / 2 / kB / T);

	// bring the total value closer to 1 for accuracy, final scaling is done later anyway
	return f_x * f_y * f_z * 1e20;
}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz, creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
	return 0;
}

void calcBlockParameters(Real* blockParams) {
   blockParams[BlockParams::Q_PER_M] = Q_P / M_P;
}

ptime get_EB_time(string EB_file_name)
{
   EB_file_name.erase(0, 6);
   EB_file_name.erase(15, 3);
   EB_file_name[8] = 'T';

   ptime time(from_iso_string(EB_file_name));
   return time;
}

const string first_EB_file_name("mstate20040218_160000.EB");
const ptime first_EB_time = get_EB_time(first_EB_file_name);
ptime loaded_EB_time(not_a_date_time);

// Returns 0.2 * the minimum allowed timestep based on given EB file for the current grid
#ifndef PARGRID
Real get_min_timestep(dccrg<SpatialCell>& mpiGrid, string filename)
#else
#error get_min_timestep not supported with PARGRID
#endif
{
	// maximum absolute values and maximum absolute values of components in current grid from EB files
	double max_E = 0, max_B = 0, E[3] = {0}, B[3] = {0};

	// load EB data from gumics
	int result;
	FILE* infile = fopen(filename.c_str(), "rb");
	if (infile == NULL) {
		cerr << "Couldn't open file " << filename << " in get_min_timestep" << endl;
		exit(EXIT_FAILURE);
	}

	uint64_t number_of_values;
	result = fread(&number_of_values, sizeof(uint64_t), 1, infile);
	if (result != 1) {
		cerr << "number_of_values read unsuccessfull" << endl;
		exit(EXIT_FAILURE);
	}

	for (uint64_t i = 0; i < number_of_values; i++) {

		double r[3];
		result = fread(r, sizeof(r), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell location" << endl;
			exit(EXIT_FAILURE);
		}

		double current_E[3];
		result = fread(current_E, sizeof(current_E), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell E" << endl;
			exit(EXIT_FAILURE);
		}

		double current_B[3];
		result = fread(current_B, sizeof(current_B), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell B" << endl;
			exit(EXIT_FAILURE);
		}

		// exclude inner boundary
		if (sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) < INNER_BOUNDARY) {
			continue;
		}

		uint64_t cell = mpiGrid.get_cell(r[0], r[1], r[2]);
		if (cell == 0) {
			continue;
		}

		double temp_E = sqrt(current_E[0] * current_E[0] + current_E[1] * current_E[1] + current_E[2] * current_E[2]);
		if (max_E < temp_E) {
			max_E = temp_E;
		}

		double temp_B = sqrt(current_B[0] * current_B[0] + current_B[1] * current_B[1] + current_B[2] * current_B[2]);
		if (max_B < temp_B) {
			max_B = temp_B;
		}

		if (E[0] < fabs(current_E[0])) {
			E[0] = fabs(current_E[0]);
		}
		if (E[1] < fabs(current_E[1])) {
			E[1] = fabs(current_E[1]);
		}
		if (E[2] < fabs(current_E[2])) {
			E[2] = fabs(current_E[2]);
		}

		if (B[0] < fabs(current_B[0])) {
			B[0] = fabs(current_B[0]);
		}
		if (B[1] < fabs(current_B[1])) {
			B[1] = fabs(current_B[1]);
		}
		if (B[2] < fabs(current_B[2])) {
			B[2] = fabs(current_B[2]);
		}
	}

	// get minimum dx and dvx
	typedef Parameters P;
	Real min_dx = numeric_limits<Real>::max();
	Real min_dvx = numeric_limits<Real>::max();
	vector<uint64_t> cells = mpiGrid.get_cells();
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		creal dx = mpiGrid.get_cell_x_size(*cell);
		creal dy = mpiGrid.get_cell_y_size(*cell);
		creal dz = mpiGrid.get_cell_z_size(*cell);

		min_dx = min(min_dx, min(dx, min(dy, dz)));

		// TODO: take into account smaller than initial velocity blocks
		creal block_dvx = (P::vxmax - P::vxmin) / P::vxblocks_ini;
		creal block_dvy = (P::vymax - P::vymin) / P::vyblocks_ini;
		creal block_dvz = (P::vzmax - P::vzmin) / P::vzblocks_ini;
		creal cell_dvx = block_dvx / WID;
		creal cell_dvy = block_dvy / WID;
		creal cell_dvz = block_dvz / WID;

		min_dvx = min(min_dvx, min(cell_dvx, min(cell_dvy, cell_dvz)));
	}

	// TODO: use maximum Q_PER_M from the simulation, in case there are more than one species of particles
	Real max_v = max(fabs(P::vxmin), max(fabs(P::vymin), max(fabs(P::vzmin), max(fabs(P::vxmax), max(fabs(P::vymax), fabs(P::vzmax))))));
	Real a_E = Q_P / M_P * max_E;
	Real a_B = Q_P / M_P * max_v * max_B;
	return 0.4 * min(min_dx / max_v, min(min_dvx / a_E, min_dvx / a_B));
}

/*!
Returns the minimum allowed timestep based on all EB files in the current directory for the current grid.
Also initializes the solar wind.
*/
#ifndef PARGRID
Real get_min_timestep(dccrg<SpatialCell>& mpiGrid)
#else
#error get_min_timestep not supported with PARGRID
#endif
{
	sw.load(Parameters::solar_wind_file);

	Real min_dt = numeric_limits<Real>::max();

	logger << "Reading EB files." << endl;

	// loop from 16:00 to 23:59
	for (int minute = 0; minute <= 1 /*(23 - 16) * 60 + 59*/; minute++) {
		ptime current_time = first_EB_time + minutes(minute);

		string EB_filename("mstate");
		EB_filename.append(to_iso_string(current_time));
		EB_filename[14] = '_';
		EB_filename.append(".EB");

		Real dt = get_min_timestep(mpiGrid, EB_filename);
		if (min_dt > dt) {
			min_dt = dt;
		}
	}

	logger << "Reading of EB files finished." << endl;

	return min_dt;
}

#ifndef PARGRID
void load_EB(dccrg<SpatialCell>& mpiGrid, string filename)
#else
#error get_min_timestep not supported with PARGRID
#endif
{
	// zero out fields before loading
	vector<uint64_t> cells = mpiGrid.get_cells();
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		SpatialCell* cell_data = mpiGrid[*cell];
		cell_data->cpu_cellParams[CellParams::EX] = 0;
		cell_data->cpu_cellParams[CellParams::EY] = 0;
		cell_data->cpu_cellParams[CellParams::EZ] = 0;
		cell_data->cpu_cellParams[CellParams::BX] = 0;
		cell_data->cpu_cellParams[CellParams::BY] = 0;
		cell_data->cpu_cellParams[CellParams::BZ] = 0;
	}

	// load EB data from gumics
	int result;
	FILE* infile = fopen(filename.c_str(), "rb");
	if (infile == NULL) {
		cerr << "Couldn't open file " << filename << " in load_EB" << endl;
		exit(EXIT_FAILURE);
	}

	// store fields in memory for fast searching, in case given grid is finer than the one used in gumics
	vector<boost::array<double, 9> > EB;

	uint64_t number_of_values;
	result = fread(&number_of_values, sizeof(uint64_t), 1, infile);
	if (result != 1) {
		cerr << "number_of_values read unsuccessfull" << endl;
		exit(EXIT_FAILURE);
	}

	for (uint64_t i = 0; i < number_of_values; i++) {

		double r[3];
		result = fread(r, sizeof(r), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell location" << endl;
			exit(EXIT_FAILURE);
		}

		double E[3];
		result = fread(E, sizeof(E), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell E" << endl;
			exit(EXIT_FAILURE);
		}

		double B[3];
		result = fread(B, sizeof(B), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell B" << endl;
			exit(EXIT_FAILURE);
		}

		// exclude inner boundary
		if (sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) < INNER_BOUNDARY) {
			continue;
		}

		uint64_t cell = mpiGrid.get_cell(r[0], r[1], r[2]);
		if (cell == 0) {
			continue;
		}

		boost::array<double, 9> data = {{r[0], r[1], r[2], E[0], E[1], E[2], B[0], B[1], B[2]}};
		EB.push_back(data);
	}

	// put gumics fields into given grid
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		const double x = mpiGrid.get_cell_x(*cell);
		const double y = mpiGrid.get_cell_y(*cell);
		const double z = mpiGrid.get_cell_z(*cell);
		const double dx = mpiGrid.get_cell_x_size(*cell);
		const double dy = mpiGrid.get_cell_y_size(*cell);
		const double dz = mpiGrid.get_cell_z_size(*cell);

		SpatialCell* cell_data = mpiGrid[*cell];

		// exclude inner boundary
		if (sqrt(x * x + y * y + z * z) < INNER_BOUNDARY) {
			cell_data->cpu_cellParams[CellParams::EX] = 0;
			cell_data->cpu_cellParams[CellParams::EY] = 0;
			cell_data->cpu_cellParams[CellParams::EZ] = 0;
			cell_data->cpu_cellParams[CellParams::BX] = 0;
			cell_data->cpu_cellParams[CellParams::BY] = 0;
			cell_data->cpu_cellParams[CellParams::BZ] = 0;
			continue;
		}

		// if gumics didn't have the fields in this cell, use the closest fields that gumics had
		int number_of_values = 0;
		double closest = numeric_limits<double>::max();
		boost::array<double, 9> closest_EB;

		for (vector<boost::array<double, 9> >::const_iterator item = EB.begin(); item != EB.end(); item++) {
			const double distance_x = fabs(x - (*item)[0]);
			const double distance_y = fabs(y - (*item)[1]);
			const double distance_z = fabs(z - (*item)[2]);
			const double distance = sqrt(distance_x * distance_x + distance_y * distance_y + distance_z * distance_z);
			if (closest > distance) {
				closest = distance;
				closest_EB = (*item);
			}

			if (distance_x > dx / 2) {
				continue;
			}
			if (distance_y > dy / 2) {
				continue;
			}
			if (distance_z > dz / 2) {
				continue;
			}
			number_of_values++;

			cell_data->cpu_cellParams[CellParams::EX] += (*item)[3];
			cell_data->cpu_cellParams[CellParams::EY] += (*item)[4];
			cell_data->cpu_cellParams[CellParams::EZ] += (*item)[5];
			cell_data->cpu_cellParams[CellParams::BX] += (*item)[6];
			cell_data->cpu_cellParams[CellParams::BY] += (*item)[7];
			cell_data->cpu_cellParams[CellParams::BZ] += (*item)[8];
		}

		if (number_of_values == 0) {
			cell_data->cpu_cellParams[CellParams::EX] += closest_EB[3];
			cell_data->cpu_cellParams[CellParams::EY] += closest_EB[4];
			cell_data->cpu_cellParams[CellParams::EZ] += closest_EB[5];
			cell_data->cpu_cellParams[CellParams::BX] += closest_EB[6];
			cell_data->cpu_cellParams[CellParams::BY] += closest_EB[7];
			cell_data->cpu_cellParams[CellParams::BZ] += closest_EB[8];
		} else {
			cell_data->cpu_cellParams[CellParams::EX] /= number_of_values;
			cell_data->cpu_cellParams[CellParams::EY] /= number_of_values;
			cell_data->cpu_cellParams[CellParams::EZ] /= number_of_values;
			cell_data->cpu_cellParams[CellParams::BX] /= number_of_values;
			cell_data->cpu_cellParams[CellParams::BY] /= number_of_values;
			cell_data->cpu_cellParams[CellParams::BZ] /= number_of_values;
		}
	}

	fclose(infile);
}

#ifndef PARGRID
void apply_boundaries(dccrg<SpatialCell>& mpiGrid, ptime time)
#else
#error get_min_timestep not supported with PARGRID
#endif
{
	solar_wind::solar_wind_t sw_data = sw.get_solar_wind(time);

	typedef Parameters P;

	vector<uint64_t> cells = mpiGrid.get_cells();
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		double x = mpiGrid.get_cell_x(*cell);
		double y = mpiGrid.get_cell_y(*cell);
		double z = mpiGrid.get_cell_z(*cell);

		double dx = mpiGrid.get_cell_x_size(*cell);
		double dy = mpiGrid.get_cell_y_size(*cell);
		double dz = mpiGrid.get_cell_z_size(*cell);

		// solar wind
		// TODO: doesn't support variable sized velocity blocks
		if (fabs(x - 16 * R_E) < dx) {

			const Real* const blockParams = mpiGrid[*cell]->cpu_blockParams;
			const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];

			// preliminary distribution
			Real total_amount = 0;

			creal block_dvx = (P::vxmax - P::vxmin) / P::vxblocks_ini;	// Size of a block in vx-direction
			creal block_dvy = (P::vymax - P::vymin) / P::vyblocks_ini;
			creal block_dvz = (P::vzmax - P::vzmin) / P::vzblocks_ini;
			creal cell_dvx = block_dvx / WID;	// Size of one cell in a block in vx-direction
			creal cell_dvy = block_dvy / WID;
			creal cell_dvz = block_dvz / WID;

			for (uint block_zi = 0; block_zi < P::vzblocks_ini; block_zi++)
			for (uint block_yi = 0; block_yi < P::vyblocks_ini; block_yi++)
			for (uint block_xi = 0; block_xi < P::vxblocks_ini; block_xi++) {

				cuint velIndex = block_xi + block_yi * P::vxblocks_ini + block_zi * P::vyblocks_ini * P::vxblocks_ini;
				creal block_vx = P::vxmin + block_xi * block_dvx; // vx-coordinate of the lower left corner
				creal block_vy = P::vymin + block_yi * block_dvy;
				creal block_vz = P::vzmin + block_zi * block_dvz;

				for (uint cell_zi = 0; cell_zi < WID; cell_zi++)
				for (uint cell_yi = 0; cell_yi < WID; cell_yi++)
				for (uint cell_xi = 0; cell_xi < WID; cell_xi++) {
					creal cell_vx = block_vx + cell_xi * cell_dvx;
					creal cell_vy = block_vy + cell_yi * cell_dvy;
					creal cell_vz = block_vz + cell_zi * cell_dvz;

					mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi] = get_initial_velocity_cell_value(cell_vx, cell_vy, cell_vz, sw_data[solar_wind::vx], sw_data[solar_wind::vy], sw_data[solar_wind::vz], sw_data[solar_wind::T]);
					total_amount += mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi];
				}
			}

			Real current_density = total_amount / (dx * dy * dz);

			// final distribution
			for (uint block_zi = 0; block_zi < P::vzblocks_ini; block_zi++)
			for (uint block_yi = 0; block_yi < P::vyblocks_ini; block_yi++)
			for (uint block_xi = 0; block_xi < P::vxblocks_ini; block_xi++) {

				cuint velIndex = block_xi + block_yi * P::vxblocks_ini + block_zi * P::vyblocks_ini * P::vxblocks_ini;
				creal block_vx = P::vxmin + block_xi * block_dvx; // vx-coordinate of the lower left corner
				creal block_vy = P::vymin + block_yi * block_dvy;
				creal block_vz = P::vzmin + block_zi * block_dvz;

				for (uint cell_zi = 0; cell_zi < WID; cell_zi++)
				for (uint cell_yi = 0; cell_yi < WID; cell_yi++)
				for (uint cell_xi = 0; cell_xi < WID; cell_xi++) {
					creal cell_vx = block_vx + cell_xi * cell_dvx;
					creal cell_vy = block_vy + cell_yi * cell_dvy;
					creal cell_vz = block_vz + cell_zi * cell_dvz;

					mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi] *= sw_data[solar_wind::density] / current_density / DV3;
				}
			}
		}

		// "ionosphere"
		if (sqrt(x * x + y * y + z * z) < INNER_BOUNDARY) {
			for (uint block_zi = 0; block_zi < P::vzblocks_ini; block_zi++) {
			for (uint block_yi = 0; block_yi < P::vyblocks_ini; block_yi++) {
			for (uint block_xi = 0; block_xi < P::vxblocks_ini; block_xi++) {

				cuint velIndex = block_xi + block_yi * P::vxblocks_ini + block_zi * P::vyblocks_ini * P::vxblocks_ini;
				for (uint cell_zi = 0; cell_zi < WID; cell_zi++) {
				for (uint cell_yi = 0; cell_yi < WID; cell_yi++) {
				for (uint cell_xi = 0; cell_xi < WID; cell_xi++) {
					mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi] = 0;
				}}}
			}}}
		}
	}
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
// Possibly changes given maximum timestep based on loaded E & B.
#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t, Real& dt)
#else
#error get_min_timestep not supported with PARGRID
#endif
{
	ptime current_time = first_EB_time + seconds(t);

	apply_boundaries(mpiGrid, current_time);

	// check whether EB have to updated
	if (!loaded_EB_time.is_not_a_date_time()
	&& (current_time - loaded_EB_time).total_seconds() < 60) {
		return;
	}

	bool min_timestep_loaded = true;
	if (loaded_EB_time.is_not_a_date_time()) {
		min_timestep_loaded = false;
	}

	// figure out the file name to use, assume files are saved every minute
	loaded_EB_time = first_EB_time + hours((current_time - first_EB_time).hours()) + minutes((current_time - first_EB_time).minutes());
	string EB_filename("mstate");
	EB_filename.append(to_iso_string(loaded_EB_time));
	EB_filename[14] = '_';
	EB_filename.append(".EB");
	//logger << "Time given: " << t << ", loading EB from file " << EB_file_name << endl;

	load_EB(mpiGrid, EB_filename);

	if (!min_timestep_loaded) {
		dt = get_min_timestep(mpiGrid);
		logger << "Timestep: " << dt << " s" << endl;
	}
}

// initialize EB to sane values
void calcCellParameters(Real* cellParams, creal& /*t*/) {
	cellParams[CellParams::EX] = 0;
	cellParams[CellParams::EY] = 0;
	cellParams[CellParams::EZ] = 0;
	cellParams[CellParams::BX] = 0;
	cellParams[CellParams::BY] = 0;
	cellParams[CellParams::BZ] = 0;
}

