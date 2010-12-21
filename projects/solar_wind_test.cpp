#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/unordered_map.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>

#include "cell_spatial.h"
#include "common.h"
#include "logger.h"
#include "parameters.h"
#include "project.h"

using namespace std;
using namespace boost::posix_time;

extern Logger logger;

// tell to main.cpp that cell parameters should always be recalculated, figure out in calcSimParameters whether to actually
bool cellParametersChanged(creal& t) {return true;}

// in SI units
#define T 1e5
#define SOLAR_WIND_DENSITY 7.3e6
#define SOLAR_WIND_V_X (-4e5)
#define SOLAR_WIND_V_Y 0
#define SOLAR_WIND_V_Z 0
#define kB 1.38e-23

#define Q_P 1.602176487e-19
#define M_P 1.672621637e-27

#define R_E 6.3712e6

Real get_initial_velocity_cell_value(const Real vx, const Real vy, const Real vz)
{
	Real f_x, f_y, f_z;

	f_x = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vx - SOLAR_WIND_V_X) * (vx - SOLAR_WIND_V_X) / 2 / kB / T);
	f_y = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vy - SOLAR_WIND_V_Y) * (vy - SOLAR_WIND_V_Y) / 2 / kB / T);
	f_z = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vz - SOLAR_WIND_V_Z) * (vz - SOLAR_WIND_V_Z) / 2 / kB / T);

	// bring the total value closer to 1 for accuracy, final scaling is done later anyway
	return f_x * f_y * f_z * 1e50;
}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz, creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {

	if (fabs(x - 28 * R_E) > dx) {
		return 0;
	}

	return get_initial_velocity_cell_value(vx, vy, vz);
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

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t) {

	// apply the solar wind boundary
	vector<uint64_t> cells = mpiGrid.get_cells();
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		double x = mpiGrid.get_cell_x(*cell);
		double dx = mpiGrid.get_cell_x_size(*cell);

		if (fabs(x - 28 * R_E) > dx) {
			continue;
		}

		// preliminary distribution

		Real total_amount = 0;

		double y = mpiGrid.get_cell_y(*cell);
		double z = mpiGrid.get_cell_z(*cell);
		double dy = mpiGrid.get_cell_y_size(*cell);
		double dz = mpiGrid.get_cell_z_size(*cell);

		typedef Parameters P;

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

				mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi] = get_initial_velocity_cell_value(cell_vx, cell_vy, cell_vz);
				total_amount += mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi];
			}
		}

		Real density = total_amount / (dx * dy * dz);

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

				mpiGrid[*cell]->cpu_avgs[velIndex * SIZE_VELBLOCK + cell_zi * WID2 + cell_yi * WID + cell_xi] *= SOLAR_WIND_DENSITY / density;
			}
		}
	}

	ptime current_time = first_EB_time + seconds(t);

	// check whether EB have to updated
	if (!loaded_EB_time.is_not_a_date_time()
	&& (current_time - loaded_EB_time).total_seconds() < 60) {
		return;
	}

	// zero out fields before loading
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		mpiGrid[*cell]->cpu_cellParams[CellParams::EX] = 0;
		mpiGrid[*cell]->cpu_cellParams[CellParams::EY] = 0;
		mpiGrid[*cell]->cpu_cellParams[CellParams::EZ] = 0;
		mpiGrid[*cell]->cpu_cellParams[CellParams::BX] = 0;
		mpiGrid[*cell]->cpu_cellParams[CellParams::BY] = 0;
		mpiGrid[*cell]->cpu_cellParams[CellParams::BZ] = 0;
	}

	// figure out the file name to use
	loaded_EB_time = first_EB_time + hours((current_time - first_EB_time).hours()) + minutes((current_time - first_EB_time).minutes());
	string EB_file_name("mstate");
	EB_file_name.append(to_iso_string(loaded_EB_time));
	EB_file_name[14] = '_';
	EB_file_name.append(".EB");
	//logger << "Time given: " << t << ", loading EB from file " << EB_file_name << endl;

	// load EB data from gumics
	boost::unordered_map<uint64_t, unsigned int> values_per_cell;
	int result;

	FILE* infile = fopen(EB_file_name.c_str(), "rb");
	if (infile == NULL) {
		cerr << "Couldn't open file " << EB_file_name << endl;
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
			cerr << "r read unsuccessfull" << endl;
			exit(EXIT_FAILURE);
		}

		uint64_t cell = mpiGrid.get_cell(r[0], r[1], r[2]);
		if (cell == 0) {
			continue;
		}

		if (values_per_cell.count(cell) == 0) {
			values_per_cell[cell] = 1;
		} else {
			values_per_cell[cell]++;
		}

		result = fseek(infile, sizeof(double) * 6, SEEK_CUR);
		if (result != 0) {
			cerr << "SEEK_CUR unsuccessfull" << endl;
			exit(EXIT_FAILURE);
		}
	}

	// Read values into cells
	result = fseek(infile, sizeof(uint64_t), SEEK_SET);
	if (result != 0) {
		cerr << "SEEK_SET unsuccessfull" << endl;
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
		if (E[0] > 1) {
			cerr << "Too large Ex" << endl;
			exit(EXIT_FAILURE);
		}
		if (E[1] > 1) {
			cerr << "Too large Ey" << endl;
			exit(EXIT_FAILURE);
		}
		if (E[2] > 1) {
			cerr << "Too large Ez" << endl;
			exit(EXIT_FAILURE);
		}

		double B[3];
		result = fread(B, sizeof(B), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell B" << endl;
			exit(EXIT_FAILURE);
		}
		if (B[0] > 1) {
			cerr << "Too large Bx" << endl;
			exit(EXIT_FAILURE);
		}
		if (B[1] > 1) {
			cerr << "Too large By" << endl;
			exit(EXIT_FAILURE);
		}
		if (B[2] > 1) {
			cerr << "Too large Bz" << endl;
			exit(EXIT_FAILURE);
		}

		uint64_t cell = mpiGrid.get_cell(r[0], r[1], r[2]);
		if (cell == 0) {
			continue;
		}

		SpatialCell* cell_data = mpiGrid[cell];
		if (cell_data == NULL) {
			continue;
		}

		if (values_per_cell.count(cell) == 0 || values_per_cell[cell] == 0) {
			// gumics didn't have the fields
			mpiGrid[cell]->cpu_cellParams[CellParams::EX] = 0;
			mpiGrid[cell]->cpu_cellParams[CellParams::EY] = 0;
			mpiGrid[cell]->cpu_cellParams[CellParams::EZ] = 0;
			mpiGrid[cell]->cpu_cellParams[CellParams::BX] = 0;
			mpiGrid[cell]->cpu_cellParams[CellParams::BY] = 0;
			mpiGrid[cell]->cpu_cellParams[CellParams::BZ] = 0;
		} else {
			mpiGrid[cell]->cpu_cellParams[CellParams::EX] += E[0] / values_per_cell[cell];
			mpiGrid[cell]->cpu_cellParams[CellParams::EY] += E[1] / values_per_cell[cell];
			mpiGrid[cell]->cpu_cellParams[CellParams::EZ] += E[2] / values_per_cell[cell];
			mpiGrid[cell]->cpu_cellParams[CellParams::BX] += B[0] / values_per_cell[cell];
			mpiGrid[cell]->cpu_cellParams[CellParams::BY] += B[1] / values_per_cell[cell];
			mpiGrid[cell]->cpu_cellParams[CellParams::BZ] += B[2] / values_per_cell[cell];
		}
	}

	// print maximum field in cells
	logger << "Maximum field in current grid from file " << EB_file_name << ": ";
	Real B_max[3] = {numeric_limits<Real>::min(), numeric_limits<Real>::min(), numeric_limits<Real>::min()};
	Real E_max[3] = {numeric_limits<Real>::min(), numeric_limits<Real>::min(), numeric_limits<Real>::min()};
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		if (E_max[0] < mpiGrid[*cell]->cpu_cellParams[CellParams::EX]) {
			E_max[0] = mpiGrid[*cell]->cpu_cellParams[CellParams::EX];
		}
		if (E_max[1] < mpiGrid[*cell]->cpu_cellParams[CellParams::EY]) {
			E_max[1] = mpiGrid[*cell]->cpu_cellParams[CellParams::EY];
		}
		if (E_max[2] < mpiGrid[*cell]->cpu_cellParams[CellParams::EZ]) {
			E_max[2] = mpiGrid[*cell]->cpu_cellParams[CellParams::EZ];
		}
		if (B_max[0] < mpiGrid[*cell]->cpu_cellParams[CellParams::BX]) {
			B_max[0] = mpiGrid[*cell]->cpu_cellParams[CellParams::BX];
		}
		if (B_max[1] < mpiGrid[*cell]->cpu_cellParams[CellParams::BY]) {
			B_max[1] = mpiGrid[*cell]->cpu_cellParams[CellParams::BY];
		}
		if (B_max[2] < mpiGrid[*cell]->cpu_cellParams[CellParams::BZ]) {
			B_max[2] = mpiGrid[*cell]->cpu_cellParams[CellParams::BZ];
		}
	}
	logger << "Ex " << E_max[0] << ", Ey " << E_max[1] << ", Ez " << E_max[2] << ", Bx " << B_max[0] << ", By " << B_max[1] << ", Bz " << B_max[2] << endl;

	fclose(infile);
}
#else
#error Gumics EB tracing not supported with PARGRID
#endif

// initialize EB to sane values
void calcCellParameters(Real* cellParams, creal& /*t*/) {
	cellParams[CellParams::EX] = 0;
	cellParams[CellParams::EY] = 0;
	cellParams[CellParams::EZ] = 0;
	cellParams[CellParams::BX] = 0;
	cellParams[CellParams::BY] = 0;
	cellParams[CellParams::BZ] = 0;
}

