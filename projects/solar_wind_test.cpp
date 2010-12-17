#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/unordered_map.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"

using namespace std;
using namespace boost::posix_time;

// tell to main.cpp that cell parameters should always be recalculated, figure out whether to actually in calcSimParameters
bool cellParametersChanged(creal& t) {return true;}

#define Q_P 1.602176487e-19
#define M_P 1.672621637e-27

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz, creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {

	// in SI units
	#define T 1e5
	#define SOLAR_WIND_DENSITY 7.3e6
	#define SOLAR_WIND_V_X (-4e5)
	#define SOLAR_WIND_V_Y 0
	#define SOLAR_WIND_V_Z 0
	#define kB 1.38e-23

	#define R_E 6.3712e6
	if (x < 24 * R_E) {
		return 0;
	}

	double f_x = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vx - SOLAR_WIND_V_X) * (vx - SOLAR_WIND_V_X) / 2 / kB / T);
	double f_y = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vy - SOLAR_WIND_V_Y) * (vy - SOLAR_WIND_V_Y) / 2 / kB / T);
	double f_z = sqrt(M_P / (2 * M_PI * kB * T)) * exp(-M_P * (vz - SOLAR_WIND_V_Z) * (vz - SOLAR_WIND_V_Z) / 2 / kB / T);

	return f_x * f_y * f_z * SOLAR_WIND_DENSITY;
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

	ptime current_time = first_EB_time + seconds(t);

	// check whether EB have to updated
	if (!loaded_EB_time.is_not_a_date_time()
	&& (current_time - loaded_EB_time).total_seconds() < 60) {
		return;
	}

	// figure out the file name to use
	loaded_EB_time = first_EB_time + hours((current_time - first_EB_time).hours()) + minutes((current_time - first_EB_time).minutes());
	string EB_file_name("mstate");
	EB_file_name.append(to_iso_string(loaded_EB_time));
	EB_file_name[14] = '_';
	EB_file_name.append(".EB");
	cout << "Time given: " << t << ", loading EB from file " << EB_file_name << endl;

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

		double B[3];
		result = fread(B, sizeof(B), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell B" << endl;
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

	fclose(infile);
}
#else
#error Gumics EB tracing not supported with PARGRID
#endif

void calcCellParameters(Real* /*cellParams*/, creal& /*t*/) {
}

