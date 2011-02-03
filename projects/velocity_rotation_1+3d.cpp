#include <cstdlib>
#include <iostream>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"

using namespace std;

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz, creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
	#define AMP (1.0 / dx / dy / dz / dvx / dvy / dvz)
	#define DELTAX 50000
	#define rad0 (1e6 / 2)

	double phi = atan2(vy, vx);
	double rad = sqrt(vx * vx + vy * vy);

	if (phi < 3 * M_PI / 4 && phi > M_PI / 4) {
		return AMP * exp(-(rad - rad0) * (rad - rad0) / DELTAX / DELTAX);
	} else {
		return 0;
	}
}

void calcBlockParameters(Real* blockParams) {
   #define Q_P 1.602176487e-19
   #define M_P 1.672621637e-27
   //blockParams[BlockParams::Q_PER_M] = Q_P / M_P;
}

void calcCellParameters(Real* cellParams, creal& /*t*/) {
   cellParams[CellParams::EX] = 0.0;
   cellParams[CellParams::EY] = 0.0;
   cellParams[CellParams::EZ] = 0.0;
   cellParams[CellParams::BX] = 0.0;
   cellParams[CellParams::BY] = 0.0;
   cellParams[CellParams::BZ] = 1e-9;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
#else
void calcSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
#endif
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->cpu_cellParams, t);
   }
}

