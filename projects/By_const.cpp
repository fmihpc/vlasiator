#include <cstdlib>
#include <iostream>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"
#include "parameters.h"

using namespace std;

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   // Every spatial cell gets the same velocity distribution.
   // The velocity distribution is a shell distribution centered at v = V0.
   
   const Real V0 = 0.5; // Center point of Gaussian distrib.
   const Real DV = 0.1; // Sigma of Gaussian distrib.   
   
   const Real DV2 = DV*DV;                 // Square of sigma
   const Real NORM = convert<Real>(1.0) / pow(convert<Real>(2.0)*M_PI*DV2,convert<Real>(3.0/2.0)); // Normalization const.
   const Real PHI_MIN = -20.0*M_PI/180.0;
   const Real PHI_MAX = +20.0*M_PI/180.0;
   
   // Sample the distribution with N points in each coordinate direction:
   const uint N = 2;
   Real result = 0.0;
   for (uint k=0; k<N; ++k) {
      const Real VZ = vz + (k+0.5)*dvz/N;
      for (uint j=0; j<N; ++j) {
	 const Real VY = vy + (j+0.5)*dvy/N;
	 for (uint i=0; i<N; ++i) {
	    const Real VX = vx + (i+0.5)*dvx/N;
	    const Real V = sqrt(VX*VX+VZ*VZ);
	    Real PHI = acos(VZ/V);
	    if (VX < 0.0) PHI *= -1.0;
	    
	    if (PHI < PHI_MIN || PHI > PHI_MAX) return 0.0;
	    result += exp(-convert<Real>(0.5)*(V-V0)*(V-V0)/DV2)*exp(-convert<Real>(0.5)*VY*VY/DV2);
	 }
      }
   }
   return result * NORM / (N*N*N);
}

void calcBlockParameters(Real* blockParams) {
   blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 1.0;
   cellParams[CellParams::BZ   ] = 0.0;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
#else
void calcSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t) {
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
#endif
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->cpu_cellParams, t);
   }
}

