#include <cstdlib>
#include <iostream>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"
#include "parameters.h"

enum cases {BXCASE,BYCASE,BZCASE};

// SET THIS TO BXCASE,BYCASE,OR BZCASE TO SELECT ONE OF THE THREE CASES

static int CASE = BYCASE;

using namespace std;

#ifndef PARGRID
bool initializeProject(dccrg<SpatialCell>& mpiGrid) {
#else
bool initializeProject(ParGrid<SpatialCell>& mpiGrid) {
#endif
   return true;
}


bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   const Real T  = 1e-6;
   const Real n  = 1.0e7;

   Real VX,VY,VZ;
   switch (CASE) {
    case BXCASE:
      VX = 0.0;
      VY = 0.5;
      VZ = 0.5;
      break;
    case BYCASE:
      VX = 0.5;
      VY = 0.0;
      VZ = 0.5;
      break;
    case BZCASE:
      VX = 0.5;
      VY = 0.5;
      VZ = 0.0;
      break;
   }
   
   creal VX2 = (vx+0.5*dvx-VX)*(vx+0.5*dvx-VX);
   creal VY2 = (vy+0.5*dvy-VY)*(vy+0.5*dvy-VY);
   creal VZ2 = (vz+0.5*dvz-VZ)*(vz+0.5*dvz-VZ);
   
   creal CONST = physicalconstants::MASS_PROTON / 2.0 / physicalconstants::K_B / T;
   Real NORM = (physicalconstants::MASS_PROTON / 2.0 / M_PI / physicalconstants::K_B / T);
   NORM = n * pow(NORM,1.5);
   
   return NORM*exp(-CONST*(VX2+VY2+VZ2));
}

void calcBlockParameters(Real* blockParams) { }

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = 0.0;
   
   typedef Parameters P;
   creal x = cellParams[CellParams::XCRD];
   creal y = cellParams[CellParams::YCRD];
   creal z = cellParams[CellParams::ZCRD];

   creal B0 = 1.0e-9;
   
   switch (CASE) {
    case BXCASE:
      if (y >= -0.25 && y <= 0.2)
	if (z >= -0.25 && z <= 0.2)
	  cellParams[CellParams::BX] = B0;
      break;
    case BYCASE:
      if (x >= -0.25 && x <= 0.2)
	if (z >= -0.25 && z <= 0.2)
	  cellParams[CellParams::BY] = B0;
      break;
    case BZCASE:
      if (x >= -0.25 && x <= 0.2)
	if (y >= -0.25 && y <= 0.2)
	  cellParams[CellParams::BZ] = B0;
      break;
   }
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

