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
   return 0.0;
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
   
   // Bx
   if (y >= -0.25 && y <= 0.2)
     if (z >= -0.25 && z <= 0.2)
       cellParams[CellParams::BX   ] = 1.0;
    
   /*
   // By
   if (x >= -0.25 && x <= 0.2)
     if (z >= -0.25 && z <= 0.2)
       cellParams[CellParams::BY   ] = 1.0;
    */
   /*
   // Bz
   if (x >= -0.25 && x <= 0.2) 
     if (y >= -0.25 && y <= 0.2)
       cellParams[CellParams::BZ   ] = 1.0;
   */
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

