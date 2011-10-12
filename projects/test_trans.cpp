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
   if (z >= 0.2 && z <= 0.2+dz) {//|| (z >= -0.38 && z <= -0.2)) {
      if (x >= -0.38 && x <= -0.2) {
	 // Lower left corner
	 if (y >= -0.38 && y <= -0.2) {
	    if (vx >= -0.5-dvx && vx < -0.5)
	      if (vy >= -0.5-dvy && vy < -0.5)
		if (vz >= 0.5 && vz < 0.5+dvz)
		  return 1.0;
	 }
	 // Upper left corner
	 if (y >= 0.2 && y <= 0.2+dy) {
	    if (vx >= -0.5-dvx && vx < -0.5)
	      if (vy >= 0.5 && vy < 0.5+dvy)
		if (vz >= 0.5 && vz < 0.5+dvz)
		  return 1.0;
	 }
      }
      
      if (x >= 0.2 && x <= 0.2+dx) {
	 // Lower right corner
	 if (y >= -0.38 && y <= -0.2) {
	    if (vx >= 0.5 && vx < 0.5+dvx)
	      if (vy >= -0.5-dvx && vy < -0.5)
		if (vz >= 0.5 && vz < 0.5+dvz)
		  return 1.0;
	 }      
	 // Upper right corner
	 if (y >= 0.2 && y <= 0.2+dy) {
	    if (vx >= 0.5 && vx < 0.5+dvx) 
	      if (vy >= 0.5 && vy < 0.5+dvy)
		if (vz >= 0.5 && vz < 0.5+dvz)
		  return 1.0;
	 }
      }
   }
   
   if (z >= -0.38 && z <= -0.2) {
      if (x >= -0.38 && x <= -0.2) {
	 // Lower left corner
	 if (y >= -0.38 && y <= -0.2) {
	    if (vx >= -0.5-dvx && vx < -0.5)
	      if (vy >= -0.5-dvy && vy < -0.5)
		if (vz >=-0.5-dvz && vz < -0.5)
		  return 1.0;
	 }
	 // Upper left corner
	 if (y >= 0.2 && y <= 0.2+dy) {
	    if (vx >= -0.5-dvx && vx < -0.5)
	      if (vy >= 0.5 && vy < 0.5+dvy)
		if (vz >=-0.5-dvz && vz < -0.5)
		  return 1.0;
	 }
      }
      
      if (x >= 0.2 && x <= 0.2+dx) {
	 // Lower right corner
	 if (y >= -0.38 && y <= -0.2) {
	    if (vx >= 0.5 && vx < 0.5+dvx)
	      if (vy >= -0.5-dvx && vy < -0.5)
		if (vz >=-0.5-dvz && vz < -0.5)
		  return 1.0;
	 }
	 // Upper right corner
	 if (y >= 0.2 && y <= 0.2+dy) {
	    if (vx >= 0.5 && vx < 0.5+dvx)
	      if (vy >= 0.5 && vy < 0.5+dvy)
		if (vz >=-0.5-dvz && vz < -0.5)
		  return 1.0;
	 }
      }
   }
   
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
   
   cellParams[CellParams::EXVOL] = 0.0;
   cellParams[CellParams::EYVOL] = 0.0;
   cellParams[CellParams::EZVOL] = 0.0;
   cellParams[CellParams::BXVOL] = 0.0;
   cellParams[CellParams::BYVOL] = 0.0;
   cellParams[CellParams::BZVOL] = 0.0;
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

