#include <cstdlib>
#include <iostream>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"

using namespace std;

typedef testTrans2Parameters TTP;
Real TTP::radLimitInf = NAN;
Real TTP::radLimitSup = NAN;
Real TTP::vRad = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("test_trans2.radLimitInf", "radial inner radius", 0.2);
   RP::add("test_trans2.radLimitSup", "radial outer radius", 0.5);
   RP::add("test_trans2.vRad", "radial velocity", 0.5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("test_trans2.radLimitInf", TTP::radLimitInf);
   RP::get("test_trans2.radLimitSup", TTP::radLimitSup);
   RP::get("test_trans2.vRad", TTP::vRad);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal X = x + 0.5*dx;
   creal Y = y + 0.5*dy;
   creal vX = vx + 0.5*dvx;
   creal vY = vy + 0.5*dvy;
   
   creal radius = sqrt(X*X + Y*Y);
   creal cell2dID = (int) (X / dx) +
   (int) (Y / dy) * Parameters::xcells_ini;
   creal vel2dID = (int) (vX / dvx) +
   (int) (vY / dvy) * Parameters::vxblocks_ini * 4;
   
   if (cell2dID == vel2dID &&
       vz + 0.5 * dvz >= -0.4 * dvz &&
       vz + 0.5 * dvz < 0.4 * dvz &&
       radius <= TTP::radLimitSup &&
       radius >= TTP::radLimitInf
       ) {
      return 1.0;
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
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

