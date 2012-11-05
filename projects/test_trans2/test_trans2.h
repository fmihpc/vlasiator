#ifndef TESTTRANS2_H
#define TESTTRANS2_H

#include "definitions.h"
#include "spatial_cell.hpp"
#include "projects/projects_common.h"
#include "projects/projects_vlasov_acceleration.h"

struct testTrans2Parameters {
   static Real radLimitInf;
   static Real radLimitSup;
} ;

/**
 * Initialize project. Can be used, e.g., to read in parameters from the input file
 */
bool initializeProject(void);

/** Register parameters that should be read in
 */
bool addProjectParameters(void);
/** Get the value that was read in
 */
bool getProjectParameters(void);

/*!\brief Set the fields and distribution of a cell according to the default simulation settings.
 * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
 * \param cell Pointer to the cell to set.
 */
void setProjectCell(SpatialCell* cell);
// WARNING Lorentz force not fixed in this project (cf. JXB term in the acceleration)!!!
template<typename UINT,typename REAL> 
void calcAccFaceX(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL HALF = 0.5;
   const REAL VX = blockParams[BlockParams::VXCRD] + I*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+HALF)*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+HALF)*blockParams[BlockParams::DVZ];
   ax = Parameters::q_per_m * (cellParams[CellParams::EX] +
   VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
   VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = Parameters::q_per_m * (cellParams[CellParams::EY] +
   VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
   VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = Parameters::q_per_m * (cellParams[CellParams::EZ] +
   VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
   VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}
// WARNING Lorentz force not fixed in this project (cf. JXB term in the acceleration)!!!
template<typename UINT,typename REAL> 
void calcAccFaceY(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL HALF = 0.5;
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+HALF)*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + J*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+HALF)*blockParams[BlockParams::DVZ];
   ax = Parameters::q_per_m * (cellParams[CellParams::EX] +
   VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
   VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = Parameters::q_per_m * (cellParams[CellParams::EY] +
   VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
   VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = Parameters::q_per_m * (cellParams[CellParams::EZ] +
   VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
   VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}
// WARNING Lorentz force not fixed in this project (cf. JXB term in the acceleration)!!!
template<typename UINT,typename REAL> 
void calcAccFaceZ(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL HALF = 0.5;
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+HALF)*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+HALF)*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + K*blockParams[BlockParams::DVZ];
   ax = Parameters::q_per_m * (cellParams[CellParams::EX] +
   VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
   VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = Parameters::q_per_m * (cellParams[CellParams::EY] +
   VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
   VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = Parameters::q_per_m * (cellParams[CellParams::EZ] +
   VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
   VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}

#endif
