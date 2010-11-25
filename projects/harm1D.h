#include "definitions.h"

/** Calculate parameters for the given velocity grid block.
 * Here you need to set a value for array indices: BlockParams::Q_PER_M.
 * @param blockParams Array containing the block parameters
 */
void calcBlockParameters(Real* blockParams);

/** Calculate parameters for the given spatial cell.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
 * 
 * The following array indices contain the coordinates of the "lower left corner" of the cell: 
 * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
 * 
 * The correct cell size is given in the following array indices: 
 * CellParams::DX, CellParams::DY, and CellParams::DZ.
 * 
 * @param cellParams Array containing cell parameters.
 */
void calcCellParameters(Real* cellParams);

/** Calculate the phase space density at the given phase space coordinates.
 * @param x X-position.
 * @param y Y-position.
 * @param z Z-position.
 * @param vx VX-position.
 * @param vy VY-position.
 * @param vz VZ-position.
 * @return The value of the distribution function at the given phase space coordinates.
 */
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& vx,creal& vy,creal& vz);

template<typename T> T spatialFluxX(cuint& i,const T& avg_neg,const T& avg_pos,const T* const blockParams) {
   //creal VX = 0.25;
   creal VX = blockParams[BlockParams::VXCRD] + i*blockParams[BlockParams::DVX];
   return convert<T>(0.5)*VX*(avg_neg+avg_pos) - convert<T>(0.5)*fabs(VX)*(avg_pos-avg_neg);
}

template<typename T> T spatialFluxY(cuint& j,const T& avg_neg,const T& avg_pos,const T* const blockParams) {
   //creal VY = 0.0;
   creal VY = blockParams[BlockParams::VYCRD] + j*blockParams[BlockParams::DVY];
   return convert<T>(0.5)*VY*(avg_neg+avg_pos) - convert<T>(0.5)*fabs(VY)*(avg_pos-avg_neg);
}

template<typename T> T spatialFluxZ(cuint& k,const T& avg_neg,const T& avg_pos,const T* const blockParams) {
   //creal VZ = 0.0;
   creal VZ = blockParams[BlockParams::VZCRD] + k*blockParams[BlockParams::DVZ];
   return convert<T>(0.5)*VZ*(avg_neg+avg_pos) - convert<T>(0.5)*fabs(VZ)*(avg_pos-avg_neg);
}

template<typename T> T velocityFluxX(const T& j,const T& k,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VY = blockParams[BlockParams::VYCRD] + (j+convert<T>(0.5))*blockParams[BlockParams::DVY];
   const T VZ = blockParams[BlockParams::VZCRD] + (k+convert<T>(0.5))*blockParams[BlockParams::DVZ];
   const T EX = cellParams[CellParams::EX];
   const T BY = cellParams[CellParams::BY];
   const T BZ = cellParams[CellParams::BZ];
   const T AX = blockParams[BlockParams::Q_PER_M]*(EX + VY*BZ - VZ*BY);
   return convert<T>(0.5)*AX*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AX)*(avg_pos-avg_neg);
}

template<typename T> T velocityFluxX(cuint& i,cuint& j,cuint& k,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VY = blockParams[BlockParams::VYCRD] + convert<T>((j+0.5))*blockParams[BlockParams::DVY];
   const T VZ = blockParams[BlockParams::VZCRD] + convert<T>((k+0.5))*blockParams[BlockParams::DVZ];
   //const T VY = blockParams[BlockParams::VYCRD] + (k+0.5)*blockParams[BlockParams::DVY];
   //const T VZ = blockParams[BlockParams::VZCRD] + (j+0.5)*blockParams[BlockParams::DVZ];
   const T EX = cellParams[CellParams::EX];
   const T BY = cellParams[CellParams::BY];
   const T BZ = cellParams[CellParams::BZ];
   
   //const T AX = 0.25;
   const T AX = blockParams[BlockParams::Q_PER_M]*(EX + VY*BZ - VZ*BY);
   return convert<T>(0.5)*AX*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AX)*(avg_pos-avg_neg);
}

template<typename T> T velocityFluxY(const T& I,const T& K,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   return convert<T>(0.0);
}

template<typename T> T velocityFluxY(cuint& i,cuint& j,cuint& k,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   /*
   const T VX = blockParams[BlockParams::VXCRD] + (i+0.5)*blockParams[BlockParams::DVX];
   //const T VZ = blockParams[BlockParams::VZCRD] + (j+0.5)*blockParams[BlockParams::DVZ];
   const T VZ = blockParams[BlockParams::VZCRD] + (k+0.5)*blockParams[BlockParams::DVZ];
   const T EY = cellParams[CellParams::EY];
   const T BX = cellParams[CellParams::BX];
   const T BZ = cellParams[CellParams::BZ];
   
   const T AY = blockParams[BlockParams::Q_PER_M]*(EY + VZ*BX - VX*BZ);
   return 0.5*AY*(avg_neg + avg_pos) - 0.5*fabs(AY)*(avg_pos-avg_neg);*/
   return 0.0;
}

template<typename T> T velocityFluxZ(const T& I,const T& J,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   return convert<T>(0.0);
}

template<typename T> T velocityFluxZ(cuint& i,cuint& j,cuint& k,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   /*
   const T VX = blockParams[BlockParams::VXCRD] + (i+0.5)*blockParams[BlockParams::DVX];
   const T VY = blockParams[BlockParams::VYCRD] + (j+0.5)*blockParams[BlockParams::DVY];
   const T EZ = cellParams[CellParams::EZ];
   const T BX = cellParams[CellParams::BX];
   const T BY = cellParams[CellParams::BY];
   
   const T AZ = blockParams[BlockParams::Q_PER_M]*(EZ + VX*BY - VY*BX);
   return 0.5*AZ*(avg_neg + avg_pos) - 0.5*fabs(AZ)*(avg_pos-avg_neg);*/
   return 0.0;
}


