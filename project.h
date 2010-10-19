#include "definitions.h"

/** Calculate parameters for the given velocity grid block.
 * Here you need to set a value for array indices: BlockParams::Q_PER_M.
 * @param blockParams Array containing the block parameters
 */
void calcBlockParameters(real* blockParams);

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
void calcCellParameters(real* cellParams);

/** Calculate the phase space density at the given phase space coordinates.
 * @param x X-position.
 * @param y Y-position.
 * @param z Z-position.
 * @param vx VX-position.
 * @param vy VY-position.
 * @param vz VZ-position.
 * @return The value of the distribution function at the given phase space coordinates.
 */
real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& vx,creal& vy,creal& vz);



