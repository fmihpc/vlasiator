#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "definitions.h"

cuint MAX_SPA_CELLS = 240;
cuint MAX_VEL_BLOCKS = 30000;

namespace Transmit {
   cuint CELL_PARAMS  = 1;
   cuint BLOCK_PARAMS = CELL_PARAMS << 1;
   cuint AVGS         = CELL_PARAMS << 2;
   cuint FLUXES       = CELL_PARAMS << 3;
   cuint DERIV1       = CELL_PARAMS << 4;
   cuint DERIV2       = CELL_PARAMS << 5;
   cuint NBRSVEL      = CELL_PARAMS << 6;
};

struct Parameters {
   static real xmin;  /**< X-coordinate of the lower left corner of the spatial grid. */
   static real xmax;  /**< Y-coordinate of the lower left corner of the spatial grid. */
   static real ymin;  /**< Z-coordinate of the lower left corner of the spatial grid. */
   static real ymax;  /**< X-coordinate of the upper right corner of the spatial grid. */
   static real zmin;  /**< Y-coordinate of the upper right corner of the spatial grid. */
   static real zmax;  /**< Z-coordinate of the upper right corner of the spatial grid. */
   static real dx_ini; /**< Initial size of spatial cell in x-direction. */
   static real dy_ini; /**< Initial size of spatial cell in y-direction. */
   static real dz_ini; /**< Initial size of spatial cell in z-direction. */
   
   static real vxmin; /**< VX-coordinate of the lower left corner of velocity grid. */
   static real vxmax; /**< VY-coordinate of the lower left corner of velocity grid. */
   static real vymin; /**< VZ-coordinate of the lower left corner of velocity grid. */
   static real vymax; /**< VX-coordinate of the upper right corner of velocity grid. */
   static real vzmin; /**< VY-coordinate of the upper right corner of velocity grid. */
   static real vzmax; /**< VZ-coordinate of the upper right corner of velocity grid. */
   
   static uint xcells_ini; /**< Initial number of spatial cells in x-direction. */
   static uint ycells_ini; /**< Initial number of spatial cells in y-direction. */
   static uint zcells_ini; /**< Initial number of spatial cells in z-direction. */
   static uint vxblocks_ini; /**< Initial number of velocity grid blocks in vx-direction. */
   static uint vyblocks_ini; /**< Initial number of velocity grid blocks in vy-direction. */
   static uint vzblocks_ini; /**< Initial number of velocity grid blocks in vz-direction. */

   static real dt;     /**< The value of the timestep to use in propagation. */
   static uint tstep;  /**< The number of the current timestep. 0=initial state. */
   static uint tsteps; /**< Total number of timesteps to calculate. */
   static uint saveInterval;
   static uint diagnInterval;
   
   static uint transmit; /**< Indicates the data that needs to be transmitted to remote nodes.
			  * This is created with bitwise or from the values defined in 
			  * namespace Transmit.*/
   
   Parameters();
};

#endif
