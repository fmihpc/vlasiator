#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "definitions.h"
#include <vector>
#include <string>

cuint MAX_VEL_BLOCKS = 460000;
//cuint MAX_VEL_BLOCKS = 3000000;
namespace Transmit {
   cuint CELL_PARAMS  = 1;
   cuint BLOCK_PARAMS = CELL_PARAMS << 1;
   cuint AVGS         = CELL_PARAMS << 2;
   cuint FLUXES       = CELL_PARAMS << 3;
   cuint DERIV1       = CELL_PARAMS << 4;
   cuint DERIV2       = CELL_PARAMS << 5;
   cuint NBRSVEL      = CELL_PARAMS << 6;
}

struct Parameters {
   static Real xmin;  /**< X-coordinate of the lower left corner of the spatial grid. */
   static Real xmax;  /**< Y-coordinate of the lower left corner of the spatial grid. */
   static Real ymin;  /**< Z-coordinate of the lower left corner of the spatial grid. */
   static Real ymax;  /**< X-coordinate of the upper right corner of the spatial grid. */
   static Real zmin;  /**< Y-coordinate of the upper right corner of the spatial grid. */
   static Real zmax;  /**< Z-coordinate of the upper right corner of the spatial grid. */
   static Real dx_ini; /**< Initial size of spatial cell in x-direction. */
   static Real dy_ini; /**< Initial size of spatial cell in y-direction. */
   static Real dz_ini; /**< Initial size of spatial cell in z-direction. */
   
   static Real vxmin; /**< VX-coordinate of the lower left corner of velocity grid. */
   static Real vxmax; /**< VY-coordinate of the lower left corner of velocity grid. */
   static Real vymin; /**< VZ-coordinate of the lower left corner of velocity grid. */
   static Real vymax; /**< VX-coordinate of the upper right corner of velocity grid. */
   static Real vzmin; /**< VY-coordinate of the upper right corner of velocity grid. */
   static Real vzmax; /**< VZ-coordinate of the upper right corner of velocity grid. */
   
   static uint xcells_ini; /**< Initial number of spatial cells in x-direction. */
   static uint ycells_ini; /**< Initial number of spatial cells in y-direction. */
   static uint zcells_ini; /**< Initial number of spatial cells in z-direction. */
   static uint vxblocks_ini; /**< Initial number of velocity grid blocks in vx-direction. */
   static uint vyblocks_ini; /**< Initial number of velocity grid blocks in vy-direction. */
   static uint vzblocks_ini; /**< Initial number of velocity grid blocks in vz-direction. */

   static Real t;      /**< Current simulation time. */
   static Real dt;     /**< The value of the timestep to use in propagation. */
   static uint tstep;  /**< The number of the current timestep. 0=initial state. */
   static uint tsteps; /**< Total number of timesteps to calculate. */
   static uint saveRestartInterval;
   static uint diagnInterval;

   static std::string solar_wind_file;	/**< Read solar wind data from this file. */

   static bool save_spatial_grid;	/**< Save spatial cell averages for the whole simulation. */
   static bool save_velocity_grid;	/**< Save the velocity grid of every spatial cell in the simulation. */
   static std::vector<Real> save_spatial_cells_x;	/**< Save the velocity grid of spatial cells at these locations. */
   static std::vector<Real> save_spatial_cells_y;
   static std::vector<Real> save_spatial_cells_z;
   
   static uint transmit; /**< Indicates the data that needs to be transmitted to remote nodes.
			  * This is created with bitwise or from the values defined in 
			  * namespace Transmit.*/

   Parameters(int argc, char* argv[]);
};

#endif
