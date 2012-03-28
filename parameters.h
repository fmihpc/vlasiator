/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <limits>
#include <mpi.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "definitions.h"

const uint64_t INVALID_CELLID = 0;

struct Parameters {
   static Real xmin;  /**< X-coordinate of the lower left corner of the spatial grid. */
//   static Real xmax;  /**< Y-coordinate of the lower left corner of the spatial grid. */
   static Real ymin;  /**< Z-coordinate of the lower left corner of the spatial grid. */
   //  static Real ymax;  /**< X-coordinate of the upper right corner of the spatial grid. */
   static Real zmin;  /**< Y-coordinate of the upper right corner of the spatial grid. */
   //  static Real zmax;  /**< Z-coordinate of the upper right corner of the spatial grid. */
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

   static Real q;                    /**< Charge of simulated particle species.*/
   static Real m;                    /**< Mass of simulated particle species.*/
   static Real q_per_m;              /**< Charge-to-mass ratio of simulated particle species,
				      * calculated from Parameters::q and Parameters::m.*/
   static Real t;                    /**< Current simulation time. */
   static Real dt;                   /**< The value of the timestep to use in propagation. */
   static luint tstep_min;           /**< Timestep when simulation starts, needed for restarts.*/
   static luint tstep;               /**< The number of the current timestep. 0=initial state. */
   static luint tsteps;              /**< Total number of timesteps to calculate. */
   static luint saveRestartInterval;
   static luint diagnInterval;

   static bool save_spatial_grid;	/**< Save spatial cell averages for the whole simulation. */
   static bool save_velocity_grid;	/**< Save the velocity grid of every spatial cell in the simulation. */
   
   static uint transmit; /**< Indicates the data that needs to be transmitted to remote nodes.
			  * This is created with bitwise or from the values defined in 
			  * namespace Transmit.*/

   static bool recalculateStencils; /**< If true, MPI stencils should be recalculated because of 
				     * load balancing.*/
   
   static bool propagateField;      /**< If true, magnetic field is propagated during the simulation.*/
   static bool propagateVlasov;     /**< If true, distribution function is propagated during the simulation.*/
   static uint splitMethod;          /**< Split method for splitting spatial/velocity space solvers. 0: first order, 1: strang splitting with half-steps for spatial space, 2: strang splitting with half-steps for velocity space **/
   static bool periodic_x, periodic_y, periodic_z; /**< Whether spatial vlasov grid is periodic */
   static Real sparseMinValue; /**< Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents */
   static Real sparseMinAvgValue; /**< Minimum value of the average of distribution function within a velocity block for the block to be considered to have contents */
   static uint blockAdjustmentInterval; /**< Block adjustment interval (steps). */
   static std::string loadBalanceAlgorithm; /**< Algorithm to be used for load balance.*/
   static std::string loadBalanceTolerance; /**< Load imbalance tolerance. */
   static uint rebalanceInterval; /**< Load rebalance interval (steps). */

   static bool addParameters();
   static bool getParameters();
   
};

#endif
