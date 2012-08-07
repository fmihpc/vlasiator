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
   static Real xmin;  /*!< X-coordinate of the lower left corner of the spatial grid. */
   static Real xmax;  /*!< X-coordinate of the upper right corner of the spatial grid. */
   static Real ymin;  /*!< Y-coordinate of the lower left corner of the spatial grid. */
   static Real ymax;  /*!< Y-coordinate of the upper right corner of the spatial grid. */
   static Real zmin;  /*!< Z-coordinate of the lower left corner of the spatial grid. */
   static Real zmax;  /*!< Z-coordinate of the upper right corner of the spatial grid. */
   static Real dx_ini; /*!< Initial size of spatial cell in x-direction. */
   static Real dy_ini; /*!< Initial size of spatial cell in y-direction. */
   static Real dz_ini; /*!< Initial size of spatial cell in z-direction. */
   
   static Real vxmin; /*!< VX-coordinate of the lower left corner of velocity grid. */
   static Real vxmax; /*!< VY-coordinate of the lower left corner of velocity grid. */
   static Real vymin; /*!< VZ-coordinate of the lower left corner of velocity grid. */
   static Real vymax; /*!< VX-coordinate of the upper right corner of velocity grid. */
   static Real vzmin; /*!< VY-coordinate of the upper right corner of velocity grid. */
   static Real vzmax; /*!< VZ-coordinate of the upper right corner of velocity grid. */
   
   static uint xcells_ini; /*!< Initial number of spatial cells in x-direction. */
   static uint ycells_ini; /*!< Initial number of spatial cells in y-direction. */
   static uint zcells_ini; /*!< Initial number of spatial cells in z-direction. */
   static uint vxblocks_ini; /*!< Initial number of velocity grid blocks in vx-direction. */
   static uint vyblocks_ini; /*!< Initial number of velocity grid blocks in vy-direction. */
   static uint vzblocks_ini; /*!< Initial number of velocity grid blocks in vz-direction. */

   static Real q;                    /*!< Charge of simulated particle species.*/
   static Real m;                    /*!< Mass of simulated particle species.*/
   static Real q_per_m;              /*!< Charge-to-mass ratio of simulated particle species,
				      * calculated from Parameters::q and Parameters::m.*/
   static Real t;                    /*!< Current simulation time. */
   static Real t_min;                    /*!< Initial simulation time. */
   static Real t_max;                    /*!< Maximum simulation time. */
   static Real dt;                   /*!< The value of the timestep to use in propagation. If CflLimit defined then it is dynamically updated during simulation*/
   static Real CFL;                  /*!< The maximum CFL limit for propagation. Used to set timestep if useCFLlimit is true. Also used to set number of acceleration steps if substepAcceleration is true */
   
   static luint tstep_min;           /*!< Timestep when simulation starts, needed for restarts.*/
   static luint tstep_max;           /*!< Maximum timestep. */
   static luint tstep;               /*!< The number of the current timestep. 0=initial state. */

   static luint diagnosticInterval;
   static Real saveRestartTimeInterval;
   static Real saveSystemTimeInterval;
   
   
   static uint transmit; /*!< Indicates the data that needs to be transmitted to remote nodes.
			  * This is created with bitwise or from the values defined in 
			  * namespace Transmit.*/

   static bool recalculateStencils; /*!< If true, MPI stencils should be recalculated because of 
				     * load balancing.*/
   
   static bool propagateField;      /*!< If true, magnetic field is propagated during the simulation.*/
   static bool propagateVlasov;     /*!< If true, distribution function is propagated during the simulation.*/
   static uint splitMethod;          /*!< Split method for splitting spatial/velocity space solvers. 0: first order, 1: strang splitting with half-steps for spatial space, 2: strang splitting with half-steps for velocity space **/
   static bool periodic_x, periodic_y, periodic_z; /*!< Whether spatial vlasov grid is periodic */
   
   static Real RK_alpha; /*!< Parameter of the second-order Runge-Kutta method employed in the field solver. **/
   
   static Real sparseMinValue; /*!< Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents */
   static Real sparseMinAvgValue; /*!< Minimum value of the average of distribution function within a velocity block for the block to be considered to have contents */
   static uint blockAdjustmentInterval; /*!< Block adjustment interval (steps). */
   static std::string loadBalanceAlgorithm; /*!< Algorithm to be used for load balance.*/
   static std::string loadBalanceTolerance; /*!< Load imbalance tolerance. */
   static uint rebalanceInterval; /*!< Load rebalance interval (steps). */
   
   static std::vector<std::string> outputVariableList; /*!< List of data reduction operators (DROs) to add to the grid file output.*/
   static std::vector<std::string> diagnosticVariableList; /*!< List of data reduction operators (DROs) to add to the diagnostic runtime output.*/
   
   static std::vector<std::string> boundaryCondList; /*!< List of boundary conditions (BC) to be used. */
   static std::vector<std::string> outflowFaceList; /*!< List of faces on which outflow boundary conditions are to be applied ([+-][xyz]). */
   static std::vector<std::string> solarWindFaceList; /*!< List of faces on which solar wind boundary conditions are to be applied ([+-][xyz]). */
   static std::string solarWindFiles[6]; /*!< Input files for the solar wind boundary conditions. */
   static bool isSolarWindDynamic; /*!< Is the solar wind inflow dynamic in time or not. */
   static Real ionoCenter[3]; /*!< Coordinates of the centre of the ionosphere. */
   static Real ionoRadius; /*!< Radius of the ionosphere. */
   
   static uint maxAccelerationSubsteps; /*!< Maximum number of substeps that is allowed */
   static bool dynamicTimestep; /*!< If true, timestep is set based on  CFL limit */
                                       
   
   static bool addParameters();
   static bool getParameters();
   
};

#endif
