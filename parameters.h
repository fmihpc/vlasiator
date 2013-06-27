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
                                      *   calculated from Parameters::q and Parameters::m.*/
   static Real t;                    /*!< Current simulation time. */
   static Real t_min;                    /*!< Initial simulation time. */
   static Real t_max;                    /*!< Maximum simulation time. */
   static Real dt;                   /*!< The value of the timestep to use in propagation. If CflLimit defined then it is dynamically updated during simulation*/
   static Real vlasovSolverMaxCFL;                  /*!< The maximum CFL limit for propagation of distribution function. Used to set timestep if useCFLlimit is true. Also used to set number of acceleration steps if substepAcceleration is true */
   static Real vlasovSolverMinCFL;                  /*!< The minimum CFL limit for propagation of distribution function. Used to set timestep if useCFLlimit is true. Also used to set number of acceleration steps if substepAcceleration is true */
   static Real fieldSolverMinCFL;     /*!< The minimum CFL limit for propagation of fields. Used to set timestep if useCFLlimit is true.*/
   static Real fieldSolverMaxCFL;     /*!< The maximum CFL limit for propagation of fields. Used to set timestep if useCFLlimit is true.*/

   static luint tstep_min;           /*!< Timestep when simulation starts, needed for restarts.*/
   static luint tstep_max;           /*!< Maximum timestep. */
   static luint tstep;               /*!< The number of the current timestep. 0=initial state. */

   static luint diagnosticInterval;
   static std::vector<std::string> systemWriteName; /*!< Names for the different classes of grid output*/
   static std::vector<Real> systemWriteTimeInterval;/*!< Interval in simusecond for output for each class*/
   static std::vector<int> systemWriteDistributionWriteStride; /*!< Every this many cells write out their velocity space in each class. */
   static std::vector<int> systemWriteDistributionWriteXlineStride; /*!< Every this many lines of cells along the x direction write out their velocity space in each class. */
   static std::vector<int> systemWriteDistributionWriteYlineStride; /*!< Every this many lines of cells along the y direction write out their velocity space in each class. */
   static std::vector<int> systemWriteDistributionWriteZlineStride; /*!< Every this many lines of cells along the z direction write out their velocity space in each class. */
   static std::vector<int> systemWrites; /*!< How many files have been written of each class*/

   static bool writeInitialState;           /*!< If true, initial state is written. This is useful for debugging as the restarts are always written out after propagation of 0.5dt in real space.*/
   static Real saveRestartWalltimeInterval; /*!< Interval in walltime seconds for restart data*/
   static uint exitAfterRestarts;           /*!< Exit after this many restarts*/
   static int restartStripeFactor;           /*!< stripe_factor for restart writing*/
   
   static uint transmit;
   /*!< Indicates the data that needs to be transmitted to remote nodes.
    * This is created with bitwise or from the values defined in 
    * namespace Transmit.*/

   static bool recalculateStencils; /*!< If true, MPI stencils should be recalculated because of load balancing.*/
   
   static bool propagateField;      /*!< If true, magnetic field is propagated during the simulation.*/
   static bool propagateVlasov;     /*!< If true, distribution function is propagated during the simulation.*/
   static bool periodic_x, periodic_y, periodic_z; /*!< Whether spatial vlasov grid is periodic */
   
   static Real maxAlfvenVelocity; /*!< Maximum Alfven velocity allowed in fastMS computation in LDZ. */
   static Real resistivity; /*!< Resistivity in Ohm's law eta*J term. */
   static bool ohmHallTerm; /*!< Hall term in Ohm's law JXB term. */
   static bool fieldSolverDiffusiveEterms; /*!< Enable resitive terms in the computation of E*/

   static bool vlasovSemiLagAcceleration; /*!< If true, a semilagrangian solver is used for acceleration, otherwise we use a the Leveque solver*/
   static bool semiLagAccFrame; /*!< Do we keep an additional dist func for acceleratio and map spatial changes to it*/
   static bool lorentzHallTerm;  /*!< Add JxB hall term to Lorentz force. Assumes there is not hall term in the fieldsolver if we use fieldSolverE*/
   static Real lorentzHallMinimumRho;  /*!< Minimum rho value used in Hall term in Lorentz force.*/
   static Real lorentzHallMaximumB;  /*!< Maximum B value used in Hall term in Lorentz force. Does not affect J, only latter total B in JxB*/
   static bool lorentzUseFieldSolverE;  /*!< If true, the E from fieldolver is used. Otherwise -V x B_vol is used*/
   static Real sparseMinValue; /*!< Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents */
   static int sparseBlockAddWidthV; /*!< Number of layers of blocks that are kept in velocity space around the blocks with content */
   static uint blockAdjustmentInterval; /*!< Block adjustment interval (steps). */
   static std::string loadBalanceAlgorithm; /*!< Algorithm to be used for load balance.*/
   static std::string loadBalanceTolerance; /*!< Load imbalance tolerance. */ 
   static uint rebalanceInterval; /*!< Load rebalance interval (steps). */
   static Real loadBalanceAlpha; /*!< alpha in cell load balance weight = blocks * (alpha + beta*substeps) */
   static Real loadBalanceBeta; /*!< alpha in cell load balance weight = blocks * (alpha + beta*substeps) */
   
   static std::vector<std::string> outputVariableList; /*!< List of data reduction operators (DROs) to add to the grid file output.*/
   static std::vector<std::string> diagnosticVariableList; /*!< List of data reduction operators (DROs) to add to the diagnostic runtime output.*/
   
   static std::string restartFileName; /*!< If defined, restart from this file*/
   static bool isRestart; /*!< true if this is a restart, false otherwise */
   static uint maxAccelerationSubsteps; /*!< Maximum number of substeps that is allowed */
   static bool dynamicTimestep; /*!< If true, timestep is set based on  CFL limit */
   
   
   static std::string projectName; /*!< Project to be used in this run. */
   
   /*! \brief Add the global parameters.
    * 
    * This function adds all the parameters that are loaded at a global level.
    * More are being loaded e.g. in the projects and in the system boundary conditions.
    * 
    * Note that due to the large number of parameters added here, no code is added to check
    * for consistency when they are read later. Please make sure when coding new parameters
    * here that the options in getParameters match the ones added here.
    * 
    * \sa getParameters
    */
   static bool addParameters();
   
   /*! \brief Get the global parameters.
    * 
    * This function gets all the parameters loaded at a global level.
    * More are being loaded e.g. in the projects and in the system boundary conditions.
    * 
    * Note that due to the large number of parameters read here, no code is added to check
    * for consistency with the loaded options, or the code here would become much less
    * readable. Please make sure when coding new parameters here that the options in
    * addParameters match the ones read here.
    * 
    * \sa addParameters
    */
   static bool getParameters();
};

#endif
