/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
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
   static int geometry; /**< Simulation geometry, one of the values defined in 
                         * geometry::Setup. Defaults to geometry::XYZ6D.*/
   static Real xmin;  /*!< X-coordinate of the lower left corner of the spatial grid. */
   static Real xmax;  /*!< X-coordinate of the upper right corner of the spatial grid. */
   static Real ymin;  /*!< Y-coordinate of the lower left corner of the spatial grid. */
   static Real ymax;  /*!< Y-coordinate of the upper right corner of the spatial grid. */
   static Real zmin;  /*!< Z-coordinate of the lower left corner of the spatial grid. */
   static Real zmax;  /*!< Z-coordinate of the upper right corner of the spatial grid. */
   static Real dx_ini; /*!< Initial size of spatial cell in x-direction. */
   static Real dy_ini; /*!< Initial size of spatial cell in y-direction. */
   static Real dz_ini; /*!< Initial size of spatial cell in z-direction. */

   static uint xcells_ini; /*!< Initial number of spatial cells in x-direction. */
   static uint ycells_ini; /*!< Initial number of spatial cells in y-direction. */
   static uint zcells_ini; /*!< Initial number of spatial cells in z-direction. */

   static Real backstreamradius; /*!< Radius of the maxwellian distribution. Used for calculating rho of the backstream population. */
   static Real backstreamvx; /*!< X coordinate of the origin of the maxwellian distribution. Used for calculating rho of the backstream population. */
   static Real backstreamvy; /*!< Y coordinate of the origin of the maxwellian distribution. Used for calculating rho of the backstream population. */
   static Real backstreamvz; /*!< Z coordinate of the origin of the maxwellian distribution. Used for calculating rho of the backstream population. */
   
   static Real t;                    /*!< Current simulation time. */
   static Real t_min;                    /*!< Initial simulation time. */
   static Real t_max;                    /*!< Maximum simulation time. */
   static Real dt;                   /*!< The value of the timestep to use in propagation. If CflLimit defined then it is dynamically updated during simulation*/
   static Real vlasovSolverMaxCFL;   /*!< The maximum CFL limit for propagation of distribution function. Used to set timestep if useCFLlimit is true. */
   static Real vlasovSolverMinCFL;   /*!< The minimum CFL limit for propagation of distribution function. Used to set timestep if useCFLlimit is true. */
   static Real fieldSolverMinCFL;     /*!< The minimum CFL limit for propagation of fields. Used to set timestep if useCFLlimit is true.*/
   static Real fieldSolverMaxCFL;     /*!< The maximum CFL limit for propagation of fields. Used to set timestep if useCFLlimit is true.*/
   static int fieldSolverSubcycles;     /*!< The number of field solver subcycles to compute.*/

   static uint tstep_min;           /*!< Timestep when simulation starts, needed for restarts.*/
   static uint tstep_max;           /*!< Maximum timestep. */
   static uint tstep;               /*!< The number of the current timestep. 0=initial state. */

   static bool meshRepartitioned;         /*!< If true, mesh was repartitioned on this time step.*/
   static std::vector<CellID> localCells; /*!< Cached copy of spatial cell IDs on this process.*/

   static uint diagnosticInterval;
   static std::vector<std::string> systemWriteName; /*!< Names for the different classes of grid output*/
   static std::vector<std::string> systemWritePath; /*!< Save this series in this location. Default is ./ */
   static std::vector<Real> systemWriteTimeInterval;/*!< Interval in simusecond for output for each class*/
   static std::vector<int> systemWriteDistributionWriteStride; /*!< Every this many cells write out their velocity space in each class. */
   static std::vector<int> systemWriteDistributionWriteXlineStride; /*!< Every this many lines of cells along the x direction write out their velocity space in each class. */
   static std::vector<int> systemWriteDistributionWriteYlineStride; /*!< Every this many lines of cells along the y direction write out their velocity space in each class. */
   static std::vector<int> systemWriteDistributionWriteZlineStride; /*!< Every this many lines of cells along the z direction write out their velocity space in each class. */
   static std::vector<int> systemWrites; /*!< How many files have been written of each class*/
   
   static bool writeInitialState;           /*!< If true, initial state is written. This is useful for debugging as the restarts are always written out after propagation of 0.5dt in real space.*/
   static Real saveRestartWalltimeInterval; /*!< Interval in walltime seconds for restart data*/
   static uint exitAfterRestarts;           /*!< Exit after this many restarts*/
   static int restartStripeFactor;          /*!< stripe_factor for restart writing*/
   static std::string restartWritePath;          /*!< Path to the location where restart files should be written. Defaults to the local directory, also if the specified destination is not writeable. */
   
   static uint transmit;
   /*!< Indicates the data that needs to be transmitted to remote nodes.
    * This is created with bitwise or from the values defined in 
    * namespace Transmit.*/
   
   static bool recalculateStencils; /*!< If true, MPI stencils should be recalculated because of load balancing.*/
   
   static bool propagateField;      /*!< If true, magnetic field is propagated during the simulation.*/
   static bool propagatePotential;  /*!< If true, electrostatic potential is solved during the simulation.*/
   static bool propagateVlasovAcceleration;     /*!< If true, distribution function is propagated in velocity space during the simulation.*/
   static bool propagateVlasovTranslation;      /*!< If true, distribution function is propagated in ordinary space during the simulation.*/

   static Real maxWaveVelocity; /*!< Maximum wave velocity allowed in LDZ. */
   static int maxFieldSolverSubcycles; /*!< Maximum allowed field solver subcycles. */
   static Real resistivity; /*!< Resistivity in Ohm's law eta*J term. */
   static uint ohmHallTerm; /*!< Enable/choose spatial order of Hall term in Ohm's law JXB term. 0: off, 1: 1st spatial order, 2: 2nd spatial order. */
   static uint ohmGradPeTerm; /*!< Enable/choose spatial order of the electron pressure gradient term in Ohm's law. 0: off, 1: 1st spatial order. */
   static Real electronTemperature; /*!< Constant electron temperature to be used for the electron pressure gradient term (K). */
   static bool fieldSolverDiffusiveEterms; /*!< Enable resistive terms in the computation of E*/
   
   static Real maxSlAccelerationRotation; /*!< Maximum rotation in acceleration for semilagrangian solver*/
   static int maxSlAccelerationSubcycles; /*!< Maximum number of subcycles in acceleration*/

   static Real hallMinimumRho;  /*!< Minimum rho value used for the Hall and electron pressure gradient terms in the Lorentz force and in the field solver.*/
   static Real sparseMinValue; /*!< (DEPRECATED) Minimum value of distribution function in any cell of a velocity 
                                * block for the block to be considered to have content.
                                * This value is only used for default particle species.*/
   static int sparseBlockAddWidthV; /*!< Number of layers of blocks that are kept in velocity space around the blocks with content */
   static bool sparse_conserve_mass; /*!< If true, density is scaled to conserve mass when removing blocks*/
   static int  sparseDynamicAlgorithm; /*!< Type of algorithm used for calculating the dynamic minValue; 0 = none, 1 = linear algorithm based on minValue and rho, 2 = linear algorithm based on minValue and Blocks, (Example linear algorithm: minValue = rho / sparse.dynamicValue * sparse.minValue)*/
   static Real sparseDynamicBulkValue1; /*!< Minimum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicMinValue = 1e3, sparse.dynamicMaxValue=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5*/
   static Real sparseDynamicBulkValue2; /*!< Maximum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicMinValue = 1e3, sparse.dynamicMaxValue=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5*/
   static Real sparseDynamicMinValue1; /*!< The minimum value for the minValue*/
   static Real sparseDynamicMinValue2; /*!< The maximum value for the minValue*/
   
   static std::string loadBalanceAlgorithm; /*!< Algorithm to be used for load balance.*/
   static std::string loadBalanceTolerance; /*!< Load imbalance tolerance. */ 
   static uint rebalanceInterval; /*!< Load rebalance interval (steps). */
   static bool prepareForRebalance; /**< If true, propagators should measure their time consumption in preparation
                                     * for mesh repartitioning.*/
   
   static std::vector<std::string> outputVariableList; /*!< List of data reduction operators (DROs) to add to the grid file output.*/
   static std::vector<std::string> diagnosticVariableList; /*!< List of data reduction operators (DROs) to add to the diagnostic runtime output.*/
   
   static std::string restartFileName; /*!< If defined, restart from this file*/
   static bool isRestart; /*!< true if this is a restart, false otherwise */
   static int writeAsFloat; /*!< true if writing into VLSV in floats instead of doubles, false otherwise */
   static bool dynamicTimestep; /*!< If true, timestep is set based on  CFL limit */
   
   static std::string projectName; /*!< Project to be used in this run. */
   
   static bool bailout_write_restart; /*!< If true, write a restart file on bailout. Gets reset when sending a STOP (true) or a KILL (false). */
   static Real bailout_min_dt; /*!< Minimum time step below which bailout occurs (s). */

   static uint amrMaxVelocityRefLevel;    /**< Maximum velocity mesh refinement level, defaults to 0.*/
   static Realf amrCoarsenLimit;          /**< If the value of refinement criterion is below this value, block can be coarsened.
                                           * The value must be smaller than amrRefineLimit.*/
   static Realf amrRefineLimit;           /**< If the value of refinement criterion is larger than this value, block should be refined.
                                           * The value must be larger than amrCoarsenLimit.*/
   static std::string amrVelRefCriterion; /**< Name of the velocity block refinement criterion function.*/

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
