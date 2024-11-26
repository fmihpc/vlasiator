/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <limits>
#include <mpi.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <map>
#include "fsgrid.hpp"

#include "definitions.h"

const uint64_t INVALID_CELLID = 0;

struct Parameters {
   static int geometry; /**< Simulation geometry, one of the values defined in
                         * geometry::Setup. Defaults to geometry::XYZ6D.*/
   static Real xmin;    /*!< X-coordinate of the lower left corner of the spatial grid. */
   static Real xmax;    /*!< X-coordinate of the upper right corner of the spatial grid. */
   static Real ymin;    /*!< Y-coordinate of the lower left corner of the spatial grid. */
   static Real ymax;    /*!< Y-coordinate of the upper right corner of the spatial grid. */
   static Real zmin;    /*!< Z-coordinate of the lower left corner of the spatial grid. */
   static Real zmax;    /*!< Z-coordinate of the upper right corner of the spatial grid. */
   static Real dx_ini;  /*!< Initial size of spatial cell in x-direction. */
   static Real dy_ini;  /*!< Initial size of spatial cell in y-direction. */
   static Real dz_ini;  /*!< Initial size of spatial cell in z-direction. */

   static uint xcells_ini; /*!< Initial number of spatial cells in x-direction. */
   static uint ycells_ini; /*!< Initial number of spatial cells in y-direction. */
   static uint zcells_ini; /*!< Initial number of spatial cells in z-direction. */

   static Real t;     /*!< Current simulation time. */
   static Real t_min; /*!< Initial simulation time. */
   static Real t_max; /*!< Maximum simulation time. */
   static Real dt;    /*!< The value of the timestep to use in propagation. If CflLimit defined then it is dynamically
                         updated during simulation*/
   static Real vlasovSolverMaxCFL;   /*!< The maximum CFL limit for propagation of distribution function. Used to set
                                        timestep if useCFLlimit is true. */
   static Real vlasovSolverMinCFL;   /*!< The minimum CFL limit for propagation of distribution function. Used to set
                                        timestep if useCFLlimit is true. */
   static bool vlasovSolverGhostTranslate;   /*!< Flag for activating all-local ghost translation. */
   static uint vlasovSolverGhostTranslateExtent;   /*!< Define extent of ghost-translated region in all-local ghost translation. */
   static Real fieldSolverMinCFL;    /*!< The minimum CFL limit for propagation of fields. Used to set timestep if
                                        useCFLlimit is true.*/
   static Real fieldSolverMaxCFL;    /*!< The maximum CFL limit for propagation of fields. Used to set timestep if
                                        useCFLlimit is true.*/
   static uint fieldSolverSubcycles; /*!< The number of field solver subcycles to compute.*/

   static uint tstep_min; /*!< Timestep when simulation starts, needed for restarts.*/
   static uint tstep_max; /*!< Maximum timestep. */
   static uint tstep;     /*!< The number of the current timestep. 0=initial state. */

   static bool meshRepartitioned;         /*!< If true, mesh was repartitioned on this time step.*/
   static std::vector<CellID> localCells; /*!< Cached copy of spatial cell IDs on this process.*/

   static uint diagnosticInterval;
   static std::vector<std::string> systemWriteName;  /*!< Names for the different classes of grid output*/
   static std::vector<std::string> systemWritePath;  /*!< Save this series in this location. Default is ./ */
   static std::vector<Real> systemWriteTimeInterval; /*!< Interval in simusecond for output for each class*/
   static std::vector<int>
       systemWriteDistributionWriteStride; /*!< Every this many cells write out their velocity space in each class. */
   static std::vector<int>
       systemWriteDistributionWriteXlineStride; /*!< Every this many lines of cells along the x direction write out
                                                   their velocity space in each class. */
   static std::vector<int>
       systemWriteDistributionWriteYlineStride; /*!< Every this many lines of cells along the y direction write out
                                                   their velocity space in each class. */
   static std::vector<int>
       systemWriteDistributionWriteZlineStride; /*!< Every this many lines of cells along the z direction write out
                                                   their velocity space in each class. */
   static std::vector<Real>
       systemWriteDistributionWriteShellRadius; /*!< At cells intersecting spheres with those radii centred at the
                                                   origin write out their velocity space in each class. */
   static std::vector<int>
       systemWriteDistributionWriteShellStride; /*!< Every this many cells for those on selected shells write out their
                                                   velocity space in each class. */
   static std::vector<bool> systemWriteFsGrid; /*!< Write fg_ variables in this file class or not.*/
   static bool systemWriteAllDROs; /*!< Write all output DROs or not.*/
   static bool diagnosticWriteAllDROs; /*!< Write all diagnostic DROs or not.*/
   static std::vector<int> systemWrites;        /*!< How many files have been written of each class*/
   static std::vector<std::pair<std::string, std::string>>
       systemWriteHints; /*!< Collection of MPI-IO hints passed for non-restart IO. Pairs of key-value strings. */
   static std::vector<std::pair<std::string, std::string>>
       restartReadHints; /*!< Collection of MPI-IO hints passed for restart IO. Pairs of key-value strings. */
   static std::vector<std::pair<std::string, std::string>>
       restartWriteHints; /*!< Collection of MPI-IO hints passed for restart IO. Pairs of key-value strings. */

   static bool writeInitialState; /*!< If true, initial state is written. This is useful for debugging as the restarts
                                     are always written out after propagation of 0.5dt in real space.*/
   static bool writeFullBGB; /*!< If true, write full BGB components and derivatives in a dedicated file, then exit.*/
   static Real saveRestartWalltimeInterval; /*!< Interval in walltime seconds for restart data*/
   static uint exitAfterRestarts;           /*!< Exit after this many restarts*/
   static uint64_t vlsvBufferSize;          /*!< Buffer size in bytes passed to VLSV writer. */
   static int restartStripeFactor;          /*!< stripe_factor for restart writing*/
   static int systemStripeFactor;             /*!< stripe_factor for bulk and initial grid writing*/
   static std::string restartWritePath; /*!< Path to the location where restart files should be written. Defaults to the
                                           local directory, also if the specified destination is not writeable. */

   static uint transmit;
   /*!< Indicates the data that needs to be transmitted to remote nodes.
    * This is created with bitwise or from the values defined in
    * namespace Transmit.*/

   static bool recalculateStencils; /*!< If true, MPI stencils should be recalculated because of load balancing.*/

   static bool propagateField;              /*!< If true, magnetic field is propagated during the simulation.*/
   static bool propagateVlasovAcceleration; /*!< If true, distribution function is propagated in velocity space during
                                               the simulation.*/
   static bool propagateVlasovTranslation;  /*!< If true, distribution function is propagated in ordinary space during
                                               the simulation.*/

   static Real maxWaveVelocity;         /*!< Maximum wave velocity allowed in LDZ. */
   static uint maxFieldSolverSubcycles; /*!< Maximum allowed field solver subcycles. */
   static Real resistivity;             /*!< Resistivity in Ohm's law eta*J term. */
   static uint ohmHallTerm; /*!< Enable/choose spatial order of Hall term in Ohm's law JXB term. 0: off, 1: 1st spatial
                               order, 2: 2nd spatial order. */
   static uint ohmGradPeTerm; /*!< Enable/choose spatial order of the electron pressure gradient term in Ohm's law. 0:
                                 off, 1: 1st spatial order. */
   static Real electronTemperature; /*!< Upstream electron temperature to be used for the electron pressure gradient
                                       term (K). */
   static Real
       electronDensity; /*!< Upstream electron density to be used for the electron pressure gradient term (m^-3). */
   static Real electronPTindex; /*!> Polytropic index for electron pressure gradient term. 0 is isobaric, 1 is
                                   isothermal, 1.667 is adiabatic electrons */

   static bool fieldSolverDiffusiveEterms; /*!< Enable resistive terms in the computation of E*/

   static Real maxSlAccelerationRotation; /*!< Maximum rotation in acceleration for semilagrangian solver*/
   static int maxSlAccelerationSubcycles; /*!< Maximum number of subcycles in acceleration*/
   static bool vlasovAccelerateMaxwellianBoundaries; /*!< Accelerate also Maxwellian boundary cells*/

   static Real hallMinimumRhom; /*!< Minimum mass density value used in the field solver.*/
   static Real hallMinimumRhoq; /*!< Minimum charge density value used for the Hall and electron pressure gradient terms
                                   in the Lorentz force and in the field solver.*/

   static std::string loadBalanceAlgorithm; /*!< Algorithm to be used for load balance.*/
   static std::map<std::string, std::string> loadBalanceOptions;  // Other Load balancing options
   static uint rebalanceInterval;           /*!< Load rebalance interval (steps). */
   static bool prepareForRebalance; /**< If true, propagators should measure their time consumption in preparation
                                     * for mesh repartitioning.*/

   static std::vector<std::string>
       outputVariableList; /*!< List of data reduction operators (DROs) to add to the grid file output.*/
   static std::vector<std::string>
       diagnosticVariableList; /*!< List of data reduction operators (DROs) to add to the diagnostic runtime output.*/

   static std::string restartFileName; /*!< If defined, restart from this file*/
   static bool isRestart;              /*!< true if this is a restart, false otherwise */
   static int writeAsFloat;            /*!< true if writing into VLSV in floats instead of doubles, false otherwise */
   static int
       writeRestartAsFloat;     /*!< true if writing into restart files in floats instead of doubles, false otherwise */
   static bool dynamicTimestep; /*!< If true, timestep is set based on  CFL limit */

   static std::string projectName; /*!< Project to be used in this run. */

   static bool bailout_write_restart; /*!< If true, write a restart file on bailout. Gets reset when sending a STOP
                                         (true) or a KILL (false). */
   static Real bailout_min_dt;        /*!< Minimum time step below which bailout occurs (s). */
   static Real bailout_max_memory;    /*!< Maximum amount of memory used per node (in GiB) over which bailout occurs. */
   static uint bailout_velocity_space_wall_margin; /*!< Safety margin in number of blocks off the v-space wall beyond which bailout occurs. */

   static uint vamrMaxVelocityRefLevel; /**< Maximum velocity mesh refinement level, defaults to 0.*/
   static Realf vamrCoarsenLimit; /**< If the value of refinement criterion is below this value, block can be coarsened.
                                  * The value must be smaller than vamrRefineLimit.*/
   static Realf vamrRefineLimit;  /**< If the value of refinement criterion is larger than this value, block should be
                                  * refined.  The value must be larger than vamrCoarsenLimit.*/
   static std::string vamrVelRefCriterion; /**< Name of the velocity block refinement criterion function.*/

   static int amrMaxSpatialRefLevel; /*!< Absolute maximum refinement level (conditions the fsgrid resolution), cannot be exceeded after initial setup of the grids. */
   static int amrMaxAllowedSpatialRefLevel; /*!< Maximum currently allowed refinement level for restart or dynamic refinement. */
   static bool adaptRefinement;
   static bool refineOnRestart;
   static bool forceRefinement;
   static bool shouldFilter;
   static bool useAlpha1;
   static Real alpha1RefineThreshold;
   static Real alpha1CoarsenThreshold;
   static bool useAlpha2;
   static Real alpha2RefineThreshold;
   static Real alpha2CoarsenThreshold;
   // TODO: consider renaming to alpha3 or something to that effect
   static bool useVorticity;
   static Real vorticityRefineThreshold;
   static Real vorticityCoarsenThreshold;
   static bool useAnisotropy;
   static Real anisotropyRefineThreshold;
   static Real anisotropyCoarsenThreshold;
   static int anisotropyMaxReflevel;
   static uint refineCadence;
   static Real refineAfter;
   static Real refineRadius;
   static Real alphaDRhoWeight;
   static Real alphaDUWeight;
   static Real alphaDPSqWeight;
   static Real alphaDBSqWeight;
   static Real alphaDBWeight;
   static Real refinementMinX; /*!< Do not refine at x coordinates below this value. */
   static Real refinementMinY; /*!< Do not refine at y coordinates below this value. */
   static Real refinementMinZ; /*!< Do not refine at z coordinates below this value. */
   static Real refinementMaxX; /*!< Do not refine at x coordinates above this value. */
   static Real refinementMaxY; /*!< Do not refine at y coordinates above this value. */
   static Real refinementMaxZ; /*!< Do not refine at z coordinates above this value. */
   static int maxFilteringPasses;
   static int amrBoxNumber;
   static std::vector<uint> amrBoxHalfWidthX;
   static std::vector<uint> amrBoxHalfWidthY;
   static std::vector<uint> amrBoxHalfWidthZ;
   static std::vector<Realf> amrBoxCenterX;
   static std::vector<Realf> amrBoxCenterY;
   static std::vector<Realf> amrBoxCenterZ;
   static std::vector<int> amrBoxMaxLevel;
   static bool amrTransShortPencils;        /*!< Use short or long pencils in AMR translation.*/
   static std::vector<std::string> blurPassString;
   static std::vector<int> numPasses;

   static std::array<FsGridTools::Task_t,3> manualFsGridDecomposition;
   static std::array<FsGridTools::Task_t,3> overrideReadFsGridDecomposition;
   
   static bool computeCurvature; /*<! Boolean flag, if true the curvature of magnetic field is computed. */

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
    * \param firstPass determines whether to only parse out particle populations (first pass),
    *                  to parse all other parameters.
    * \sa addParameters
    */
   static void getParameters();
};

#endif
