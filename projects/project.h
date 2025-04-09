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
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef PROJECT_H
#define PROJECT_H

#include <random>
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "fsgrid.hpp"

namespace projects {

   /** Returns the phase-space density of a Maxwellian distribution function
    * NOTE: This function is called inside parallel region so it must be declared as const.
    * @param vx The vx-coordinate relative to the Maxwellian centre
    * @param vy The vy-coordinate relative to the Maxwellian centre
    * @param vz The vz-coordinate relative to the Maxwellian centre
    * @param T The Maxwellian temperature
    * @param rho The total number density of the Maxwellian distribution
    * @param mass The mass of the particle in question
    * @return The distribution function value at the given velocity coordinates
    * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
    */
   ARCH_HOSTDEV inline Realf MaxwellianPhaseSpaceDensity(
      creal& vx, creal& vy, creal& vz,
      creal& T, creal& rho, creal& mass
      ) {
      return rho * pow(mass / (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
         exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * T));
   }

   /** Returns the phase-space density of a Tri-Maxwellian distribution function
    * NOTE: This function is called inside parallel region so it must be declared as const.
    * @param vx The vx-coordinate relative to the Tri-Maxwellian centre
    * @param vy The vy-coordinate relative to the Tri-Maxwellian centre
    * @param vz The vz-coordinate relative to the Tri-Maxwellian centre
    * @param Tx The Tri-Maxwellian x-directional temperature
    * @param Ty The Tri-Maxwellian y-directional temperature
    * @param Tz The Tri-Maxwellian z-directional temperature
    * @param rho The total number density of the Tri-Maxwellian distribution
    * @param mass The mass of the particle in question
    * @return The distribution function value at the given velocity coordinates
    * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
    */
   ARCH_HOSTDEV inline Realf TriMaxwellianPhaseSpaceDensity(
      creal& vx, creal& vy, creal& vz,
      creal& Tx, creal& Ty, creal& Tz,
      creal& rho, creal& mass) {
      return rho * pow(mass / (2.0 * M_PI * physicalconstants::K_B), 1.5) *
         exp(- mass * (vx*vx/Tx + vy*vy/Ty + vz*vz/Tz) / (2.0 * physicalconstants::K_B)) /
         sqrt(Tx*Ty*Tz);
   }

   class Project {
    public:
      Project();
      virtual ~Project();
      
      /*! Register parameters that should be read in. */
      static void addParameters();
      
      virtual Real getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const uint popID) const;
      
      /*! Get the value that was read in. */
      virtual void getParameters();
      
      /*! Initialize project. Can be used, e.g., to read in parameters from the input file. */
      virtual bool initialize();
      
      /*! Perform some operation at each time step in the main program loop. */
      virtual void hook(
         cuint& stage,
         const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid
      ) const;
      
      bool initialized();
      
      /** Set the background and perturbed magnetic fields for this project.
       * \param perBGrid Grid on which values of the perturbed field can be set if needed.
       * \param BgBGrid Grid on which values for the background field can be set if needed, e.g. using the background field functions.
       * \param technicalGrid Technical fsgrid, available if some of its data is necessary.
       * 
       * \sa setBackgroundField, setBackgroundFieldToZero
       */
      virtual void setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
      );
      
      /*! Setup data structures for subsequent setCell calls.
       * This will most likely be empty for most projects, except for some advanced
       * data juggling ones (like restart from a subset of a larger run)
       * \param cells Local cellIDs of this task.
       */
      virtual void setupBeforeSetCell(const std::vector<CellID>& cells);

      /*!\brief Set the perturbed fields and distribution of a cell according to the default simulation settings.
       * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
       * \param cell Pointer to the cell to set.
       */
      void setCell(spatial_cell::SpatialCell* cell);
         
      virtual bool canRefine(spatial_cell::SpatialCell* cell) const;

      virtual bool shouldRefineCell(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID id, Real r_max2) const;

      virtual bool shouldUnrefineCell(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID id, Real r_max2) const;

      virtual bool refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const;

      /*!\brief Adapts refinement by one level according to the project. Returns true if any cells were refined, false if not.
       * \param mpiGrid grid to refine
       * @return The amount of cells set to refine
       */
      virtual uint64_t adaptRefinement( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const;


      /*!\brief Refine/unrefine spatial cells one level to the static criteria in the config
      * \param mpiGrid Spatial grid
      * \param n Static refinement pass. 0th pass refines level 0 cells and unrefines max level cells, 1st pass refines level 1 and unrefines max level -1 etc.
      */
      virtual bool forceRefinement( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, int n ) const;

      /*!\brief Boxcar filters spatial cells that were recently refined
       * \param mpiGrid grid to filter
       */
      virtual bool filterRefined( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const;

    protected:
      /*! \brief Prepares a  list of blocks to loop through when initialising.
       * 
       * The base class version just prepares all blocks, which amounts to looping through the whole velocity space.
       * This is very expensive and becomes prohibitive in cases where a large velocity space is needed with only
       * small portions actually containing something. Use with care.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * The function stores the prepared blocks into cell->velocity_block_with_content_list and returns the count.
       */
      virtual uint findBlocksToInitialize(spatial_cell::SpatialCell* cell,const uint popID) const;
      
      /*! \brief Sets the distribution function in a cell.
       * 
       * Uses the function findBlocksToInitialize and loops through the list returned by it to initialize the cells' velocity space.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * 
       * \sa findBlocksToInitialize
       */
      void setVelocitySpace(const uint popID,spatial_cell::SpatialCell* cell) const;
         
      /** Calculate potentially needed parameters for the given spatial cell at the given time.
       * 
       * Currently this function is only called during initialization.
       * 
       * The following array indices contain the coordinates of the "lower left corner" of the cell: 
       * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
       * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
       * @param cell Pointer to the spatial cell to be handled.
       * @param t The current value of time. This is passed as a convenience. If you need more detailed information 
       * of the state of the simulation, you can read it from Parameters.
       */
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);

      /** Calculates the distribution function contents for the spatial cell and population in question
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * This function will contain a loop over nRequested velocity blocks,
       * with the block GlobalIDs provided in the GIDlist buffer,
       * storing phase-space densities in the bufferData buffer.
       *
       * @param popID Particle species ID.
       * @return The total stored number density of the cell for this population
       * The physical unit of this quantity is 1/m^3.
       */
      virtual Realf fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                  const uint popID,
                                  const uint nRequested) const = 0;

      void printPopulations();
      
      virtual bool rescalesDensity(const uint popID) const;
      void rescaleDensity(spatial_cell::SpatialCell* cell,const uint popID) const;
      
      /** Get random number between 0 and 1.0. One should always first initialize the rng.
       * @param rngDataBuffer struct of type random_data
       * @return Uniformly distributed random number between 0 and 1.*/
      Real getRandomNumber(std::default_random_engine& randGen) const;

      /** Set random seed (thread-safe). Seed is based on the seed read
       *  in from cfg + the seedModifier parameter
       * @param seedModifier CellID value to use as seed modifier
       * @param rngStateBuffer buffer where random number values are kept
       * @param rngDataBuffer struct of type random_data
       */
      void setRandomSeed(uint64_t seedModifier, std::default_random_engine& randGen) const;

      /** Set random seed (thread-safe) that is always the same for
       * this particular cellID. Can be used to make reproducible
       * simulations that do not depend on number of processes or threads.
       * @param cell SpatialCell used to infer CellID value to use as seed modifier
       * @param rngStateBuffer buffer where random number values are kept
       * @param rngDataBuffer struct of type random_data
       */
      void setRandomCellSeed(spatial_cell::SpatialCell* cell, std::default_random_engine& randGen) const;
      
    private:
      uint seed;
      static char rngStateBuffer[256];

      bool baseClassInitialized;                      /**< If true, base class has been initialized.*/
   };
   
   Project* createProject();
} // namespace projects


#endif

