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
#include "../spatial_cell_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "fsgrid.hpp"

namespace projects {
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
         
      Real setVelocityBlock(spatial_cell::SpatialCell* cell,const vmesh::LocalID& blockLID,const uint popID) const;

      virtual bool canRefine(spatial_cell::SpatialCell* cell) const;

      virtual bool shouldRefineCell(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID id, Real r_max2) const;

      virtual bool shouldUnrefineCell(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID id, Real r_max2) const;

      virtual bool refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const;

      /*!\brief Adapts refinement by one level according to the project. Returns true if any cells were refined, false if not.
       * \param mpiGrid grid to refine
       * @return The amount of cells set to refine
       */
      virtual int adaptRefinement( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const;


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
      /*! \brief Returns a list of blocks to loop through when initialising.
       * 
       * The base class version just returns all blocks, which amounts to looping through the whole velocity space.
       * This is very expensive and becomes prohibitive in cases where a large velocity space is needed with only
       * small portions actually containing something. Use with care.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       */
      virtual std::vector<vmesh::GlobalID> findBlocksToInitialize(spatial_cell::SpatialCell* cell,const uint popID) const;
      
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

      /** Integrate the distribution function over the given six-dimensional phase-space cell.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * @param x Starting value of the x-coordinate of the cell.
       * @param y Starting value of the y-coordinate of the cell.
       * @param z Starting value of the z-coordinate of the cell.
       * @param dx The size of the cell in x-direction.
       * @param dy The size of the cell in y-direction.
       * @param dz The size of the cell in z-direction.
       * @param vx Starting value of the vx-coordinate of the cell.
       * @param vy Starting value of the vy-coordinate of the cell.
       * @param vz Starting value of the vz-coordinate of the cell.
       * @param dvx The size of the cell in vx-direction.
       * @param dvy The size of the cell in vy-direction.
       * @param dvz The size of the cell in vz-direction.
       * @param popID Particle species ID.
       * @return The volume average of the distribution function in the given phase space cell.
       * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
       */
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const uint popID) const = 0;
      
      /*!
       Get random number between 0 and 1.0. One should always first initialize the rng.
       */
      Real getRandomNumber() const;
         
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

