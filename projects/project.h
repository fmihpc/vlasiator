/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2011-2015 Finnish Meteorological Institute
 * 
 */

#ifndef PROJECT_H
#define PROJECT_H

#include "../spatial_cell.hpp"

namespace projects {
   class Project {
    public:
      Project();
      virtual ~Project();
      
      /*! Register parameters that should be read in. */
      static void addParameters();
      
      virtual Real getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const int& popID) const;
      
      /*! Get the value that was read in. */
      virtual void getParameters();
      
      /*! Initialize project. Can be used, e.g., to read in parameters from the input file. */
      virtual bool initialize();
      
      bool initialized();
      
      /*! set background field, should set it for all cells.
       * Currently this function is only called during the initialization.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * @param cell Pointer to the spatial cell.*/
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
      
      /*!\brief Set the perturbed fields and distribution of a cell according to the default simulation settings.
       * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * \param cell Pointer to the cell to set.
       */
      void setCell(spatial_cell::SpatialCell* cell);
         
      Real setVelocityBlock(spatial_cell::SpatialCell* cell,const vmesh::LocalID& blockLID,const int& popID) const;

    protected:
      /*! \brief Returns a list of blocks to loop through when initialising.
       * 
       * The base class version just returns all blocks, which amounts to looping through the whole velocity space.
       * This is very expensive and becomes prohibitive in cases where a large velocity space is needed with only
       * small portions actually containing something. Use with care.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       */
      virtual std::vector<vmesh::GlobalID> findBlocksToInitialize(spatial_cell::SpatialCell* cell,const int& popID) const;
      
      /*! \brief Sets the distribution function in a cell.
       * 
       * Uses the function findBlocksToInitialize and loops through the list returned by it to initialize the cells' velocity space.
       * NOTE: This function is called inside parallel region so it must be declared as const.
       * 
       * \sa findBlocksToInitialize
       */
      void setVelocitySpace(const int& popID,spatial_cell::SpatialCell* cell) const;
         
      /** Calculate parameters for the given spatial cell at the given time.
       * Here you need to set values for the following array indices:
       * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
       * 
       * Currently this function is only called during initialization.
       * 
       * The following array indices contain the coordinates of the "lower left corner" of the cell: 
       * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
       * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
       * @param cellParams Array containing cell parameters.
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
                                         const int& popID) const;
      
      /*!
       Get random number between 0 and 1.0. One should always first initialize the rng.
       */
      Real getRandomNumber(spatial_cell::SpatialCell* cell) const;
         
      void printPopulations();
      
      virtual bool rescalesDensity(const int& popID) const;
      void rescaleDensity(spatial_cell::SpatialCell* cell,const int& popID) const;
      
      /*!  Set random seed (thread-safe). Seed is based on the seed read
       in from cfg + the seedModifier parameter
       * 
       \param seedModified d. Seed is based on the seed read in from cfg + the seedModifier parameter                                   
       */
      void setRandomSeed(spatial_cell::SpatialCell* cell,uint64_t seedModifier) const;
      /*!
       Set random seed (thread-safe) that is always the same for
       this particular cellID. Can be used to make reproducible
       simulations that do not depend on number of processes or threads.
       * 
       \param  cellParams The cell parameters list in each spatial cell
       */
      void setRandomCellSeed(spatial_cell::SpatialCell* cell,const Real* const cellParams) const;

    private:
      uint seed;
      static char rngStateBuffer[256];
      static random_data rngDataBuffer;
      #pragma omp threadprivate(rngStateBuffer,rngDataBuffer)

      bool baseClassInitialized;                      /**< If true, base class has been initialized.*/
      std::vector<std::string> popNames;              /**< Name(s) of particle population(s), read from configuration file.*/
      std::vector<int> popCharges;                    /**< Particle population charge(s), read from configuration file.*/
      std::vector<std::string> popMassUnits;          /**< Units in which particle population mass(es) were given,
                                                       * read from configuration file.*/
      std::vector<double> popMasses;                  /**< Particle population mass(es) in chosen units.
                                                       * Read from configuration file.*/
      std::vector<std::string> popMeshNames;          /**< Name of the velocity mesh each species should use.*/
      std::vector<double> popSparseMinValue;          /**< Sparse mesh threshold value for the population.*/
   };
   
   Project* createProject();
} // namespace projects


#endif

