#ifndef PROJECT_H
#define PROJECT_H

#include "projects_common.h"
#include <cstdlib>


namespace projects {
   class Project {
      public:
         Project();
         virtual ~Project();
         
         /*! Register parameters that should be read in. */
         static void addParameters();
         /*! Get the value that was read in. */
         virtual void getParameters();
         
          /*! Initialize project. Can be used, e.g., to read in parameters from the input file. */
         virtual bool initialize();

         /*! set background field, should set it for all cells.
          * Currently this function is only called during the initialization.
          * @param cell Pointer to the spatial cell.*/
         virtual void setCellBackgroundField(SpatialCell* cell);
      
         /*!\brief Set the perturbed fields and distribution of a cell according to the default simulation settings.
         * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
         * \param cell Pointer to the cell to set.
         */
         void setCell(SpatialCell* cell);
         
      protected:
         /*! \brief Returns a list of blocks to loop through when initialising.
          * 
          * The base class version just returns all blocks, which amounts to looping through the whole velocity space.
          * This is very expensive and becomes prohibitive in cases where a large velocity space is needed with only
          * small portions actually containing something. Use with care.
          */
         virtual std::vector<uint> findBlocksToInitialize(SpatialCell* cell);
         
         /*! \brief Sets the distribution function in a cell.
          * 
          * Uses the function findBlocksToInitialize and loops through the list returned by it to initialize the cells' velocity space.
          * 
          * \sa findBlocksToInitialize
          */
         void setVelocitySpace(SpatialCell* cell);
         
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
         virtual void calcCellParameters(Real* cellParams,creal& t);
         
         /** Integrate the distribution function over the given six-dimensional phase-space cell.
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
          * @return The volume average of the distribution function in the given phase space cell.
          * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
          */
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz);

      /*!
        Get random number between 0 and 1.0. One should always first initialize the rng.
      */
         Real getRandomNumber();
      /*!  Set random seed (thread-safe). Seed is based on the seed read
        in from cfg + the seedModifier parameter

        \param seedModified d. Seed is based on the seed read in from cfg + the seedModifier parameter                                   
      */
         void setRandomSeed(uint64_t seedModifier);
      /*!
        Set random seed (thread-safe) that is always the same for
        this particular cellID. Can be used to make reproducible
        simulations that do not depend on number of processes or threads.

        \param  cellParams The cell parameters list in each spatial cell
      */
      void setRandomCellSeed(const Real* const cellParams);
      
      private:
         uint seed;
         static char rngStateBuffer[256];
         static random_data rngDataBuffer;
#pragma omp threadprivate(rngStateBuffer,rngDataBuffer) 
   };
   
   Project* createProject();
} // namespace projects


#endif

