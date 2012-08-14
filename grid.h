#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include <string>

//Init parallel grid
bool initializeGrid(int argn,
                    char **argc,
                    dccrg::Dccrg<SpatialCell>& mpiGrid,
                    SysBoundary& sysBoundaries,
                    bool& isSysBoundaryDynamic
                   );
//Balance load
void balanceLoad(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);


/*!
  \brief Adjust which blocks are existing based on the current state.

  This function is complete, existence of blocks is adjusted based on
  both spatial and velocity data, and all remote velocity blocks are
  updated. Not thread-safe, uses OpenMP-parallelization.

  \param mpiGrid   The DCCRG grid with spatial cells  
*/
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell>& mpiGrid);


/*!

Updates velocity block lists between remote neighbors and
prepares local copies of remote neighbors to receive velocity block
data. This is needed if one has locally adjusted velocity blocks

\param mpiGrid   The DCCRG grid with spatial cells    
*/
void updateRemoteVelocityBlockLists(dccrg::Dccrg<SpatialCell>& mpiGrid);

//write out system
bool writeGrid(const dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const std::string& name,
               const uint& index,
               const bool& writeRestart);

// Write out diagnostic
bool computeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid, DataReducer& dataReducer, luint tstep);

/*!
  \brief Subroutine for setting up a single spatial cell.

 * @param cell The spatial cell which is to be initialized.
 * @param xmin x-coordinate of the lower left corner of the cell.
 * @param ymin y-coordinate of the lower left corner of the cell.
 * @param zmin z-coordinate of the lower left corner of the cell.
 * @param dx Size of the cell in x-direction.
 * @param dy Size of the cell in y-direction.
 * @param dz Size of the cell in z-direction.
 * @param isRemote If true, the given cell is a remote cell (resides on another process) 
 * and its initial state need not be calculated.
 * @return If true, the cell was initialized successfully. Otherwise an error has 
 * occurred and the simulation should be aborted.
 */
bool initSpatialCell(SpatialCell& cell,
                     SysBoundary& sysBoundaries,
                     creal& xmin, creal& ymin, creal& zmin,
                     creal& dx, creal& dy, creal& dz,
                     const bool& isRemote);

#endif
