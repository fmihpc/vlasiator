#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include "datareducer.h"
#include <string>

//Init parallel grid
bool initializeGrid(int argn, char **argc,dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);
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

/*!

\brief Write out system into a vlsv file

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param name       File name prefix, file will be called "name.index.vlsv"
\param index      File index, file will be called "name.index.vlsv"
\param writeRestart If true, the full velocity distribution will be written out.
*/


bool writeGrid(const dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const std::string& name,
               const uint& index,
               const bool& writeRestart);


/*!

\brief Read in state from a vlsv file in order to restart simulations
*/


bool readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
              const std::string& name);



/*!

\brief Write out simulation diagnostics into diagnostic.txt

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute diagnostic data
\param tstep Current simulation step, first colmun in file
*/

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
bool initSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,creal& zmin,creal& dx,creal& dy,creal& dz,
		     const bool& isRemote);

#endif
