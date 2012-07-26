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

//Adjust which blocks are existing based on the current state. 
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

/** Subroutine for setting up a single spatial cell.
 */
bool initSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,creal& zmin,creal& dx,creal& dy,creal& dz,
		     const bool& isRemote);

#endif
