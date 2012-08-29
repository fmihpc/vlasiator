#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include <string>

//Init parallel grid
void initializeGrid(int argn,
                    char **argc,
                    dccrg::Dccrg<SpatialCell>& mpiGrid,
                    SysBoundary& sysBoundaries);

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
bool computeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                       DataReducer& dataReducer);

void setNotSysBoundaryCell(SpatialCell* cell);

#endif
