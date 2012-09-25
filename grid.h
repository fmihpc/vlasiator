#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include <string>

/*!
  \brief Initialize parallel grid
*/
void initializeGrid(int argn,
                    char **argc,
                    dccrg::Dccrg<SpatialCell>& mpiGrid,
                    SysBoundary& sysBoundaries);

/*!
  \brief Balance load

    \param[in,out] mpiGrid The DCCRG grid with spatial cells
*/
  void balanceLoad(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);


/*!
  \brief Adjust which blocks are existing based on the current state.

  Before calling this function it is assumed that remote cells have
  the correct velocity block lists. If one has done local changes to
  these lists, then they should be updated using
  updateRemoteVelocityBlockLists() before calling this function.

  This function will update the global view of the block-structure so that it is consistent and up-to-date on all processes:  
  - Computes which blocks have contents according to threshold
  - Adds/keeps blocks if they have content, or if their neighbors in spatial or velocity space have contents
  - Removes blocks which do not fulfill above mentioned criteria
  - Updates velocity block lists of remote cells using updateRemoteVelocityBlockLists(). Note that data in blocks is not transferred!
  - Updates movers so that their block-dependent internal datastructures are up-to-date, if reInitMover is true.

  
  \param[in,out] mpiGrid The DCCRG grid with spatial cells
  \param[in] reInitMover Does it also re-initialize the solvers based on the new
  block structure? The default is true, and that should always be used
  unless it is called before solvers have been initialized.
  
  
*/
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell>& mpiGrid, bool reInitMover=true);


/*!

Updates velocity block lists between remote neighbors and
prepares local copies of remote neighbors to receive velocity block
data. This is needed if one has locally adjusted velocity blocks

\param mpiGrid   The DCCRG grid with spatial cells    
*/
void updateRemoteVelocityBlockLists(dccrg::Dccrg<SpatialCell>& mpiGrid);





void setNotSysBoundaryCell(SpatialCell* cell);

#endif
