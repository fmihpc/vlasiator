#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include "projects/project.h"
#include <string>


//neighborhoods, these are initialized in initializeGrid

#define FIELD_SOLVER_NEIGHBORHOOD_ID 1
#define VLASOV_SOLVER_NEIGHBORHOOD_ID 2   //up to third(PPM) neighbor in each face direction
#define VLASOV_SOLVER_X_NEIGHBORHOOD_ID 3 //up to third(PPM) neighbor in x face directions
#define VLASOV_SOLVER_Y_NEIGHBORHOOD_ID 4 //up to third(PPM) neighbor in y face directions
#define VLASOV_SOLVER_Z_NEIGHBORHOOD_ID 5 //up to third(PPM) neighbor in z face directions
#define VLASOV_SOLVER_SOURCE_X_NEIGHBORHOOD_ID 6 //nearest neighbor in X face direction, f() can propagate to local cells in X dir
#define VLASOV_SOLVER_SOURCE_Y_NEIGHBORHOOD_ID 7 //nearest neighbor in Y face direction, f() can propagate to local cells in Y dir
#define VLASOV_SOLVER_SOURCE_Z_NEIGHBORHOOD_ID 8 //nearest neighbor in Z face direction, f() can propagate to local cells in Z dir
#define SYSBOUNDARIES_NEIGHBORHOOD_ID 9 // When classifying sysboundaries, all 26 nearest neighbors are included,
#define SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID 10 //Up to second nearest neighbors in all directions (also diagonals)
#define NEAREST_NEIGHBORHOOD_ID 11  //nearest neighbors
#define FULL_NEIGHBORHOOD_ID 12      //Up to second nearest neighbors in all directions (also diagonals) + vlasov solver neighborhood
#define DIST_FUNC_NEIGHBORHOOD_ID 13 //nearest neighbors in all directions (also diagonals) + vlasov solver neighborhood
#define SHIFT_P_X_NEIGHBORHOOD_ID 14 //Shift in +x direction
#define SHIFT_P_Y_NEIGHBORHOOD_ID 15 //Shift in +y direction
#define SHIFT_P_Z_NEIGHBORHOOD_ID 16 //Shift in +z direction
#define SHIFT_M_X_NEIGHBORHOOD_ID 17 //Shift in -x direction
#define SHIFT_M_Y_NEIGHBORHOOD_ID 18 //Shift in -y direction
#define SHIFT_M_Z_NEIGHBORHOOD_ID 19 //Shift in -z direction

//fieldsolver stencil.
#define FS_STENCIL_WIDTH 2
//Vlasov propagator stencils in ordinary space, velocity space may be
//higher. Assume H4 (or H5) for PPM, H6 for PQM
#ifdef TRANS_SEMILAG_PLM
#define  VLASOV_STENCIL_WIDTH 1
#endif
#ifdef TRANS_SEMILAG_PPM
#define  VLASOV_STENCIL_WIDTH 2
#endif
#ifdef TRANS_SEMILAG_PQM
#define  VLASOV_STENCIL_WIDTH 3
#endif

/*!
  \brief Initialize parallel grid
*/
void initializeGrid(
   int argn,
   char **argc,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   Project& project
);

/*!
  \brief Balance load

    \param[in,out] mpiGrid The DCCRG grid with spatial cells
*/
void balanceLoad(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);


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
  
*/
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);


/*!

Updates velocity block lists between remote neighbors and
prepares local copies of remote neighbors to receive velocity block
data. This is needed if one has locally adjusted velocity blocks

\param mpiGrid   The DCCRG grid with spatial cells
*/
void updateRemoteVelocityBlockLists(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
);

/*! Deallocates all blocks in remote cells in order to save
 *  memory. 
 * \param mpiGrid Spatial grid
 */
void deallocateRemoteCellBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);




//subroutine to adjust blocks of local cells; remove/add based on user-defined limits
bool adjust_local_velocity_blocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);




/*! Estimates memory consumption and writes it into logfile. Collective operation on MPI_COMM_WORLD
 * \param mpiGrid Spatial grid
 */
void report_grid_memory_consumption(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

/*! Shrink to fit velocity space data to save memory.
 * \param mpiGrid Spatial grid
 */
void shrink_to_fit_grid_data(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);



#endif
