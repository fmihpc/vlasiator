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
#define VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID 6 //nearest neighbor in X face direction, f() can propagate to local cells in X dir, and are target for local cells
#define VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID 7 //nearest neighbor in Y face direction, f() can propagate to local cells in Y dir, and are target for local cells
#define VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID 8 //nearest neighbor in Z face direction, f() can propagate to local cells in Z dir, and are target for local cells
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
#define POISSON_NEIGHBORHOOD_ID 20   // Nearest face neighbors 

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
void balanceLoad(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, SysBoundary& sysBoundaries);



/*!

Updates velocity block lists between remote neighbors and
prepares local copies of remote neighbors to receive velocity block
data. This is needed if one has locally adjusted velocity blocks

\param mpiGrid   The DCCRG grid with spatial cells
*/
void updateRemoteVelocityBlockLists(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const int& popID
);

/*! Deallocates all blocks in remote cells in order to save
 *  memory. 
 * \param mpiGrid Spatial grid
 */
void deallocateRemoteCellBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);



/*! Adjust sparse velocity space to make it consistent in all 6 dimensions

 1) Compute which blocks have content (done for all cells in mpiGrid)
 2) Adjust local velocity blocks. That is, make sure blocks exist which have content, or have
 neighbors with content in all 6-dimensions. This is done for cells in cellsToAdjust list.
 3) Make sure remote cells are up-to-date and ready to receive data, if doPrepareToReceiveBlocks is true.

 Note that block existence does not use vlasov stencil as it is important to also include diagonals to avoid massloss

 \param mpiGrid  Parallel grid with spatial cells
 \param cellsToAdjust  List of cells that are adjusted, that is cells which blocks are added or removed. 
 \param doPrepareToReceiveBlocks If true, then remote cells are set up so that velocity space data can be received. Global operation, value has to be the same for all processes.
 
*/
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<uint64_t>& cellsToAdjust,
                          bool doPrepareToReceiveBlocks,
                            const int& popID);



/*! Estimates memory consumption and writes it into logfile. Collective operation on MPI_COMM_WORLD
 * \param mpiGrid Spatial grid
 */
void report_grid_memory_consumption(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

/*! Shrink to fit velocity space data to save memory.
 * \param mpiGrid Spatial grid
 */
void shrink_to_fit_grid_data(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

/** Validate the velocity mesh structure. This function is only relevant for 
 * the AMR mesh. It makes sure that the mesh structure is valid for all spatial cells, 
 * i.e., that each velocity block has at most one refinement level difference to 
 * its neighbors (in spatial and velocity space).
 * @param mpiGrid Parallel grid.
 * @return If true, the mesh is valid. Otherwise an error has occurred and the simulation 
 * should be aborted.*/
bool validateMesh(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const int& popID);

#endif
