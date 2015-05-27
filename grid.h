#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "sysboundary/sysboundary.h"
#include "projects/project.h"
#include <string>

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
                          const std::vector<CellID>& cellsToAdjust,
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
