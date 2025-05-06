/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef GRID_H
#define GRID_H

#include "definitions.h"
#include "fsgrid.hpp"
#include "projects/project.h"
#include "spatial_cells/block_adjust_wrapper.hpp"
#include "spatial_cells/spatial_cell_wrapper.hpp"
#include "sysboundary/sysboundary.h"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <span>
#include <string>

struct FieldSolverData {
   FieldSolverGrid& fsgrid;

   std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perB;
   std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perBDt2;
   std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> E;
   std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> EDt2;
   std::span<const std::array<Real, fsgrids::ehall::N_EHALL>> EHall;
   std::span<const std::array<Real, fsgrids::egradpe::N_EGRADPE>> EGradPe;
   std::span<const std::array<Real, fsgrids::egradpe::N_EGRADPE>> EGradPeDt2;
   std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments;
   std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> momentsDt2;
   std::span<const std::array<Real, fsgrids::dperb::N_DPERB>> dPerB;
   std::span<const std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dMoments;
   std::span<const std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dMomentsDt2;
   std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> BgB;
   std::span<const std::array<Real, fsgrids::volfields::N_VOL>> vol;
   std::span<const fsgrids::technical> technical;

   FieldSolverData(const fsgrid::FsData<std::array<Real, fsgrids::bfield::N_BFIELD>>& perb,
           const fsgrid::FsData<std::array<Real, fsgrids::bfield::N_BFIELD>>& perbdt2,
           const fsgrid::FsData<std::array<Real, fsgrids::efield::N_EFIELD>>& e,
           const fsgrid::FsData<std::array<Real, fsgrids::efield::N_EFIELD>>& edt2,
           const fsgrid::FsData<std::array<Real, fsgrids::ehall::N_EHALL>>& ehall,
           const fsgrid::FsData<std::array<Real, fsgrids::egradpe::N_EGRADPE>>& egradpe,
           const fsgrid::FsData<std::array<Real, fsgrids::egradpe::N_EGRADPE>>& egradpedt2,
           const fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& moments,
           const fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& momentsdt2,
           const fsgrid::FsData<std::array<Real, fsgrids::dperb::N_DPERB>>& dperb,
           const fsgrid::FsData<std::array<Real, fsgrids::dmoments::N_DMOMENTS>>& dmoments,
           const fsgrid::FsData<std::array<Real, fsgrids::dmoments::N_DMOMENTS>>& dmomentsdt2,
           const fsgrid::FsData<std::array<Real, fsgrids::bgbfield::N_BGB>>& bgb,
           const fsgrid::FsData<std::array<Real, fsgrids::volfields::N_VOL>>& vol,
           const fsgrid::FsData<fsgrids::technical>& technical, FieldSolverGrid& fsgrid)
       : fsgrid(fsgrid), perB(perb.view()), perBDt2(perbdt2.view()), E(e.view()), EDt2(edt2.view()),
         EHall(ehall.view()), EGradPe(egradpe.view()), EGradPeDt2(egradpedt2.view()), moments(moments.view()),
         momentsDt2(momentsdt2.view()), dPerB(dperb.view()), dMoments(dmoments.view()), dMomentsDt2(dmomentsdt2.view()),
         BgB(bgb.view()), vol(vol.view()), technical(technical.view()) {}
};

/*!
  \brief Initialize DCCRG and fsgrids
*/
void initializeGrids(int argn, char** argc, dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                     fsgrid::FsData<std::array<Real, fsgrids::bfield::N_BFIELD>>& perb,
                     fsgrid::FsData<std::array<Real, fsgrids::bgbfield::N_BGB>>& bgb,
                     fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& moments,
                     fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& momentsdt2,
                     fsgrid::FsData<std::array<Real, fsgrids::dmoments::N_DMOMENTS>>& dmoments,
                     fsgrid::FsData<std::array<Real, fsgrids::efield::N_EFIELD>>& e,
                     fsgrid::FsData<std::array<Real, fsgrids::egradpe::N_EGRADPE>>& egradpe,
                     fsgrid::FsData<std::array<Real, fsgrids::volfields::N_VOL>>& vol,
                     fsgrid::FsData<fsgrids::technical>& technical, FieldSolverGrid& fsgrid,
                     SysBoundary& sysBoundaries, Project& project);

/*!
  \brief Balance load

    \param[in,out] mpiGrid The DCCRG grid with spatial cells
*/
void balanceLoad(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, SysBoundary& sysBoundaries, std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
   bool doTranslationLists = true);

/* helper for calculating AMR flags and cell lists and building pencils
 */
void prepareAMRLists(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);


/*!

Updates velocity block lists between remote neighbors and
prepares local copies of remote neighbors to receive velocity block
data. This is needed if one has locally adjusted velocity blocks

\param mpiGrid   The DCCRG grid with spatial cells
*/
void updateRemoteVelocityBlockLists(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint popID,
   const uint neighborhood=Neighborhoods::DIST_FUNC
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
                          const uint popID);

/*! Shrink to fit velocity space data to save memory.
 * \param mpiGrid Spatial grid
 */
void shrink_to_fit_grid_data(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

void setFaceNeighborRanks( dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid );

/*! Map grid refinement to FsGrid
 */
void mapRefinement(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid);

/*! Refine spatial cells and update necessary information
 * \param mpiGrid Spatial grid
 * \param fsgrid Technical grid
 * \param sysBoundaries System boundaries
 * \param project Project used
 * \param useStatic Used for forcing static refinement on restart. Negative values use adaptive refinement, non-negative values correspond to static refinement pass in Project::forceRefinement
 */
bool adaptRefinement(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid, SysBoundary& sysBoundaries, Project& project, int useStatic = -1);

void recalculateLocalCellsCache(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

#endif
