#pragma once
#include "../definitions.h"
#include <fsgrid.hpp>
#include <vector>
#include <array>

enum FieldsToCommunicate {
   PERBXVOL,
   PERBYVOL,
   PERBZVOL,
   dPERBXVOLdx,
   dPERBXVOLdy,
   dPERBXVOLdz,
   dPERBYVOLdx,
   dPERBYVOLdy,
   dPERBYVOLdz,
   dPERBZVOLdx,
   dPERBZVOLdy,
   dPERBZVOLdz,
   BGBXVOL,
   BGBYVOL,
   BGBZVOL,
   EXGRADPE,
   EYGRADPE,
   EZGRADPE,
   EXVOL,
   EYVOL,
   EZVOL,
   CURVATUREX,
   CURVATUREY,
   CURVATUREZ,
   N_FIELDSTOCOMMUNICATE
};

std::vector<CellID> mapDccrgIdToFsGridGlobalID(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
					       CellID dccrgID);

/*! Take input moments from DCCRG grid and put them into the Fieldsolver grid
 * \param mpiGrid The DCCRG grid carrying rho, rhoV and P
 * \param cells List of local cells
 * \param momentsGrid Fieldsolver grid for these quantities
 * \param dt2 Whether to copy base moments, or _DT2 moments
 *
 * This function assumes that proper grid coupling has been set up.
 */
void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                           FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                           bool dt2=false);

/*! Copy field solver result (VOLB, VOLE, VOLPERB derivatives, gradpe) and store them back into DCCRG
 * \param mpiGrid The DCCRG grid carrying fields.
 * \param cells List of local cells
 * \param volumeFieldsGrid Fieldsolver grid for these quantities
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volumeFieldsGrid,
			 FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
			 FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
			 FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
			 dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			 const std::vector<CellID>& cells
			 );

/*! Copy background B fields and store them into DCCRG
 * \param mpiGrid The DCCRG grid carrying fields.
 * \param cells List of local cells
 * \param BgBGrid Background field fsgrid
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getBgFieldsAndDerivativesFromFsGrid(
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells
);

/*! Copy field derivatives from the appropriate FsGrids and store them back into DCCRG
 *
 * This should only be neccessary for debugging.
 */
void getDerivativesFromFsGrid(
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dperbGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dmomentsGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells
);


int getNumberOfCellsOnMaxRefLvl(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                const std::vector<CellID>& cells);

void feedBoundaryIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			const std::vector<CellID>& cells,
			FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid);
