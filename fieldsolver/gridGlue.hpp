#pragma once
#include <fsgrid.hpp>
#include <vector>
#include <array>

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
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                           bool dt2=false);

/*! Copy field solver result (Volume-averaged fields) and store them back into DCCRG
 * \param mpiGrid The DCCRG grid carrying fields.
 * \param cells List of local cells
 * \param volumeFieldsGrid Fieldsolver grid for these quantities
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells);

/*! Copy field derivatives from the appropriate FsGrids and store them back into DCCRG
 *
 * This should only be neccessary for debugging.
 */
void getDerivativesFromFsGrid(FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dperbGrid,
                          FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dmomentsGrid,
                          FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& bgbfieldGrid,
                          dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& cells);

/*! Transfer data into technical grid (boundary info etc.)
 * \param mpiGrid The DCCRG grid carrying rho, rhoV and P
 * \param cells List of local cells
 * \param technicalGrid the target Fieldsolver grid for this information
 *
 * This function assumes that proper grid coupling has been set up.
 */
void setupTechnicalFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, FsGrid< fsgrids::technical, 2>& technicalGrid);

/*! Transfer max timestep data from technical grid back into DCCRG.
 * \param technicalGrid the target Fieldsolver grid for this information
 * \param mpiGrid The DCCRG grid carrying rho, rhoV and P
 * \param cells List of local cells
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getFsGridMaxDt(FsGrid< fsgrids::technical, 2>& technicalGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells);

/*! Transfer background fields into the appropriate FsGrid structure
 *  This requires separate handling, since the source data is not lying
 *  continuous in memory on the DCCRG side.
 *
 * \param mpiGrid The DCCRG grid carrying fieldparam data
 * \param cells List of local cells
 * \param targetGrid Fieldsolver grid for these quantities
 *
 * This function assumes that proper grid coupling has been set up.
 */
void feedBgFieldsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
    const std::vector<CellID>& cells,
    FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid);

/*! Transfer field data from DCCRG cellparams into the appropriate FsGrid structure
 * \param mpiGrid The DCCRG grid carrying fieldparam data
 * \param cells List of local cells
 * \param index Index into the cellparams array from which to copy
 * \param targetGrid Fieldsolver grid for these quantities
 *
 * The cellparams with indices from index to index+numFields are copied over, and
 * have to be continuous in memory.
 *
 * This function assumes that proper grid coupling has been set up.
 */
template< unsigned int numFields > void feedFieldDataIntoFsGrid(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, int index,
      FsGrid< std::array<Real, numFields>, 2>& targetGrid) {

   targetGrid.setupForTransferIn(cells.size());

   for(CellID i : cells) {
      // TODO: This assumes that the field data are lying continuous in memory.
      // Check definition of CellParams in common.h if unsure.
      std::array<Real, numFields>* cellDataPointer = reinterpret_cast<std::array<Real, numFields>*>(
            &(mpiGrid[i]->get_cell_parameters()[index]));
      targetGrid.transferDataIn(i - 1, cellDataPointer);
   }

   targetGrid.finishTransfersIn();
}

/*! Transfer field data from an FsGrid back into the appropriate CellParams slot in DCCRG 
 * \param sourceGrid Fieldsolver grid for these quantities
 * \param mpiGrid The DCCRG grid carrying fieldparam data
 * \param cells List of local cells
 * \param index Index into the cellparams array into which to copy
 *
 * The cellparams with indices from index to index+numFields are copied over, and
 * have to be continuous in memory.
 *
 * This function assumes that proper grid coupling has been set up.
 */
template< unsigned int numFields > void getFieldDataFromFsGrid(
      FsGrid< std::array<Real, numFields>, 2>& sourceGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, int index) {

   sourceGrid.setupForTransferOut(cells.size());

   for(CellID i : cells) {
      // TODO: This assumes that the field data are lying continuous in memory.
      // Check definition of CellParams in common.h if unsure.
      std::array<Real, numFields>* cellDataPointer = reinterpret_cast<std::array<Real, numFields>*>(
            &(mpiGrid[i]->get_cell_parameters()[index]));
      sourceGrid.transferDataOut(i - 1, cellDataPointer);
   }

   sourceGrid.finishTransfersOut();
}

