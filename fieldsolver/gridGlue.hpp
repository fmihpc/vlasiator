#pragma once
#include <fsgrid.hpp>
#include <vector>
#include <array>

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

int getNumberOfCellsOnMaxRefLvl(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                const std::vector<CellID>& cells);

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
      const std::vector<CellID>& cells, int cellParamsIndex, 
      FsGrid< std::array<Real, numFields>, 2>& targetGrid) {
   
   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   targetGrid.setupForTransferIn(nCells);

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   for(CellID dccrgId : cells) {
     const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
     // TODO: This assumes that the field data are lying continuous in memory.
     // Check definition of CellParams in common.h if unsure.
     std::array<Real, numFields>* cellDataPointer = reinterpret_cast<std::array<Real, numFields>*>(
                                                     &(mpiGrid[dccrgId]->get_cell_parameters()[cellParamsIndex]));
     for (auto fsgridId : fsgridIds) {
       targetGrid.transferDataIn(fsgridId, cellDataPointer);
     }
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
   
   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   std::vector< std::array<Real, numFields> > transferBuffer(nCells);
   std::vector< std::array<Real, numFields>*> transferBufferPointer;
   sourceGrid.setupForTransferOut(nCells);

   int k = 0;
   for(CellID dccrgId : cells) {
      // TODO: This assumes that the field data are lying continuous in memory.
      // Check definition of CellParams in common.h if unsure.
      
      //std::array<Real, numFields>* cellDataPointer = reinterpret_cast<std::array<Real, numFields>*>(
      //      &(mpiGrid[dccrgId]->get_cell_parameters()[index]));

      transferBufferPointer.push_back(&transferBuffer[k]);
      
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
      for (auto fsgridId : fsgridIds) {
         std::array<Real, numFields>* cellDataPointer = &transferBuffer[k++];
         sourceGrid.transferDataOut(fsgridId, cellDataPointer);
      }
   }

   sourceGrid.finishTransfersOut();

   // Average data in transferBuffer
#pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {
      
      CellID dccrgId = cells[i];

      // Set cell data to 0
      for (int iField = 0; iField < numFields; ++iField) {
         mpiGrid[dccrgId]->get_cell_parameters()[index+iField] = 0.0;
      }      
      
      // Calculate the number of fsgrid cells we need to average into the current dccrg cell
      auto refLvl = mpiGrid.mapping.get_refinement_level(dccrgId);
      int nCells = pow(pow(2,mpiGrid.mapping.get_maximum_refinement_level() - mpiGrid.mapping.get_refinement_level(dccrgId)),3);

      for(int iCell = 0; iCell < nCells; ++iCell) {

         std::array<Real, numFields>* cellDataPointer = transferBufferPointer[i] + iCell;

         for (int iField = 0; iField < numFields; ++iField) {
            mpiGrid[dccrgId]->get_cell_parameters()[index+iField] += cellDataPointer->at(iField);
         }
      }

      for (int iField = 0; iField < numFields; ++iField) {
         mpiGrid[dccrgId]->get_cell_parameters()[index+iField] /= nCells;
      }      
   }

}
