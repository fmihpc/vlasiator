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
                           FsGrid< fsgrids::technical, 2>& technicalGrid,
                           bool dt2=false);

/*! Copy field solver result (VOLB, VOLE, VOLPERB derivatives, gradpe) and store them back into DCCRG
 * \param mpiGrid The DCCRG grid carrying fields.
 * \param cells List of local cells
 * \param volumeFieldsGrid Fieldsolver grid for these quantities
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
			 FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
			 FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
			 FsGrid< fsgrids::technical, 2>& technicalGrid,
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
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
   FsGrid< fsgrids::technical, 2>& technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells
);

/*! Copy field derivatives from the appropriate FsGrids and store them back into DCCRG
 *
 * This should only be neccessary for debugging.
 */
void getDerivativesFromFsGrid(
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dperbGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dmomentsGrid,
   FsGrid< fsgrids::technical, 2>& technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells
);


int getNumberOfCellsOnMaxRefLvl(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                const std::vector<CellID>& cells);


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
   FsGrid< fsgrids::technical, 2>& technicalGrid,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   int index
) {
   
   cint nCellsOnMaxRefLvl = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   std::vector< std::array<Real, numFields> > transferBufferData(nCellsOnMaxRefLvl);
   std::vector< std::array<Real, numFields>*> transferBufferPointerData;
   std::vector< fsgrids::technical > transferBufferTechnical(nCellsOnMaxRefLvl);
   std::vector< fsgrids::technical*> transferBufferPointerTechnical;
   sourceGrid.setupForTransferOut(nCellsOnMaxRefLvl);
   technicalGrid.setupForTransferOut(nCellsOnMaxRefLvl);

   int k = 0;
   for(CellID dccrgId : cells) {
      // TODO: This assumes that the field data are lying continuous in memory.
      // Check definition of CellParams in common.h if unsure.
      
      //std::array<Real, numFields>* cellDataPointer = reinterpret_cast<std::array<Real, numFields>*>(
      //      &(mpiGrid[dccrgId]->get_cell_parameters()[index]));

      transferBufferPointerData.push_back(&transferBufferData[k]);
      transferBufferPointerTechnical.push_back(&transferBufferTechnical[k]);
      
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
      for (auto fsgridId : fsgridIds) {
         std::array<Real, numFields>* cellDataPointer = &transferBufferData[k];
         sourceGrid.transferDataOut(fsgridId, cellDataPointer);
         fsgrids::technical* thisCellDataTechnical = &transferBufferTechnical[k];
         technicalGrid.transferDataOut(fsgridId, thisCellDataTechnical);
         k++;
      }
   }

   sourceGrid.finishTransfersOut();
   technicalGrid.finishTransfersOut();

   // Average data in transferBuffer
   // Disregard DO_NOT_COMPUTE cells
#pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {
      
      const CellID dccrgId = cells[i];
      
      // Set cell data to 0
      for (int iField = 0; iField < numFields; ++iField) {
         mpiGrid[dccrgId]->get_cell_parameters()[index+iField] = 0.0;
      }      
      
      // Calculate the number of fsgrid cells we loop through
      cint nCells = pow(pow(2,mpiGrid.mapping.get_maximum_refinement_level() - mpiGrid.mapping.get_refinement_level(dccrgId)),3);
      // Count the number of fsgrid cells we need to average into the current dccrg cell
      int nCellsToSum = 0;
      
      for(int iCell = 0; iCell < nCells; ++iCell) {
         if ((transferBufferPointerTechnical[i] + iCell)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         } else {
            nCellsToSum++;
            std::array<Real, numFields>* cellDataPointer = transferBufferPointerData[i] + iCell;
            
            for (int iField = 0; iField < numFields; ++iField) {
               mpiGrid[dccrgId]->get_cell_parameters()[index+iField] += cellDataPointer->at(iField);
            }
         }
      }
      
      if (nCellsToSum > 0) {
         for (int iField = 0; iField < numFields; ++iField) {
            mpiGrid[dccrgId]->get_cell_parameters()[index+iField] /= nCellsToSum;
         }
      }
   }

}
