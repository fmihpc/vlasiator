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
			FsGrid< fsgrids::technical, 2> & technicalGrid);


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
   FsGrid< std::array<Real, numFields>, FS_STENCIL_WIDTH> & sourceGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
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

/*Compute coupling DCCRG <=> FSGRID 

  onDccrgMapRemoteProcess   maps fsgrid processes (key) => set of dccrg cellIDs owned by current rank that map to  the fsgrid cells owned by fsgrid process (val)

  onFsgridMapRemoteProcess  maps dccrg processes  (key) => set of dccrg cellIDs owned by dccrg-process that map to current rank fsgrid cells 
  onFsgridMapCells          maps remote dccrg CellIDs to local fsgrid cells
*/

template <typename T, int stencil> void computeCoupling(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells,
      FsGrid< T, stencil>& momentsGrid,
      FsGridCouplingInformation* coupling) {

  //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
  std::vector<CellID> dccrgCells = cells;
  std::sort(dccrgCells.begin(), dccrgCells.end());

  //make sure the datastructures are clean
  coupling->onDccrgMapRemoteProcess.clear();
  coupling->onFsgridMapRemoteProcess.clear();
  coupling->onFsgridMapCells.clear();

  //size of fsgrid local part
  const std::array<int, 3> gridDims(momentsGrid.getLocalSize());

  //Compute what we will receive, and where it should be stored
  for (int k=0; k<gridDims[2]; k++) {
    for (int j=0; j<gridDims[1]; j++) {
      for (int i=0; i<gridDims[0]; i++) {
        const std::array<int, 3> globalIndices = momentsGrid.getGlobalIndices(i,j,k);
        const dccrg::Types<3>::indices_t  indices = {{(uint64_t)globalIndices[0],
                        (uint64_t)globalIndices[1],
                        (uint64_t)globalIndices[2]}}; //cast to avoid warnings
        CellID dccrgCell = mpiGrid.get_existing_cell(indices, 0, mpiGrid.mapping.get_maximum_refinement_level());

        int process = mpiGrid.get_process(dccrgCell);
        int64_t  fsgridLid = momentsGrid.LocalIDForCoords(i,j,k);
        //int64_t  fsgridGid = momentsGrid.GlobalIDForCoords(i,j,k);
        coupling->onFsgridMapRemoteProcess[process].insert(dccrgCell); //cells are ordered (sorted) in set
        coupling->onFsgridMapCells[dccrgCell].push_back(fsgridLid);
      }
    }
  }

  // Compute where to send data and what to send
  for(uint64_t i=0; i< dccrgCells.size(); i++) {
     //compute to which processes this cell maps
     std::vector<CellID> fsCells = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgCells[i]);

     //loop over fsgrid cells which this dccrg cell maps to
     for (auto const &fsCellID : fsCells) {
       int process = momentsGrid.getTaskForGlobalID(fsCellID).first; //process on fsgrid
       coupling->onDccrgMapRemoteProcess[process].insert(dccrgCells[i]); //add to map
     }
  }
}

