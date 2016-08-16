#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "../grid.h"
#include "../spatial_cell.hpp"
#include "../definitions.h"
#include "../common.h"
#include "gridGlue.hpp"

void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid) {

   momentsGrid.setupForTransferIn(cells.size());


   for(CellID i : cells) {
      // TODO: This assumes that RHO, RHOV and P (diagonals) are lying continuous in memory.
      // Check definition of CellParams in common.h if unsure.
      std::array<Real, fsgrids::moments::N_MOMENTS>* cellDataPointer = reinterpret_cast<std::array<Real, fsgrids::moments::N_MOMENTS>*>(
            &(mpiGrid[i]->get_cell_parameters()[CellParams::RHO]));
      momentsGrid.transferDataIn(i - 1, cellDataPointer);
   }

   momentsGrid.finishTransfersIn();
}


void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells) {

   volumeFieldsGrid.setupForTransferOut(cells.size());

   for(CellID i : cells) {
      std::array<Real, fsgrids::volfields::N_VOL>* cellDataPointer = reinterpret_cast<std::array<Real, fsgrids::volfields::N_VOL>*>(
            &(mpiGrid[i]->get_cell_parameters()[CellParams::RHO]));
      volumeFieldsGrid.transferDataOut(i - 1, cellDataPointer);
   }

   volumeFieldsGrid.finishTransfersOut();
}


void setupTechnicalFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, FsGrid< fsgrids::technical, 2>& technicalGrid) {

   technicalGrid.setupForTransferIn(cells.size());

   // Fill the transfer buffers from the spatial cell structs
   std::vector<fsgrids::technical> transferBuffer(cells.size());

   size_t count=0;
   for(CellID i : cells) {

      fsgrids::technical* thisCellData = &transferBuffer[count++];
      // Data needs to be collected from some different places for this grid.
      thisCellData->sysBoundaryFlag = mpiGrid[i]->sysBoundaryFlag;
      thisCellData->sysBoundaryLayer = mpiGrid[i]->sysBoundaryLayer;
      thisCellData->maxFsDt = mpiGrid[i]->get_cell_parameters()[CellParams::MAXFDT];
      technicalGrid.transferDataIn(i,thisCellData);
   }

   technicalGrid.finishTransfersIn();
}
