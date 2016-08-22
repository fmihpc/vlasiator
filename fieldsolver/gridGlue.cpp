#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "../grid.h"
#include "../spatial_cell.hpp"
#include "../definitions.h"
#include "../common.h"
#include "gridGlue.hpp"

void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid, bool dt2 /*=false*/) {

   momentsGrid.setupForTransferIn(cells.size());

   int startindex=CellParams::RHO;
   if(dt2) {
      startindex=CellParams::RHO_DT2;
   }

   for(CellID i : cells) {
      // TODO: This assumes that RHO, RHOV and P (diagonals) are lying continuous in memory.
      // Check definition of CellParams in common.h if unsure.
      std::array<Real, fsgrids::moments::N_MOMENTS>* cellDataPointer = reinterpret_cast<std::array<Real, fsgrids::moments::N_MOMENTS>*>(
            &(mpiGrid[i]->get_cell_parameters()[startindex]));
      momentsGrid.transferDataIn(i - 1, cellDataPointer);
   }

   momentsGrid.finishTransfersIn();
}


void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells) {

   volumeFieldsGrid.setupForTransferOut(cells.size());
   // Fill the transfer buffers from the spatial cell structs
   std::vector< std::array<Real, fsgrids::volfields::N_VOL> > transferBuffer(cells.size());

   size_t count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[count++];

      // Data needs to be collected from some different places for this grid.
      auto cellParams = mpiGrid[i]->get_cell_parameters();
      thisCellData->at(fsgrids::volfields::PERBXVOL) = cellParams[CellParams::PERBXVOL];
      thisCellData->at(fsgrids::volfields::PERBYVOL) = cellParams[CellParams::PERBYVOL];
      thisCellData->at(fsgrids::volfields::PERBZVOL) = cellParams[CellParams::PERBZVOL];
      thisCellData->at(fsgrids::volfields::EXVOL) = cellParams[CellParams::EXVOL];
      thisCellData->at(fsgrids::volfields::EYVOL) = cellParams[CellParams::EYVOL];
      thisCellData->at(fsgrids::volfields::EZVOL) = cellParams[CellParams::EZVOL];
      thisCellData->at(fsgrids::volfields::dPERBXVOLdy) = mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBXVOLdy];
      thisCellData->at(fsgrids::volfields::dPERBXVOLdz) = mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBXVOLdz];
      thisCellData->at(fsgrids::volfields::dPERBYVOLdx) = mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBYVOLdx];
      thisCellData->at(fsgrids::volfields::dPERBYVOLdz) = mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBYVOLdz];
      thisCellData->at(fsgrids::volfields::dPERBZVOLdx) = mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBZVOLdx];
      thisCellData->at(fsgrids::volfields::dPERBZVOLdy) = mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBZVOLdy];
      volumeFieldsGrid.transferDataOut(i - 1, thisCellData);
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
      technicalGrid.transferDataIn(i - 1,thisCellData);
   }

   technicalGrid.finishTransfersIn();
}

void getFsGridMaxDt(FsGrid< fsgrids::technical, 2>& technicalGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells) {

   technicalGrid.setupForTransferOut(cells.size());

   // Buffer to store contents of the grid
   std::vector<fsgrids::technical> transferBuffer(cells.size());
   size_t count=0;

   for(CellID i : cells) {
      fsgrids::technical* thisCellData = &transferBuffer[count++];
      technicalGrid.transferDataOut(i - 1, thisCellData);
   }

   technicalGrid.finishTransfersOut();

   // After the transfer is completed, stuff the recieved maxFDt into the cells.
   count=0;
   for(CellID i : cells) {
      mpiGrid[i]->get_cell_parameters()[CellParams::MAXFDT] = transferBuffer[count++].maxFsDt;
   }
}
