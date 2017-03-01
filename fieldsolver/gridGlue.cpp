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

   // Setup transfer buffers
   std::vector< std::array<Real, fsgrids::moments::N_MOMENTS> > transferBuffer(cells.size());

   // Fill from cellParams
   size_t count=0;
   for(CellID i : cells) {
      auto cellParams = mpiGrid[i]->get_cell_parameters();
      std::array<Real, fsgrids::moments::N_MOMENTS>* thisCellData = &transferBuffer[count++];
    
      if(!dt2) {
        thisCellData->at(fsgrids::moments::RHO) = cellParams[CellParams::RHO];
        thisCellData->at(fsgrids::moments::RHOVX) = cellParams[CellParams::RHOVX];
        thisCellData->at(fsgrids::moments::RHOVY) = cellParams[CellParams::RHOVY];
        thisCellData->at(fsgrids::moments::RHOVZ) = cellParams[CellParams::RHOVZ];
        thisCellData->at(fsgrids::moments::P_11) = cellParams[CellParams::P_11];
        thisCellData->at(fsgrids::moments::P_22) = cellParams[CellParams::P_22];
        thisCellData->at(fsgrids::moments::P_33) = cellParams[CellParams::P_33];
      } else {
        thisCellData->at(fsgrids::moments::RHO) = cellParams[CellParams::RHO_DT2];
        thisCellData->at(fsgrids::moments::RHOVX) = cellParams[CellParams::RHOVX_DT2];
        thisCellData->at(fsgrids::moments::RHOVY) = cellParams[CellParams::RHOVY_DT2];
        thisCellData->at(fsgrids::moments::RHOVZ) = cellParams[CellParams::RHOVZ_DT2];
        thisCellData->at(fsgrids::moments::P_11) = cellParams[CellParams::P_11_DT2];
        thisCellData->at(fsgrids::moments::P_22) = cellParams[CellParams::P_22_DT2];
        thisCellData->at(fsgrids::moments::P_33) = cellParams[CellParams::P_33_DT2];
      }

      momentsGrid.transferDataIn(i - 1, thisCellData);
   }

   // Finish the actual transfer
   momentsGrid.finishTransfersIn();
}


void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells) {

   // Setup transfer buffers
   std::vector< std::array<Real, fsgrids::volfields::N_VOL> > transferBuffer(cells.size());

   // Setup transfer pointers
   volumeFieldsGrid.setupForTransferOut(cells.size());
   size_t count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[count++];

      volumeFieldsGrid.transferDataOut(i - 1, thisCellData);
   }
   // Do the transfer
   volumeFieldsGrid.finishTransfersOut();

   // Distribute data from the transfer buffer back into the appropriate mpiGrid places
   count = 0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[count++];
      auto cellParams = mpiGrid[i]->get_cell_parameters();

      cellParams[CellParams::PERBXVOL]                          = thisCellData->at(fsgrids::volfields::PERBXVOL);
      cellParams[CellParams::PERBYVOL]                          = thisCellData->at(fsgrids::volfields::PERBYVOL);
      cellParams[CellParams::PERBZVOL]                          = thisCellData->at(fsgrids::volfields::PERBZVOL);
      cellParams[CellParams::EXVOL]                             = thisCellData->at(fsgrids::volfields::EXVOL);
      cellParams[CellParams::EYVOL]                             = thisCellData->at(fsgrids::volfields::EYVOL);
      cellParams[CellParams::EZVOL]                             = thisCellData->at(fsgrids::volfields::EZVOL);
      mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBXVOLdy] = thisCellData->at(fsgrids::volfields::dPERBXVOLdy);
      mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBXVOLdz] = thisCellData->at(fsgrids::volfields::dPERBXVOLdz);
      mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBYVOLdx] = thisCellData->at(fsgrids::volfields::dPERBYVOLdx);
      mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBYVOLdz] = thisCellData->at(fsgrids::volfields::dPERBYVOLdz);
      mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBZVOLdx] = thisCellData->at(fsgrids::volfields::dPERBZVOLdx);
      mpiGrid[i]->derivativesBVOL[bvolderivatives::dPERBZVOLdy] = thisCellData->at(fsgrids::volfields::dPERBZVOLdy);
   }

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
      //thisCellData->maxFsDt = mpiGrid[i]->get_cell_parameters()[CellParams::MAXFDT];
      thisCellData->maxFsDt = std::numeric_limits<Real>::max();
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
      mpiGrid[i]->get_cell_parameters()[CellParams::MAXFDT] = transferBuffer[count].maxFsDt;
      mpiGrid[i]->get_cell_parameters()[CellParams::FSGRID_RANK] = transferBuffer[count].fsGridRank;
      mpiGrid[i]->get_cell_parameters()[CellParams::FSGRID_BOUNDARYTYPE] = transferBuffer[count++].sysBoundaryFlag;
   }
}
