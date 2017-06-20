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
   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {
      auto cellParams = mpiGrid[cells[i]]->get_cell_parameters();
      std::array<Real, fsgrids::moments::N_MOMENTS>* thisCellData = &transferBuffer[i];
    
      if(!dt2) {
         thisCellData->at(fsgrids::moments::RHOM) = cellParams[CellParams::RHOM];
         thisCellData->at(fsgrids::moments::RHOQ) = cellParams[CellParams::RHOQ];
         thisCellData->at(fsgrids::moments::VX) = cellParams[CellParams::VX];
         thisCellData->at(fsgrids::moments::VY) = cellParams[CellParams::VY];
         thisCellData->at(fsgrids::moments::VZ) = cellParams[CellParams::VZ];
         thisCellData->at(fsgrids::moments::P_11) = cellParams[CellParams::P_11];
         thisCellData->at(fsgrids::moments::P_22) = cellParams[CellParams::P_22];
         thisCellData->at(fsgrids::moments::P_33) = cellParams[CellParams::P_33];
      } else {
         thisCellData->at(fsgrids::moments::RHOM) = cellParams[CellParams::RHOM_DT2];
         thisCellData->at(fsgrids::moments::RHOQ) = cellParams[CellParams::RHOQ_DT2];
         thisCellData->at(fsgrids::moments::VX) = cellParams[CellParams::VX_DT2];
         thisCellData->at(fsgrids::moments::VY) = cellParams[CellParams::VY_DT2];
         thisCellData->at(fsgrids::moments::VZ) = cellParams[CellParams::VZ_DT2];
         thisCellData->at(fsgrids::moments::P_11) = cellParams[CellParams::P_11_DT2];
         thisCellData->at(fsgrids::moments::P_22) = cellParams[CellParams::P_22_DT2];
         thisCellData->at(fsgrids::moments::P_33) = cellParams[CellParams::P_33_DT2];
      }
   }

   for(int i=0; i< cells.size(); i++) {
      momentsGrid.transferDataIn(cells[i] - 1, &transferBuffer[i]);
   }

   // Finish the actual transfer
   momentsGrid.finishTransfersIn();
}


void feedBgFieldsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
    const std::vector<CellID>& cells, FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& bgBGrid) {

  bgBGrid.setupForTransferIn(cells.size());

   // Setup transfer buffers
   std::vector< std::array<Real, fsgrids::bgbfield::N_BGB> > transferBuffer(cells.size());

   // Fill from cellParams
   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {
      auto cellParams = mpiGrid[cells[i]]->get_cell_parameters();
      auto derivatives = mpiGrid[cells[i]]->derivatives;
      auto volumeDerivatives = mpiGrid[cells[i]]->derivativesBVOL;
      std::array<Real, fsgrids::bgbfield::N_BGB>* thisCellData = &transferBuffer[i];

      thisCellData->at(fsgrids::bgbfield::BGBX) = cellParams[CellParams::BGBX];
      thisCellData->at(fsgrids::bgbfield::BGBY) = cellParams[CellParams::BGBY];
      thisCellData->at(fsgrids::bgbfield::BGBZ) = cellParams[CellParams::BGBZ];
      thisCellData->at(fsgrids::bgbfield::BGBXVOL) = cellParams[CellParams::BGBXVOL];
      thisCellData->at(fsgrids::bgbfield::BGBYVOL) = cellParams[CellParams::BGBYVOL];
      thisCellData->at(fsgrids::bgbfield::BGBZVOL) = cellParams[CellParams::BGBZVOL];

      thisCellData->at(fsgrids::bgbfield::dBGBxdy) = derivatives[fieldsolver::dBGBxdy];
      thisCellData->at(fsgrids::bgbfield::dBGBxdz) = derivatives[fieldsolver::dBGBxdz];
      thisCellData->at(fsgrids::bgbfield::dBGBydx) = derivatives[fieldsolver::dBGBydx];
      thisCellData->at(fsgrids::bgbfield::dBGBydz) = derivatives[fieldsolver::dBGBydz];
      thisCellData->at(fsgrids::bgbfield::dBGBzdx) = derivatives[fieldsolver::dBGBzdx];
      thisCellData->at(fsgrids::bgbfield::dBGBzdy) = derivatives[fieldsolver::dBGBzdy];

      thisCellData->at(fsgrids::bgbfield::dBGBXVOLdy) = volumeDerivatives[bvolderivatives::dBGBXVOLdy];
      thisCellData->at(fsgrids::bgbfield::dBGBXVOLdz) = volumeDerivatives[bvolderivatives::dBGBXVOLdz];
      thisCellData->at(fsgrids::bgbfield::dBGBYVOLdx) = volumeDerivatives[bvolderivatives::dBGBYVOLdx];
      thisCellData->at(fsgrids::bgbfield::dBGBYVOLdz) = volumeDerivatives[bvolderivatives::dBGBYVOLdz];
      thisCellData->at(fsgrids::bgbfield::dBGBZVOLdx) = volumeDerivatives[bvolderivatives::dBGBZVOLdx];
      thisCellData->at(fsgrids::bgbfield::dBGBZVOLdy) = volumeDerivatives[bvolderivatives::dBGBZVOLdy];
   }

   for(int i=0; i< cells.size(); i++) {
      bgBGrid.transferDataIn(cells[i] - 1, &transferBuffer[i]);
   }

   // Finish the actual transfer
   bgBGrid.finishTransfersIn();

}

void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells) {

   // Setup transfer buffers
   std::vector< std::array<Real, fsgrids::volfields::N_VOL> > transferBuffer(cells.size());

   // Setup transfer pointers
   volumeFieldsGrid.setupForTransferOut(cells.size());
   for(int i=0; i< cells.size(); i++) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[i];
      volumeFieldsGrid.transferDataOut(cells[i] - 1, thisCellData);
   }
   // Do the transfer
   volumeFieldsGrid.finishTransfersOut();

   // Distribute data from the transfer buffer back into the appropriate mpiGrid places
   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[i];
      auto cellParams = mpiGrid[cells[i]]->get_cell_parameters();

      cellParams[CellParams::PERBXVOL]                          = thisCellData->at(fsgrids::volfields::PERBXVOL);
      cellParams[CellParams::PERBYVOL]                          = thisCellData->at(fsgrids::volfields::PERBYVOL);
      cellParams[CellParams::PERBZVOL]                          = thisCellData->at(fsgrids::volfields::PERBZVOL);
      cellParams[CellParams::EXVOL]                             = thisCellData->at(fsgrids::volfields::EXVOL);
      cellParams[CellParams::EYVOL]                             = thisCellData->at(fsgrids::volfields::EYVOL);
      cellParams[CellParams::EZVOL]                             = thisCellData->at(fsgrids::volfields::EZVOL);
      mpiGrid[cells[i]]->derivativesBVOL[bvolderivatives::dPERBXVOLdy] = thisCellData->at(fsgrids::volfields::dPERBXVOLdy);
      mpiGrid[cells[i]]->derivativesBVOL[bvolderivatives::dPERBXVOLdz] = thisCellData->at(fsgrids::volfields::dPERBXVOLdz);
      mpiGrid[cells[i]]->derivativesBVOL[bvolderivatives::dPERBYVOLdx] = thisCellData->at(fsgrids::volfields::dPERBYVOLdx);
      mpiGrid[cells[i]]->derivativesBVOL[bvolderivatives::dPERBYVOLdz] = thisCellData->at(fsgrids::volfields::dPERBYVOLdz);
      mpiGrid[cells[i]]->derivativesBVOL[bvolderivatives::dPERBZVOLdx] = thisCellData->at(fsgrids::volfields::dPERBZVOLdx);
      mpiGrid[cells[i]]->derivativesBVOL[bvolderivatives::dPERBZVOLdy] = thisCellData->at(fsgrids::volfields::dPERBZVOLdy);
   }

}


void getDerivativesFromFsGrid(FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dperbGrid,
                          FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dmomentsGrid,
                          FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& bgbfieldGrid,
                          dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& cells) {

   // Setup transfer buffers
   std::vector< std::array<Real, fsgrids::dperb::N_DPERB> > dperbTransferBuffer(cells.size());
   std::vector< std::array<Real, fsgrids::dmoments::N_DMOMENTS> > dmomentsTransferBuffer(cells.size());
   std::vector< std::array<Real, fsgrids::bgbfield::N_BGB> > bgbfieldTransferBuffer(cells.size());

   // Transfer dperbGrid data
   dperbGrid.setupForTransferOut(cells.size());
   for(int i=0; i< cells.size(); i++) {
      std::array<Real, fsgrids::dperb::N_DPERB>* thisCellData = &dperbTransferBuffer[i];
      dperbGrid.transferDataOut(cells[i] - 1, thisCellData);
   }
   // Do the transfer
   dperbGrid.finishTransfersOut();

   // Transfer dmomentsGrid data
   dmomentsGrid.setupForTransferOut(cells.size());
   for(int i=0; i< cells.size(); i++) {
      std::array<Real, fsgrids::dmoments::N_DMOMENTS>* thisCellData = &dmomentsTransferBuffer[i];
      dmomentsGrid.transferDataOut(cells[i] - 1, thisCellData);
   }
   // Do the transfer
   dmomentsGrid.finishTransfersOut();

   // Transfer bgbfieldGrid data
   bgbfieldGrid.setupForTransferOut(cells.size());
   for(int i=0; i< cells.size(); i++) {
      std::array<Real, fsgrids::bgbfield::N_BGB>* thisCellData = &bgbfieldTransferBuffer[i];
      bgbfieldGrid.transferDataOut(cells[i] - 1, thisCellData);
   }
   // Do the transfer
   bgbfieldGrid.finishTransfersOut();

   // Distribute data from the transfer buffers back into the appropriate mpiGrid places
   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {
      std::array<Real, fsgrids::dperb::N_DPERB>* dperb = &dperbTransferBuffer[i];
      std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments = &dmomentsTransferBuffer[i];
      std::array<Real, fsgrids::bgbfield::N_BGB>* bgbfield = &bgbfieldTransferBuffer[i];
      auto cellParams = mpiGrid[cells[i]]->get_cell_parameters();

      mpiGrid[cells[i]]->derivatives[fieldsolver::drhomdx] = dmoments->at(fsgrids::dmoments::drhomdx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::drhomdy] = dmoments->at(fsgrids::dmoments::drhomdy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::drhomdz] = dmoments->at(fsgrids::dmoments::drhomdz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::drhoqdx] = dmoments->at(fsgrids::dmoments::drhoqdx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::drhoqdy] = dmoments->at(fsgrids::dmoments::drhoqdy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::drhoqdz] = dmoments->at(fsgrids::dmoments::drhoqdz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp11dx] = dmoments->at(fsgrids::dmoments::dp11dx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp11dy] = dmoments->at(fsgrids::dmoments::dp11dy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp11dz] = dmoments->at(fsgrids::dmoments::dp11dz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp22dx] = dmoments->at(fsgrids::dmoments::dp22dx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp22dy] = dmoments->at(fsgrids::dmoments::dp22dy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp22dz] = dmoments->at(fsgrids::dmoments::dp22dz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp33dx] = dmoments->at(fsgrids::dmoments::dp33dx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp33dy] = dmoments->at(fsgrids::dmoments::dp33dy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dp33dz] = dmoments->at(fsgrids::dmoments::dp33dz);

      mpiGrid[cells[i]]->derivatives[fieldsolver::dVxdx] = dmoments->at(fsgrids::dmoments::dVxdx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVxdy] = dmoments->at(fsgrids::dmoments::dVxdy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVxdz] = dmoments->at(fsgrids::dmoments::dVxdz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVydx] = dmoments->at(fsgrids::dmoments::dVydx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVydy] = dmoments->at(fsgrids::dmoments::dVydy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVydz] = dmoments->at(fsgrids::dmoments::dVydz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVzdx] = dmoments->at(fsgrids::dmoments::dVzdx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVzdy] = dmoments->at(fsgrids::dmoments::dVzdy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dVzdz] = dmoments->at(fsgrids::dmoments::dVzdz);

      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBxdy] = dperb->at(fsgrids::dperb::dPERBxdy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBxdz] = dperb->at(fsgrids::dperb::dPERBxdz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBydx] = dperb->at(fsgrids::dperb::dPERBydx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBydz] = dperb->at(fsgrids::dperb::dPERBydz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBzdx] = dperb->at(fsgrids::dperb::dPERBzdx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBzdy] = dperb->at(fsgrids::dperb::dPERBzdy);

      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBxdyy] = dperb->at(fsgrids::dperb::dPERBxdyy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBxdzz] = dperb->at(fsgrids::dperb::dPERBxdzz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBydxx] = dperb->at(fsgrids::dperb::dPERBydxx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBydzz] = dperb->at(fsgrids::dperb::dPERBydzz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBzdxx] = dperb->at(fsgrids::dperb::dPERBzdxx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBzdyy] = dperb->at(fsgrids::dperb::dPERBzdyy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBxdyz] = dperb->at(fsgrids::dperb::dPERBxdyz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBydxz] = dperb->at(fsgrids::dperb::dPERBydxz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dPERBzdxy] = dperb->at(fsgrids::dperb::dPERBzdxy);

      mpiGrid[cells[i]]->derivatives[fieldsolver::dBGBxdy] = bgbfield->at(fsgrids::bgbfield::dBGBxdy);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dBGBxdz] = bgbfield->at(fsgrids::bgbfield::dBGBxdz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dBGBydx] = bgbfield->at(fsgrids::bgbfield::dBGBydx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dBGBydz] = bgbfield->at(fsgrids::bgbfield::dBGBydz);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dBGBzdx] = bgbfield->at(fsgrids::bgbfield::dBGBzdx);
      mpiGrid[cells[i]]->derivatives[fieldsolver::dBGBzdy] = bgbfield->at(fsgrids::bgbfield::dBGBzdy);
   }

}
    

void setupTechnicalFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, FsGrid< fsgrids::technical, 2>& technicalGrid) {

   technicalGrid.setupForTransferIn(cells.size());

   // Fill the transfer buffers from the spatial cell structs
   std::vector<fsgrids::technical> transferBuffer(cells.size());

   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {

      fsgrids::technical* thisCellData = &transferBuffer[i];
      // Data needs to be collected from some different places for this grid.
      thisCellData->sysBoundaryFlag = mpiGrid[cells[i]]->sysBoundaryFlag;
      thisCellData->sysBoundaryLayer = mpiGrid[cells[i]]->sysBoundaryLayer;
      //thisCellData->maxFsDt = mpiGrid[i]->get_cell_parameters()[CellParams::MAXFDT];
      thisCellData->maxFsDt = std::numeric_limits<Real>::max();
   }
   for(int i=0; i< cells.size(); i++) {
      technicalGrid.transferDataIn(cells[i] - 1,&transferBuffer[i]);
   }

   technicalGrid.finishTransfersIn();
}

void getFsGridMaxDt(FsGrid< fsgrids::technical, 2>& technicalGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells) {

   technicalGrid.setupForTransferOut(cells.size());

   // Buffer to store contents of the grid
   std::vector<fsgrids::technical> transferBuffer(cells.size());

   for(int i=0; i< cells.size(); i++) {
      fsgrids::technical* thisCellData = &transferBuffer[i];
      technicalGrid.transferDataOut(cells[i] - 1, thisCellData);
   }

   technicalGrid.finishTransfersOut();

   // After the transfer is completed, stuff the recieved maxFDt into the cells.
   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {
      mpiGrid[cells[i]]->get_cell_parameters()[CellParams::MAXFDT] = transferBuffer[i].maxFsDt;
      mpiGrid[cells[i]]->get_cell_parameters()[CellParams::FSGRID_RANK] = transferBuffer[i].fsGridRank;
      mpiGrid[cells[i]]->get_cell_parameters()[CellParams::FSGRID_BOUNDARYTYPE] = transferBuffer[i].sysBoundaryFlag;
   }
}

