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
      std::array<Real, fsgrids::moments::N_MOMENTS>* thisCellData = &transferBuffer[count];
    
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

      count++;

      momentsGrid.transferDataIn(i - 1, thisCellData);
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
   size_t count=0;
   for(CellID i : cells) {
      auto cellParams = mpiGrid[i]->get_cell_parameters();
      auto derivatives = mpiGrid[i]->derivatives;
      auto volumeDerivatives = mpiGrid[i]->derivativesBVOL;
      std::array<Real, fsgrids::bgbfield::N_BGB>* thisCellData = &transferBuffer[count];

      thisCellData->at(fsgrids::bgbfield::BGBX) = cellParams[CellParams::BGBX];
      thisCellData->at(fsgrids::bgbfield::BGBY) = cellParams[CellParams::BGBY];
      thisCellData->at(fsgrids::bgbfield::BGBZ) = cellParams[CellParams::BGBZ];
      thisCellData->at(fsgrids::bgbfield::BGBXVOL) = cellParams[CellParams::BGBXVOL];
      thisCellData->at(fsgrids::bgbfield::BGBYVOL) = cellParams[CellParams::BGBYVOL];
      thisCellData->at(fsgrids::bgbfield::BGBZVOL) = cellParams[CellParams::BGBZVOL];
      thisCellData->at(fsgrids::bgbfield::BGBX_000_010) = cellParams[CellParams::BGBX_000_010];
      thisCellData->at(fsgrids::bgbfield::BGBX_100_110) = cellParams[CellParams::BGBX_100_110];
      thisCellData->at(fsgrids::bgbfield::BGBX_001_011) = cellParams[CellParams::BGBX_001_011];
      thisCellData->at(fsgrids::bgbfield::BGBX_101_111) = cellParams[CellParams::BGBX_101_111];
      thisCellData->at(fsgrids::bgbfield::BGBX_000_001) = cellParams[CellParams::BGBX_000_001];
      thisCellData->at(fsgrids::bgbfield::BGBX_100_101) = cellParams[CellParams::BGBX_100_101];
      thisCellData->at(fsgrids::bgbfield::BGBX_010_011) = cellParams[CellParams::BGBX_010_011];
      thisCellData->at(fsgrids::bgbfield::BGBX_110_111) = cellParams[CellParams::BGBX_110_111];
      thisCellData->at(fsgrids::bgbfield::BGBY_000_100) = cellParams[CellParams::BGBY_000_100];
      thisCellData->at(fsgrids::bgbfield::BGBY_010_110) = cellParams[CellParams::BGBY_010_110];
      thisCellData->at(fsgrids::bgbfield::BGBY_001_101) = cellParams[CellParams::BGBY_001_101];
      thisCellData->at(fsgrids::bgbfield::BGBY_011_111) = cellParams[CellParams::BGBY_011_111];
      thisCellData->at(fsgrids::bgbfield::BGBY_000_001) = cellParams[CellParams::BGBY_000_001];
      thisCellData->at(fsgrids::bgbfield::BGBY_100_101) = cellParams[CellParams::BGBY_100_101];
      thisCellData->at(fsgrids::bgbfield::BGBY_010_011) = cellParams[CellParams::BGBY_010_011];
      thisCellData->at(fsgrids::bgbfield::BGBY_110_111) = cellParams[CellParams::BGBY_110_111];
      thisCellData->at(fsgrids::bgbfield::BGBZ_000_100) = cellParams[CellParams::BGBZ_000_100];
      thisCellData->at(fsgrids::bgbfield::BGBZ_010_110) = cellParams[CellParams::BGBZ_010_110];
      thisCellData->at(fsgrids::bgbfield::BGBZ_001_101) = cellParams[CellParams::BGBZ_001_101];
      thisCellData->at(fsgrids::bgbfield::BGBZ_011_111) = cellParams[CellParams::BGBZ_011_111];
      thisCellData->at(fsgrids::bgbfield::BGBZ_000_010) = cellParams[CellParams::BGBZ_000_010];
      thisCellData->at(fsgrids::bgbfield::BGBZ_100_110) = cellParams[CellParams::BGBZ_100_110];
      thisCellData->at(fsgrids::bgbfield::BGBZ_001_011) = cellParams[CellParams::BGBZ_001_011];
      thisCellData->at(fsgrids::bgbfield::BGBZ_101_111) = cellParams[CellParams::BGBZ_101_111];

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

      bgBGrid.transferDataIn(i-1, thisCellData);
      count++;
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
   size_t count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[count];

      volumeFieldsGrid.transferDataOut(i - 1, thisCellData);
      count++;
   }
   // Do the transfer
   volumeFieldsGrid.finishTransfersOut();

   // Distribute data from the transfer buffer back into the appropriate mpiGrid places
   count = 0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[count];
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

      count++;
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
   size_t count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::dperb::N_DPERB>* thisCellData = &dperbTransferBuffer[count];

      dperbGrid.transferDataOut(i - 1, thisCellData);
      count++;
   }
   // Do the transfer
   dperbGrid.finishTransfersOut();

   // Transfer dmomentsGrid data
   dmomentsGrid.setupForTransferOut(cells.size());
   count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::dmoments::N_DMOMENTS>* thisCellData = &dmomentsTransferBuffer[count];

      dmomentsGrid.transferDataOut(i - 1, thisCellData);
      count++;
   }
   // Do the transfer
   dmomentsGrid.finishTransfersOut();

   // Transfer bgbfieldGrid data
   bgbfieldGrid.setupForTransferOut(cells.size());
   count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::bgbfield::N_BGB>* thisCellData = &bgbfieldTransferBuffer[count];

      bgbfieldGrid.transferDataOut(i - 1, thisCellData);
      count++;
   }
   // Do the transfer
   bgbfieldGrid.finishTransfersOut();

   // Distribute data from the transfer buffers back into the appropriate mpiGrid places
   count=0;
   for(CellID i : cells) {
      std::array<Real, fsgrids::dperb::N_DPERB>* dperb = &dperbTransferBuffer[count];
      std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments = &dmomentsTransferBuffer[count];
      std::array<Real, fsgrids::bgbfield::N_BGB>* bgbfield = &bgbfieldTransferBuffer[count];
      auto cellParams = mpiGrid[i]->get_cell_parameters();

      mpiGrid[i]->derivatives[fieldsolver::drhodx] = dmoments->at(fsgrids::dmoments::drhodx);
      mpiGrid[i]->derivatives[fieldsolver::drhody] = dmoments->at(fsgrids::dmoments::drhody);
      mpiGrid[i]->derivatives[fieldsolver::drhodz] = dmoments->at(fsgrids::dmoments::drhodz);
      mpiGrid[i]->derivatives[fieldsolver::dp11dx] = dmoments->at(fsgrids::dmoments::dp11dx);
      mpiGrid[i]->derivatives[fieldsolver::dp11dy] = dmoments->at(fsgrids::dmoments::dp11dy);
      mpiGrid[i]->derivatives[fieldsolver::dp11dz] = dmoments->at(fsgrids::dmoments::dp11dz);
      mpiGrid[i]->derivatives[fieldsolver::dp22dx] = dmoments->at(fsgrids::dmoments::dp22dx);
      mpiGrid[i]->derivatives[fieldsolver::dp22dy] = dmoments->at(fsgrids::dmoments::dp22dy);
      mpiGrid[i]->derivatives[fieldsolver::dp22dz] = dmoments->at(fsgrids::dmoments::dp22dz);
      mpiGrid[i]->derivatives[fieldsolver::dp33dx] = dmoments->at(fsgrids::dmoments::dp33dx);
      mpiGrid[i]->derivatives[fieldsolver::dp33dy] = dmoments->at(fsgrids::dmoments::dp33dy);
      mpiGrid[i]->derivatives[fieldsolver::dp33dz] = dmoments->at(fsgrids::dmoments::dp33dz);

      mpiGrid[i]->derivatives[fieldsolver::dVxdx] = dmoments->at(fsgrids::dmoments::dVxdx);
      mpiGrid[i]->derivatives[fieldsolver::dVxdy] = dmoments->at(fsgrids::dmoments::dVxdy);
      mpiGrid[i]->derivatives[fieldsolver::dVxdz] = dmoments->at(fsgrids::dmoments::dVxdz);
      mpiGrid[i]->derivatives[fieldsolver::dVydx] = dmoments->at(fsgrids::dmoments::dVydx);
      mpiGrid[i]->derivatives[fieldsolver::dVydy] = dmoments->at(fsgrids::dmoments::dVydy);
      mpiGrid[i]->derivatives[fieldsolver::dVydz] = dmoments->at(fsgrids::dmoments::dVydz);
      mpiGrid[i]->derivatives[fieldsolver::dVzdx] = dmoments->at(fsgrids::dmoments::dVzdx);
      mpiGrid[i]->derivatives[fieldsolver::dVzdy] = dmoments->at(fsgrids::dmoments::dVzdy);
      mpiGrid[i]->derivatives[fieldsolver::dVzdz] = dmoments->at(fsgrids::dmoments::dVzdz);

      mpiGrid[i]->derivatives[fieldsolver::dPERBxdy] = dperb->at(fsgrids::dperb::dPERBxdy);
      mpiGrid[i]->derivatives[fieldsolver::dPERBxdz] = dperb->at(fsgrids::dperb::dPERBxdz);
      mpiGrid[i]->derivatives[fieldsolver::dPERBydx] = dperb->at(fsgrids::dperb::dPERBydx);
      mpiGrid[i]->derivatives[fieldsolver::dPERBydz] = dperb->at(fsgrids::dperb::dPERBydz);
      mpiGrid[i]->derivatives[fieldsolver::dPERBzdx] = dperb->at(fsgrids::dperb::dPERBzdx);
      mpiGrid[i]->derivatives[fieldsolver::dPERBzdy] = dperb->at(fsgrids::dperb::dPERBzdy);

      mpiGrid[i]->derivatives[fieldsolver::dPERBxdyy] = dperb->at(fsgrids::dperb::dPERBxdyy);
      mpiGrid[i]->derivatives[fieldsolver::dPERBxdzz] = dperb->at(fsgrids::dperb::dPERBxdzz);
      mpiGrid[i]->derivatives[fieldsolver::dPERBydxx] = dperb->at(fsgrids::dperb::dPERBydxx);
      mpiGrid[i]->derivatives[fieldsolver::dPERBydzz] = dperb->at(fsgrids::dperb::dPERBydzz);
      mpiGrid[i]->derivatives[fieldsolver::dPERBzdxx] = dperb->at(fsgrids::dperb::dPERBzdxx);
      mpiGrid[i]->derivatives[fieldsolver::dPERBzdyy] = dperb->at(fsgrids::dperb::dPERBzdyy);
      mpiGrid[i]->derivatives[fieldsolver::dPERBxdyz] = dperb->at(fsgrids::dperb::dPERBxdyz);
      mpiGrid[i]->derivatives[fieldsolver::dPERBydxz] = dperb->at(fsgrids::dperb::dPERBydxz);
      mpiGrid[i]->derivatives[fieldsolver::dPERBzdxy] = dperb->at(fsgrids::dperb::dPERBzdxy);

      mpiGrid[i]->derivatives[fieldsolver::dBGBxdy] = bgbfield->at(fsgrids::bgbfield::dBGBxdy);
      mpiGrid[i]->derivatives[fieldsolver::dBGBxdz] = bgbfield->at(fsgrids::bgbfield::dBGBxdz);
      mpiGrid[i]->derivatives[fieldsolver::dBGBydx] = bgbfield->at(fsgrids::bgbfield::dBGBydx);
      mpiGrid[i]->derivatives[fieldsolver::dBGBydz] = bgbfield->at(fsgrids::bgbfield::dBGBydz);
      mpiGrid[i]->derivatives[fieldsolver::dBGBzdx] = bgbfield->at(fsgrids::bgbfield::dBGBzdx);
      mpiGrid[i]->derivatives[fieldsolver::dBGBzdy] = bgbfield->at(fsgrids::bgbfield::dBGBzdy);

      count++;
   }

}
    

void setupTechnicalFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, FsGrid< fsgrids::technical, 2>& technicalGrid) {

   technicalGrid.setupForTransferIn(cells.size());

   // Fill the transfer buffers from the spatial cell structs
   std::vector<fsgrids::technical> transferBuffer(cells.size());

   size_t count=0;
   for(CellID i : cells) {

      fsgrids::technical* thisCellData = &transferBuffer[count];
      // Data needs to be collected from some different places for this grid.
      thisCellData->sysBoundaryFlag = mpiGrid[i]->sysBoundaryFlag;
      thisCellData->sysBoundaryLayer = mpiGrid[i]->sysBoundaryLayer;
      //thisCellData->maxFsDt = mpiGrid[i]->get_cell_parameters()[CellParams::MAXFDT];
      thisCellData->maxFsDt = std::numeric_limits<Real>::max();
      technicalGrid.transferDataIn(i - 1,thisCellData);

      count++;
   }

   technicalGrid.finishTransfersIn();
   
   // Checking that spatial cells are cubic, otherwise field solver is incorrect (cf. derivatives in E, Hall term)
   if((abs((technicalGrid.DX-technicalGrid.DY)/technicalGrid.DX) > 0.001) ||
      (abs((technicalGrid.DX-technicalGrid.DZ)/technicalGrid.DX) > 0.001) ||
      (abs((technicalGrid.DY-technicalGrid.DZ)/technicalGrid.DY) > 0.001)) {
      std::cerr << "WARNING: Your spatial cells seem not to be cubic. However the field solver is assuming them to be. Use at your own risk and responsibility!" << std::endl;
      }
   
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
      mpiGrid[i]->get_cell_parameters()[CellParams::FSGRID_BOUNDARYTYPE] = transferBuffer[count].sysBoundaryFlag;

      count++;
   }
}

