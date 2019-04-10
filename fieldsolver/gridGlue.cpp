#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "../grid.h"
#include "../spatial_cell.hpp"
#include "../definitions.h"
#include "../common.h"
#include "gridGlue.hpp"

/*
Calculate the number of cells on the maximum refinement level overlapping the list of dccrg cells in cells.
*/
int getNumberOfCellsOnMaxRefLvl(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& cells) {

   int nCells = 0;
   auto maxRefLvl = mpiGrid.mapping.get_maximum_refinement_level();
   
   for (auto cellid : cells) {
      auto refLvl = mpiGrid.get_refinement_level(cellid);
      nCells += pow(pow(2,maxRefLvl-refLvl),3);
   }

   return nCells;
   
}



void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid, bool dt2 /*=false*/) {

   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   momentsGrid.setupForTransferIn(nCells);

   std::vector< std::array<Real, fsgrids::moments::N_MOMENTS> > transferBuffer(cells.size());
   
   // Fill from cellParams
#pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {
      CellID dccrgId = cells[i];
      auto cellParams = mpiGrid[dccrgId]->get_cell_parameters();
      
      std::array<Real, fsgrids::moments::N_MOMENTS>* thisCellData = &transferBuffer[i];
      
      if(dt2) {
         thisCellData->at(fsgrids::moments::RHOM) = cellParams[CellParams::RHOM_DT2];
         thisCellData->at(fsgrids::moments::RHOQ) = cellParams[CellParams::RHOQ_DT2];
         thisCellData->at(fsgrids::moments::VX)   = cellParams[CellParams::VX_DT2];
         thisCellData->at(fsgrids::moments::VY)   = cellParams[CellParams::VY_DT2];
         thisCellData->at(fsgrids::moments::VZ)   = cellParams[CellParams::VZ_DT2];
         thisCellData->at(fsgrids::moments::P_11) = cellParams[CellParams::P_11_DT2];
         thisCellData->at(fsgrids::moments::P_22) = cellParams[CellParams::P_22_DT2];
         thisCellData->at(fsgrids::moments::P_33) = cellParams[CellParams::P_33_DT2];
      } else {
         thisCellData->at(fsgrids::moments::RHOM) = cellParams[CellParams::RHOM];
         thisCellData->at(fsgrids::moments::RHOQ) = cellParams[CellParams::RHOQ];
         thisCellData->at(fsgrids::moments::VX)   = cellParams[CellParams::VX];
         thisCellData->at(fsgrids::moments::VY)   = cellParams[CellParams::VY];
         thisCellData->at(fsgrids::moments::VZ)   = cellParams[CellParams::VZ];
         thisCellData->at(fsgrids::moments::P_11) = cellParams[CellParams::P_11];
         thisCellData->at(fsgrids::moments::P_22) = cellParams[CellParams::P_22];
         thisCellData->at(fsgrids::moments::P_33) = cellParams[CellParams::P_33];
      }
   }



   for (uint i = 0;i < cells.size(); ++i) {
      CellID dccrgId = cells[i];
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
      for (auto fsgridId : fsgridIds) {
         momentsGrid.transferDataIn(fsgridId, &transferBuffer[i]);
      }
   }

   // Finish the actual transfer
   momentsGrid.finishTransfersIn();
}


void feedBgFieldsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
    const std::vector<CellID>& cells, FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& bgBGrid) {

   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   bgBGrid.setupForTransferIn(nCells);

   // Setup transfer buffers
   std::vector< std::array<Real, fsgrids::bgbfield::N_BGB> > transferBuffer(cells.size());
   
   // Fill from cellParams
   // We only need to read data once per dccrg cell here
#pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {
      CellID dccrgId = cells[i];
      auto cellParams = mpiGrid[dccrgId]->get_cell_parameters();
      auto derivatives = mpiGrid[dccrgId]->derivatives;
      auto volumeDerivatives = mpiGrid[dccrgId]->derivativesBVOL;

      //      std::cout << "I am at line " << __LINE__ << " of " << __FILE__ << std::endl;
               
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

   // Copy data into each fsgrid cell overlapping the dccrg cell
   for (uint i = 0; i < cells.size(); ++i) {
      CellID dccrgId = cells[i];
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
      for (auto fsgridId : fsgridIds) {
         bgBGrid.transferDataIn(fsgridId, &transferBuffer[i]);
      }
   }
   
   // Finish the actual transfer
   bgBGrid.finishTransfersIn();

}

void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volumeFieldsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells) {


   // Setup transfer buffers
   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   std::vector< std::array<Real, fsgrids::volfields::N_VOL> > transferBuffer(nCells);
   std::vector< std::array<Real, fsgrids::volfields::N_VOL>*> transferBufferPointer;

   // Setup transfer pointers
   volumeFieldsGrid.setupForTransferOut(nCells);
   int k = 0;
   for(auto dccrgId : cells) {
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
      // Store a pointer to the first fsgrid cell that maps to each dccrg Id
      transferBufferPointer.push_back(&transferBuffer[k]);
      for (auto fsgridId : fsgridIds) {
         std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = &transferBuffer[k++];
         volumeFieldsGrid.transferDataOut(fsgridId, thisCellData);
      }
   }
   // Do the transfer
   volumeFieldsGrid.finishTransfersOut();

   // Build a list of index pairs to cellparams and fsgrid
   std::vector<std::pair<int,int>> iCellParams;
   iCellParams.reserve(6);
   iCellParams.push_back(std::make_pair(CellParams::PERBXVOL, fsgrids::volfields::PERBXVOL));
   iCellParams.push_back(std::make_pair(CellParams::PERBYVOL, fsgrids::volfields::PERBYVOL));
   iCellParams.push_back(std::make_pair(CellParams::PERBZVOL, fsgrids::volfields::PERBZVOL));
   iCellParams.push_back(std::make_pair(CellParams::EXVOL,    fsgrids::volfields::EXVOL));
   iCellParams.push_back(std::make_pair(CellParams::EYVOL,    fsgrids::volfields::EYVOL));
   iCellParams.push_back(std::make_pair(CellParams::EZVOL,    fsgrids::volfields::EZVOL));

   // Build lists of index pairs to dccrg and fsgrid
   std::vector<std::pair<int,int>> iDerivativesBVOL;
   iDerivativesBVOL.reserve(6);
   iDerivativesBVOL.push_back(std::make_pair(bvolderivatives::dPERBXVOLdy, fsgrids::volfields::dPERBXVOLdy));
   iDerivativesBVOL.push_back(std::make_pair(bvolderivatives::dPERBXVOLdz, fsgrids::volfields::dPERBXVOLdz));
   iDerivativesBVOL.push_back(std::make_pair(bvolderivatives::dPERBYVOLdx, fsgrids::volfields::dPERBYVOLdx));
   iDerivativesBVOL.push_back(std::make_pair(bvolderivatives::dPERBYVOLdz, fsgrids::volfields::dPERBYVOLdz));
   iDerivativesBVOL.push_back(std::make_pair(bvolderivatives::dPERBZVOLdx, fsgrids::volfields::dPERBZVOLdx));
   iDerivativesBVOL.push_back(std::make_pair(bvolderivatives::dPERBZVOLdy, fsgrids::volfields::dPERBZVOLdy));
   
   // Distribute data from the transfer buffer back into the appropriate mpiGrid places
   #pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {

      int dccrgId = cells[i];
      auto cellParams = mpiGrid[dccrgId]->get_cell_parameters();

      // Calculate the number of fsgrid cells we need to average into the current dccrg cell
      int nCells = pow(pow(2,mpiGrid.mapping.get_maximum_refinement_level() - mpiGrid.mapping.get_refinement_level(dccrgId)),3);

      // TODO: Could optimize here by adding a separate branch for nCells == 1 with direct assignment of the value
      // Could also do the average in a temporary value and only access grid structure once.
      
      // Initialize values to 0
      for (auto j : iCellParams)      cellParams[j.first]                        = 0.0;
      for (auto j : iDerivativesBVOL) mpiGrid[dccrgId]->derivativesBVOL[j.first] = 0.0;
                  
      for(int iCell = 0; iCell < nCells; ++iCell) {
         // The fsgrid cells that cover the i'th dccrg cell are pointed at by
         // transferBufferPointer[i] ... transferBufferPointer[i] + nCell. We want to average
         // over all of them to get the value for the dccrg cell
         std::array<Real, fsgrids::volfields::N_VOL>* thisCellData = transferBufferPointer[i] + iCell;

         for (auto j : iCellParams)      cellParams[j.first]                        += thisCellData->at(j.second);
         for (auto j : iDerivativesBVOL) mpiGrid[dccrgId]->derivativesBVOL[j.first] += thisCellData->at(j.second);
      }

      // Divide by the number of cells to get the average
      for (auto j : iCellParams)      cellParams[j.first]                        /= nCells;
      for (auto j : iDerivativesBVOL) mpiGrid[dccrgId]->derivativesBVOL[j.first] /= nCells;

   }
}


void getDerivativesFromFsGrid(FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dperbGrid,
                          FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dmomentsGrid,
                          FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& bgbfieldGrid,
                          dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& cells) {

   // Setup transfer buffers
   int nCellsOnMaxRefLvl = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);
   std::vector< std::array<Real, fsgrids::dperb::N_DPERB> > dperbTransferBuffer(nCellsOnMaxRefLvl);
   std::vector< std::array<Real, fsgrids::dmoments::N_DMOMENTS> > dmomentsTransferBuffer(nCellsOnMaxRefLvl);
   std::vector< std::array<Real, fsgrids::bgbfield::N_BGB> > bgbfieldTransferBuffer(nCellsOnMaxRefLvl);
   
   std::vector< std::array<Real, fsgrids::dperb::N_DPERB>*> dperbTransferBufferPointer;
   std::vector< std::array<Real, fsgrids::dmoments::N_DMOMENTS>*> dmomentsTransferBufferPointer;
   std::vector< std::array<Real, fsgrids::bgbfield::N_BGB>*> bgbfieldTransferBufferPointer;
   
   dperbGrid.setupForTransferOut(nCellsOnMaxRefLvl);
   dmomentsGrid.setupForTransferOut(nCellsOnMaxRefLvl);
   bgbfieldGrid.setupForTransferOut(nCellsOnMaxRefLvl);
   
   int k = 0;
   for (auto dccrgId : cells) {

      // Assuming same local size in all fsgrids
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgId);
      // Store a pointer to the first fsgrid cell that maps to each dccrg Id
      dperbTransferBufferPointer.push_back(&dperbTransferBuffer[k]);
      dmomentsTransferBufferPointer.push_back(&dmomentsTransferBuffer[k]);
      bgbfieldTransferBufferPointer.push_back(&bgbfieldTransferBuffer[k]);

      for (auto fsgridId : fsgridIds) {
      
         std::array<Real, fsgrids::dperb::N_DPERB>* dperbCellData = &dperbTransferBuffer[k];
         dperbGrid.transferDataOut(fsgridId, dperbCellData);
         std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmomentsCellData = &dmomentsTransferBuffer[k];
         dmomentsGrid.transferDataOut(fsgridId, dmomentsCellData);
         std::array<Real, fsgrids::bgbfield::N_BGB>* bgbfieldCellData = &bgbfieldTransferBuffer[k++];
         bgbfieldGrid.transferDataOut(fsgridId, bgbfieldCellData);
      }
   }
   
   // Do the transfer
   dperbGrid.finishTransfersOut();
   dmomentsGrid.finishTransfersOut();
   bgbfieldGrid.finishTransfersOut();

   std::vector<std::pair<int,int>> iDmoments;
   std::vector<std::pair<int,int>> iDperb;
   std::vector<std::pair<int,int>> iBgbfield;
   iDmoments.reserve(24);
   iDmoments.push_back(std::make_pair(fieldsolver::drhomdx, fsgrids::dmoments::drhomdx));
   iDmoments.push_back(std::make_pair(fieldsolver::drhomdy, fsgrids::dmoments::drhomdy));
   iDmoments.push_back(std::make_pair(fieldsolver::drhomdz, fsgrids::dmoments::drhomdz));
   iDmoments.push_back(std::make_pair(fieldsolver::drhoqdx, fsgrids::dmoments::drhoqdx));
   iDmoments.push_back(std::make_pair(fieldsolver::drhoqdy, fsgrids::dmoments::drhoqdy));
   iDmoments.push_back(std::make_pair(fieldsolver::drhoqdz, fsgrids::dmoments::drhoqdz));
   iDmoments.push_back(std::make_pair(fieldsolver::dp11dx , fsgrids::dmoments::dp11dx ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp11dy , fsgrids::dmoments::dp11dy ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp11dz , fsgrids::dmoments::dp11dz ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp22dx , fsgrids::dmoments::dp22dx ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp22dy , fsgrids::dmoments::dp22dy ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp22dz , fsgrids::dmoments::dp22dz ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp33dx , fsgrids::dmoments::dp33dx ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp33dy , fsgrids::dmoments::dp33dy ));
   iDmoments.push_back(std::make_pair(fieldsolver::dp33dz , fsgrids::dmoments::dp33dz ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVxdx  , fsgrids::dmoments::dVxdx  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVxdy  , fsgrids::dmoments::dVxdy  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVxdz  , fsgrids::dmoments::dVxdz  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVydx  , fsgrids::dmoments::dVydx  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVydy  , fsgrids::dmoments::dVydy  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVydz  , fsgrids::dmoments::dVydz  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVzdx  , fsgrids::dmoments::dVzdx  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVzdy  , fsgrids::dmoments::dVzdy  ));
   iDmoments.push_back(std::make_pair(fieldsolver::dVzdz  , fsgrids::dmoments::dVzdz  ));

   iDperb.reserve(15);
   iDperb.push_back(std::make_pair(fieldsolver::dPERBxdy , fsgrids::dperb::dPERBxdy ));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBxdz , fsgrids::dperb::dPERBxdz ));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBydx , fsgrids::dperb::dPERBydx ));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBydz , fsgrids::dperb::dPERBydz ));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBzdx , fsgrids::dperb::dPERBzdx ));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBzdy , fsgrids::dperb::dPERBzdy ));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBxdyy, fsgrids::dperb::dPERBxdyy));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBxdzz, fsgrids::dperb::dPERBxdzz));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBydxx, fsgrids::dperb::dPERBydxx));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBydzz, fsgrids::dperb::dPERBydzz));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBzdxx, fsgrids::dperb::dPERBzdxx));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBzdyy, fsgrids::dperb::dPERBzdyy));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBxdyz, fsgrids::dperb::dPERBxdyz));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBydxz, fsgrids::dperb::dPERBydxz));
   iDperb.push_back(std::make_pair(fieldsolver::dPERBzdxy, fsgrids::dperb::dPERBzdxy));

   iBgbfield.reserve(6);
   iBgbfield.push_back(std::make_pair(fieldsolver::dBGBxdy, fsgrids::bgbfield::dBGBxdy));
   iBgbfield.push_back(std::make_pair(fieldsolver::dBGBxdz, fsgrids::bgbfield::dBGBxdz));
   iBgbfield.push_back(std::make_pair(fieldsolver::dBGBydx, fsgrids::bgbfield::dBGBydx));
   iBgbfield.push_back(std::make_pair(fieldsolver::dBGBydz, fsgrids::bgbfield::dBGBydz));
   iBgbfield.push_back(std::make_pair(fieldsolver::dBGBzdx, fsgrids::bgbfield::dBGBzdx));
   iBgbfield.push_back(std::make_pair(fieldsolver::dBGBzdy, fsgrids::bgbfield::dBGBzdy));
   
   // Distribute data from the transfer buffers back into the appropriate mpiGrid places
   #pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {

      int dccrgId = cells[i];

      // Calculate the number of fsgrid cells we need to average into the current dccrg cell
      auto refLvl = mpiGrid.mapping.get_refinement_level(dccrgId);
      int nCells = pow(pow(2,mpiGrid.mapping.get_maximum_refinement_level() - mpiGrid.mapping.get_refinement_level(dccrgId)),3);

      for (auto j : iDmoments) mpiGrid[dccrgId]->derivatives[j.first] = 0.0;
      for (auto j : iDperb   ) mpiGrid[dccrgId]->derivatives[j.first] = 0.0;
      for (auto j : iBgbfield) mpiGrid[dccrgId]->derivatives[j.first] = 0.0;
      
      for(int iCell = 0; iCell < nCells; ++iCell) {
         // The fsgrid cells that cover the i'th dccrg cell are pointed at by
         // transferBufferPointer[i] ... transferBufferPointer[i] + nCell. We want to average
         // over all of them to get the value for the dccrg cell
         
         std::array<Real, fsgrids::dperb::N_DPERB>* dperb    = dperbTransferBufferPointer[i] + iCell;
         std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments = dmomentsTransferBufferPointer[i] + iCell;
         std::array<Real, fsgrids::bgbfield::N_BGB>* bgbfield = bgbfieldTransferBufferPointer[i] + iCell;
      
         for (auto j : iDmoments) mpiGrid[dccrgId]->derivatives[j.first] += dmoments->at(j.second);
         for (auto j : iDperb   ) mpiGrid[dccrgId]->derivatives[j.first] += dperb   ->at(j.second);
         for (auto j : iBgbfield) mpiGrid[dccrgId]->derivatives[j.first] += bgbfield->at(j.second);
      }

      for (auto j : iDmoments) mpiGrid[dccrgId]->derivatives[j.first] /= nCells;
      for (auto j : iDperb   ) mpiGrid[dccrgId]->derivatives[j.first] /= nCells;
      for (auto j : iBgbfield) mpiGrid[dccrgId]->derivatives[j.first] /= nCells;

   }
}

bool belongsToLayer(const int layer, const int x, const int y, const int z,
                    FsGrid< fsgrids::technical, 2>& technicalGrid) {

   bool belongs = false;
   
   // loop through all neighbors (including diagonals)
   for (int ix = -1; ix <= 1; ++ix) {
      for (int iy = -1; iy <= 1; ++iy) {
         for (int iz = -1; iz <= 1; ++iz) {
            
            // not strictly necessary but logically we should not consider the cell itself
            // among its neighbors.
            if( ix == 0 && iy == 0 && iz == 0 || !technicalGrid.get(x+ix,y+iy,z+iz)) {
               continue;
            }
            
            if(layer == 1 && technicalGrid.get(x+ix,y+iy,z+iz)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               // in the first layer, boundary cell belongs if it has a non-boundary neighbor
               belongs = true;
               return belongs;
               
            } else if (layer > 1 && technicalGrid.get(x+ix,y+iy,z+iz)->sysBoundaryLayer == layer - 1) {
               // in all other layers, boundary cell belongs if it has a neighbor in the previous layer
               belongs = true;
               return belongs;
            }
         }
      }
   }

   return belongs;
}

void setupTechnicalFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells, FsGrid< fsgrids::technical, 2>& technicalGrid) {

   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);   
   technicalGrid.setupForTransferIn(nCells);

   // Setup transfer buffers
   std::vector< fsgrids::technical > transferBuffer(cells.size());
   
#pragma omp parallel for
   for(uint i = 0; i < cells.size(); ++i) {

      fsgrids::technical* thisCellData = &transferBuffer[i];
      // Data needs to be collected from some different places for this grid.
      thisCellData->sysBoundaryFlag = mpiGrid[cells[i]]->sysBoundaryFlag;
      // Remove boundary layer copy here
      // thisCellData->sysBoundaryLayer = mpiGrid[cells[i]]->sysBoundaryLayer;
      thisCellData->maxFsDt = std::numeric_limits<Real>::max();        
   }

   for(uint i = 0; i < cells.size(); ++i) {
      
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, cells[i]);
      
      for (auto fsgridId : fsgridIds) {
         // std::cout << "fsgridId: " << fsgridId << ", fsgrid Cell Coordinates:";
         // auto coords = technicalGrid.globalIDtoCellCoord(fsgridId);
         // for (auto coord : coords) std::cout << " " << coord;
         // std::cout << std::endl;
         technicalGrid.transferDataIn(fsgridId,&transferBuffer[i]);
      }
   }

   technicalGrid.finishTransfersIn();

   auto localSize = technicalGrid.getLocalSize();
   
   // Add layer calculation here. Include diagonals +-1.

   // Initialize boundary layer flags to 0.
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            technicalGrid.get(x,y,z)->sysBoundaryLayer = 0;
         }
      }
   }   

   // In dccrg initialization the max number of boundary layers is set to 3.
   const int MAX_NUMBER_OF_BOUNDARY_LAYERS = 3 * (mpiGrid.get_maximum_refinement_level() + 1);

   // loop through max number of layers
   for(uint layer = 1; layer <= MAX_NUMBER_OF_BOUNDARY_LAYERS; ++layer) {
      
      // loop through all cells in grid
#pragma omp parallel for collapse(3)
      for (int x = 0; x < localSize[0]; ++x) {
         for (int y = 0; y < localSize[1]; ++y) {
            for (int z = 0; z < localSize[2]; ++z) {
               
               // for the first layer, consider all cells that belong to a boundary, for other layers
               // consider all cells that have not yet been labeled.
               if((layer == 1 && technicalGrid.get(x,y,z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ||
                  (layer > 1 && technicalGrid.get(x,y,z)->sysBoundaryLayer == 0)) {
                  
                  if (belongsToLayer(layer, x, y, z, technicalGrid)) {
                     
                     technicalGrid.get(x,y,z)->sysBoundaryLayer = layer;
                     
                     if (layer > 1) {
                        technicalGrid.get(x,y,z)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE;
                     }
                  }
               }
            }
         }
      }
   }

   // for (int x = 0; x < localSize[0]; ++x) {
   //    for (int y = 0; y < localSize[1]; ++y) {
   //       for (int z = 0; z < localSize[2]; ++z) {
   //          std::cout << "boundary layer+flag at " << x << ", " << y << ", " << z << " = ";
   //          std::cout << technicalGrid.get(x,y,z)->sysBoundaryLayer;
   //          std::cout << " ";
   //          std::cout << technicalGrid.get(x,y,z)->sysBoundaryFlag;
   //       }
   //    }
   // }     
   //abort();
}

void getFsGridMaxDt(FsGrid< fsgrids::technical, 2>& technicalGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells) {

   int nCells = getNumberOfCellsOnMaxRefLvl(mpiGrid, cells);   
   technicalGrid.setupForTransferOut(nCells);

   // Buffer to store contents of the grid
   std::vector<fsgrids::technical> transferBuffer(nCells);
   std::vector<fsgrids::technical*> transferBufferPointer;

   int k = 0;
   for(int i=0; i< cells.size(); i++) {

      transferBufferPointer.push_back(&transferBuffer[k]);
      const auto fsgridIds = mapDccrgIdToFsGridGlobalID(mpiGrid, cells[i]);
      for (auto fsgridId : fsgridIds) {
         fsgrids::technical* thisCellData = &transferBuffer[k++];
         technicalGrid.transferDataOut(fsgridId, thisCellData);
      }
   }

   technicalGrid.finishTransfersOut();

   // After the transfer is completed, stuff the recieved maxFDt into the cells.
   #pragma omp parallel for
   for(int i=0; i< cells.size(); i++) {
      
      int dccrgId = cells[i];
      auto cellParams = mpiGrid[dccrgId]->get_cell_parameters();
      
      // Calculate the number of fsgrid cells we need to average into the current dccrg cell
      int nCells = pow(pow(2,mpiGrid.get_maximum_refinement_level() - mpiGrid.get_refinement_level(dccrgId)),3);

      cellParams[CellParams::MAXFDT] = std::numeric_limits<Real>::max();
      //cellParams[CellParams::FSGRID_RANK] = 0;
      //cellParams[CellParams::FSGRID_BOUNDARYTYPE] = 0;

      for (int iCell = 0; iCell < nCells; ++iCell) {
                 
         fsgrids::technical* thisCellData = transferBufferPointer[i] + iCell;

         if (thisCellData->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || thisCellData->sysBoundaryLayer == 1) {
         
            cellParams[CellParams::MAXFDT] = std::min(cellParams[CellParams::MAXFDT],thisCellData->maxFsDt);
            
         }
         
         //TODO: Implement something for FSGRID_RANK and FSGRID_BOUNDARYTYPE
         //cellParams[CellParams::FSGRID_RANK] = thisCellData->fsGridRank;
         //cellParams[CellParams::FSGRID_BOUNDARYTYPE] = thisCellData->sysBoundaryFlag;
      }
   }
}

/*
Map from dccrg cell id to fsgrid global cell ids when they aren't identical (ie. when dccrg has refinement).
*/

std::vector<CellID> mapDccrgIdToFsGridGlobalID(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
					       CellID dccrgID) {

   const auto maxRefLvl  = mpiGrid.get_maximum_refinement_level();
   const auto refLvl = mpiGrid.get_refinement_level(dccrgID);
   const auto cellLength = pow(2,maxRefLvl-refLvl);
   const auto topLeftIndices = mpiGrid.mapping.get_indices(dccrgID);
   std::array<int,3> indices;
   std::vector<std::array<int,3>> allIndices;

   std::array<int,3> fsgridDims;
   fsgridDims[0] = P::xcells_ini * pow(2,mpiGrid.get_maximum_refinement_level());
   fsgridDims[1] = P::ycells_ini * pow(2,mpiGrid.get_maximum_refinement_level());
   fsgridDims[2] = P::zcells_ini * pow(2,mpiGrid.get_maximum_refinement_level());
   
   for (uint k = 0; k < cellLength; ++k) {
      for (uint j = 0; j < cellLength; ++j) {
         for (uint i = 0; i < cellLength; ++i) {
            indices[0] = topLeftIndices[0] + i;
            indices[1] = topLeftIndices[1] + j;
            indices[2] = topLeftIndices[2] + k;
            allIndices.push_back(indices);
         }
      }
   }

   std::vector<CellID> fsgridIDs;  


   for (auto cellCoord: allIndices) {
     
     fsgridIDs.push_back(cellCoord[0] 
			 + cellCoord[1] * fsgridDims[0] 
			 + cellCoord[2] * fsgridDims[1] * fsgridDims[0]);

   }

   return fsgridIDs;
}

