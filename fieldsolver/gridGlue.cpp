#include <dccrg.hpp>
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
      std::array<Real, fsgrids::moments::N_MOMENTS>* cellData = 
      momentsGrid.transferDataIn(i,mpiGrid[i]->get_cell_parameters[CellParams::RHO]);
   }

   momentsGrid.finishTransfersIn();
}
