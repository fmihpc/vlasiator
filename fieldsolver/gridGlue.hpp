#pragma once
#include "../definitions.h"
#include <fsgrid.hpp>
#include <vector>
#include <array>
#include <map>
#include <set>

// Datastructure for coupling
extern std::map<int, std::set<CellID> > onDccrgMapRemoteProcessGlobal; 
extern std::map<int, std::set<CellID> > onFsgridMapRemoteProcessGlobal; 
extern std::map<CellID, std::vector<int64_t> >  onFsgridMapCellsGlobal;

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
   dVxdx,
   dVxdy,
   dVxdz,
   dVydx,
   dVydy,
   dVydz,
   dVzdx,
   dVzdy,
   dVzdz,
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
void getFieldsFromFsGrid(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volumeFieldsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
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
			FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid);

/*Compute coupling DCCRG <=> FSGRID 

  onDccrgMapRemoteProcess   maps fsgrid processes (key) => set of dccrg cellIDs owned by current rank that map to  the fsgrid cells owned by fsgrid process (val)

  onFsgridMapRemoteProcess  maps dccrg processes  (key) => set of dccrg cellIDs owned by dccrg-process that map to current rank fsgrid cells 
  onFsgridMapCells          maps remote dccrg CellIDs to local fsgrid cells
*/

// this function is declared here as it is a template function

template <typename T, int stencil> void computeCoupling(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
               const std::vector<CellID>& cells,
               FsGrid< T, stencil>& momentsGrid) {
    
   phiprof::Timer couplingTimerActual {"CouplingTimerActual"};

   //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
   std::vector<CellID> dccrgCells = cells;
   std::sort(dccrgCells.begin(), dccrgCells.end());

   //make sure the datastructures are clean
   onDccrgMapRemoteProcessGlobal.clear();
   onFsgridMapRemoteProcessGlobal.clear();
   onFsgridMapCellsGlobal.clear();
  
  
   //size of fsgrid local part
   const std::array<FsGridTools::FsIndex_t, 3> gridDims(momentsGrid.getLocalSize());
  
 
   //Compute what we will receive, and where it should be stored
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               const std::array<FsGridTools::FsIndex_t, 3> globalIndices = momentsGrid.getGlobalIndices(i,j,k);
               const dccrg::Types<3>::indices_t  indices = {{(uint64_t)globalIndices[0],
                        (uint64_t)globalIndices[1],
                        (uint64_t)globalIndices[2]}}; //cast to avoid warnings
         CellID dccrgCell = mpiGrid.get_existing_cell(indices, 0, mpiGrid.mapping.get_maximum_refinement_level());
        
         int process = mpiGrid.get_process(dccrgCell);
         FsGridTools::LocalID fsgridLid = momentsGrid.LocalIDForCoords(i,j,k);
         //int64_t  fsgridGid = momentsGrid.GlobalIDForCoords(i,j,k);
         onFsgridMapRemoteProcessGlobal[process].insert(dccrgCell); //cells are ordered (sorted) in set
         onFsgridMapCellsGlobal[dccrgCell].push_back(fsgridLid);
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
         onDccrgMapRemoteProcessGlobal[process].insert(dccrgCells[i]); //add to map
      }    
   }
}
