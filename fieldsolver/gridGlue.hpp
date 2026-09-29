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
 * \param moments fsgrid for moments quantities
 * \param technical fsgrid with technical parameters
 * \param fsgrid fsgrid container
 * \param dt2 Whether to copy base moments, or _DT2 moments
 *
 * This function assumes that proper grid coupling has been set up.
 */
void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& moments,
                           fsgrids::technicalspan technical, FieldSolverGrid &fsgrid, bool dt2 = false);

/*! Copy field solver result (VOLB, VOLE, VOLPERB derivatives, gradpe) and store them back into DCCRG
 * \param volumefields fsgrid for volume-averaged fields
 * \param bgb fsgrid for background magnetic field
 * \param egradpe fsgrid for grad(Pe) component of electric field
 * \param dmoments fsgrid for moments derivatives
 * \param technical fsgrid with technical parameters
 * \param fsgrid fsgrids container
 * \param mpiGrid The DCCRG grid
 * \param cells List of local DCCRG cells
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getFieldsFromFsGrid(fsgrids::constvolspan volumefields,
                         fsgrids::constbgbspan bgb,
                         fsgrids::constegradpespan egradpe,
                         fsgrids::constdmomentsspan dmoments,
                         fsgrids::consttechnicalspan technical, FieldSolverGrid& fsgrid,
                         dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                         const std::vector<CellID>& cells);

/*! Copy background B fields and store them into DCCRG
 * \param bgb fsgrid for background magnetic field
 * \param technical fsgrid with technical parameters
 * \param fsgrid fsgrids container 
 * \param mpiGrid The DCCRG grid
 * \param cells List of local DCCRG cells
 *
 * This function assumes that proper grid coupling has been set up.
 */
void getBgFieldsAndDerivativesFromFsGrid(fsgrid::FsData<std::array<Real, fsgrids::bgbfield::N_BGB>>& bgb,
                                         fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                                         dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                         const std::vector<CellID>& cells);

int getNumberOfCellsOnMaxRefLvl(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                const std::vector<CellID>& cells);

void feedBoundaryIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			const std::vector<CellID>& cells,
			fsgrids::technicalspan technical, FieldSolverGrid &fsgrid);

/*Compute coupling DCCRG <=> FSGRID 

  onDccrgMapRemoteProcess   maps fsgrid processes (key) => set of dccrg cellIDs owned by current rank that map to  the fsgrid cells owned by fsgrid process (val)

  onFsgridMapRemoteProcess  maps dccrg processes  (key) => set of dccrg cellIDs owned by dccrg-process that map to current rank fsgrid cells 
  onFsgridMapCells          maps remote dccrg CellIDs to local fsgrid cells
*/

// this function is declared here as it is a template function

template <int STENCIL>
void computeCoupling(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& cells,
                     fsgrid::FsGrid<STENCIL>& fsgrid, fsgrids::technicalspan technical) {

   phiprof::Timer couplingTimerActual {"CouplingTimerActual"};

   //sorted list of dccrg cells. cells is typicall already sorted, but just to make sure....
   std::vector<CellID> dccrgCells = cells;
   std::sort(dccrgCells.begin(), dccrgCells.end());

   //make sure the datastructures are clean
   onDccrgMapRemoteProcessGlobal.clear();
   onFsgridMapRemoteProcessGlobal.clear();
   onFsgridMapCellsGlobal.clear();
  
   const auto maxRefLevel = mpiGrid.mapping.get_maximum_refinement_level();

   //Compute what we will receive, and where it should be stored
   //Probably shouldn't be parallelised unless doing the inserts atomically
   fsgrid.serial_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                     phiprof::initializeTimer("Coupling loop"), technical,
                     [=, &mpiGrid](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
      const std::array<fsgrid::FsSize_t, 3> globalIndices = coordinates.localToGlobal(stencil.i, stencil.j, stencil.k);
      const dccrg::Types<3>::indices_t  indices = {
         {(uint64_t)globalIndices[0], (uint64_t)globalIndices[1], (uint64_t)globalIndices[2]}}; //cast to avoid warnings
      const CellID dccrgCell =
         mpiGrid.get_existing_cell(indices, 0, maxRefLevel);

      const int process = mpiGrid.get_process(dccrgCell);
      const fsgrid::LocalID fsgridLid = coordinates.localIDFromLocalCoordinates(stencil.i, stencil.j, stencil.k);
      onFsgridMapRemoteProcessGlobal[process].insert(dccrgCell); // cells are ordered (sorted) in set
      onFsgridMapCellsGlobal[dccrgCell].push_back(fsgridLid);
   });

   // Compute where to send data and what to send
   for(uint64_t i=0; i< dccrgCells.size(); i++) {
      //compute to which processes this cell maps
      const std::vector<CellID> fsCells = mapDccrgIdToFsGridGlobalID(mpiGrid, dccrgCells[i]);

      //loop over fsgrid cells which this dccrg cell maps to
      for (auto const &fsCellID : fsCells) {
         const int process = fsgrid.getTaskForGlobalID(fsCellID);      // process on fsgrid
         onDccrgMapRemoteProcessGlobal[process].insert(dccrgCells[i]); //add to map
      }    
   }
}
