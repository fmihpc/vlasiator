// use DCCRG version Nov 8th 2018 01482cfba8
#include "../grid.h"
using namespace std;
using namespace spatial_cell;

#include "cpu_trans_pencils.hpp"

// Cell lists for local translation
std::unordered_set<CellID> LocalSet_x;
std::unordered_set<CellID> LocalSet_y;
std::unordered_set<CellID> LocalSet_z;

std::array<setOfPencils,3> DimensionPencils;
std::array<std::unordered_set<CellID>,3> DimensionTargetCells;

//Is cell translated? It is not translated if DO_NO_COMPUTE or if it is sysboundary cell and not in first sysboundarylayer
bool do_translate_cell(SpatialCell* SC){
   if(SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
      (SC->sysBoundaryLayer != 1 && SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
      return false;
   else
      return true;
}

bool check_is_translated(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID cid, int dimension) {
   if (P::vlasovSolverGhostTranslate) {
      switch (dimension) {
         case 0:
            if (LocalSet_x.count(cid)) {
               return true;
            }
            break;
         case 1:
            if (LocalSet_y.count(cid)) {
               return true;
            }
            break;
         case 2:
            if (LocalSet_z.count(cid)) {
               return true;
            }
            break;
         default:
            cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
            abort();
      }
      return false;
   } else {
      if (mpiGrid.is_local(cid)) return true;
      return false;
   }
}

bool check_is_written_to(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID cid, int dimension) {
   if (P::vlasovSolverGhostTranslate) {
      SpatialCell *cell = mpiGrid[cid];
      // Order is z -> x -> y
      switch (dimension) {
         // checks for (cell) && (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) are before call
         case 0: // Second direction (x): Write into all cells which are used in y-translation
            if (LocalSet_y.count(cid)) {
               return true;
            }
            break;
         case 1: // Last direction (y): Write only into local cells
            if (mpiGrid.is_local(cid)) {
               return true;
            }
            break;
         case 2: // First direction (z): Write into all cells which are used in x-translation
            if (LocalSet_x.count(cid)) {
               return true;
            }
            break;
         default:
            cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
            abort();
      }
      return false;
   } else {
      if (mpiGrid.is_local(cid)) return true;
      return false;
   }
}

/* Get the one-dimensional neighborhood index for a given direction and neighborhood size.
 *
 * @param dimension spatial dimension of neighborhood
 * @param stencil neighborhood size in cells
 * @return neighborhood index that can be passed to DCCRG functions
 */
int getNeighborhood(const uint dimension, const uint stencil) {

   int neighborhood = 0;

   if (stencil == 1) {
      switch (dimension) {
      case 0:
         neighborhood = VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID;
         break;
      default:
         cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
         abort();
      }
   }
   if (stencil > 1) {
      switch (dimension) {
      case 0:
         neighborhood = VLASOV_SOLVER_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = VLASOV_SOLVER_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = VLASOV_SOLVER_Z_NEIGHBORHOOD_ID;
         break;
      default:
         cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
         abort();
      }
   }
   return neighborhood;
}

/** 
    Helper function for locating unique cells in a given direction
*/
void findNeighborhoodCells(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const CellID startingCellID,
                           uint dimension,
                           uint searchLength,
                           std::vector<CellID>& foundCells) {

   int neighborhood = getNeighborhood(dimension,searchLength);
   foundCells.clear();

   SpatialCell *ccell = mpiGrid[startingCellID];
   if (!ccell) return;

   std::set< int > distancesplus;
   std::set< int > distancesminus;
   std::unordered_set<CellID> foundNeighbors;
   std::unordered_set<CellID> foundSet;

   const auto* NbrPairs = mpiGrid.get_neighbors_of(startingCellID, neighborhood);

   // Create list of unique distances
   // Check for already found cells is required as one cell can be listed at several different distances
   for (const auto& nbrPair : *NbrPairs) {
      if(nbrPair.second[dimension] > 0) {
         if (foundNeighbors.find(nbrPair.first) == foundNeighbors.end()) {
            distancesplus.insert(nbrPair.second[dimension]);
            foundNeighbors.insert(nbrPair.first);
         }
      }
      if(nbrPair.second[dimension] < 0) {
         if (foundNeighbors.find(nbrPair.first) == foundNeighbors.end()) {
            distancesminus.insert(-nbrPair.second[dimension]);
            foundNeighbors.insert(nbrPair.first);
         }
      }
   }

   int iSrc = searchLength-1;
   // Iterate through positive distances starting from the smallest distance.
   for (auto it = distancesplus.begin(); it != distancesplus.end(); ++it) {
      if (iSrc < 0) break; // found enough elements
      // Check all neighbors at distance *it
      for (const auto& nbrPair : *NbrPairs) {
         SpatialCell *ncell = mpiGrid[nbrPair.first];
         if (!ncell) continue;
         int distanceInRefinedCells = nbrPair.second[dimension];
         if (distanceInRefinedCells == *it) {
            foundSet.insert(nbrPair.first);
         }
      } // end loop over neighbors
      iSrc--;
   } // end loop over positive distances

   iSrc = searchLength-1;
   // Iterate through negtive distances starting from the smallest distance.
   for (auto it = distancesminus.begin(); it != distancesminus.end(); ++it) {
      if (iSrc < 0) break; // found enough elements
      // Check all neighbors at distance *it
      for (const auto& nbrPair : *NbrPairs) {
         SpatialCell *ncell = mpiGrid[nbrPair.first];
         if (!ncell) continue;
         int distanceInRefinedCells = -nbrPair.second[dimension];
         if (distanceInRefinedCells == *it) {
            foundSet.insert(nbrPair.first);
         }
      } // end loop over neighbors
      iSrc--;
   } // end loop over negative distances
   foundCells.assign(foundSet.begin(), foundSet.end());
}

void prepareGhostTranslationCellLists(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const vector<CellID>& localPropagatedCells) {
   // return if there's no cells to start with
   if(localPropagatedCells.size() == 0) {
      std::cerr<<"No cells!"<<std::endl;
      return;
   }

   // Clear existing lists
   LocalSet_x.clear();
   LocalSet_y.clear();
   LocalSet_z.clear();
   std::vector<CellID> foundCells;

   // Cell lists include:
   // 1) Any local (translated) non-sysboundary cells
   // 2) per-dimension, in translation order, one layer of face neighbours, remote or local (including translated sysboundary cells)
   // Done only at LB so not threaded for now

   // do_translate_cell is in cpu_trans_map.cpp:
   // It is not translated if DO_NO_COMPUTE or if it is sysboundary cell and not in first sysboundarylayer

   /** Translation order (dimensions) is 1: z 2: x 3: y
       Prepare in reverse order
       First: y-direction
   */
   for (uint i=0; i<localPropagatedCells.size(); i++) {
      CellID c = localPropagatedCells[i];
      SpatialCell *ccell = mpiGrid[c];
      if (!ccell) continue;

      // Is the cell translated?
      if (!do_translate_cell(ccell)) continue;

      // Is the cell a non-sysboundary cell?
      // (translated sysboundarycells are only included via the neighbourhood loop)
      if (ccell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;

      LocalSet_y.insert(c);
      ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_Y] = true;
      // y-translation
      findNeighborhoodCells(mpiGrid, c, 1, 1, foundCells);
      for (uint j=0; j<foundCells.size(); j++) {
         CellID n = foundCells[j];
         SpatialCell *ncell = mpiGrid[n];
         if (!ncell) continue;
         // Is the cell translated?
         if (!do_translate_cell(ncell)) continue;
         LocalSet_y.insert(n);
         ncell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_Y] = true;
      }
   }

   /** Now use y-translation source cells as starting points
       and evaluate x-direction
    */
   for (auto it=LocalSet_y.begin(); it!=LocalSet_y.end(); ++it) {
      CellID c = (*it);
      SpatialCell *ccell = mpiGrid[c];
      if (!ccell) continue;

      // Is the cell translated?
      if (!do_translate_cell(ccell)) continue;

      // Is the cell a non-sysboundary cell?
      // (translated sysboundarycells are only included via the neighbourhood loop)
      if (ccell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;

      LocalSet_x.insert(c);
      ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X] = true;
      // x-translation
      findNeighborhoodCells(mpiGrid, c, 0, 1, foundCells);
      for (uint j=0; j<foundCells.size(); j++) {
         CellID n = foundCells[j];
         SpatialCell *ncell = mpiGrid[n];
         if (!ncell) continue;
         // Is the cell translated?
         if (!do_translate_cell(ncell)) continue;
         LocalSet_x.insert(n);
         ncell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X] = true;
      }
   }


   /** Now use x-translation source cells as starting points
       and evaluate the (last) z-direction
    */
   for (auto it=LocalSet_x.begin(); it!=LocalSet_x.end(); ++it) {
      CellID c = (*it);
      SpatialCell *ccell = mpiGrid[c];
      if (!ccell) continue;

      // Is the cell translated?
      if (!do_translate_cell(ccell)) continue;

      // Is the cell a non-sysboundary cell?
      // (translated sysboundarycells are only included via the neighbourhood loop)
      if (ccell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;

      LocalSet_z.insert(c);
      ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_Z] = true;

      // z-translation
      findNeighborhoodCells(mpiGrid, c, 2, 1, foundCells);
      for (uint j=0; j<foundCells.size(); j++) {
         CellID n = foundCells[j];
         SpatialCell *ncell = mpiGrid[n];
         if (!ncell) continue;
         // Is the cell translated?
         if (!do_translate_cell(ncell)) continue;
         LocalSet_z.insert(n);
         ncell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_Z] = true;
      }
   }
   return;
}


void flagSpatialCellsForAmrCommunication(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                         const vector<CellID>& localPropagatedCells) {

   // Only flag/unflag cells if AMR is active
   if (mpiGrid.get_maximum_refinement_level()==0) return;

   // return if there's no cells to flag
   if(localPropagatedCells.size() == 0) {
      return;
   }

   for (int dimension=0; dimension<3; dimension++) {
      // These neighborhoods now include the AMR addition beyond the regular vlasov stencil
      int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

      // Set flags: loop over local cells
#pragma omp parallel for
      for (uint i=0; i<localPropagatedCells.size(); i++) {
         CellID c = localPropagatedCells[i];
         SpatialCell *ccell = mpiGrid[c];
         if (!ccell) continue;

         // Translated cells also need to be included in order to communicate boundary cell VDFs.
         // Attempting to leave these out for the x or y dimensions also resulted in diffs.
         // if (!do_translate_cell(ccell)) continue;

         // Start with false
         ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] = false;

         // In dimension, check iteratively if any neighbors up to VLASOV_STENCIL_WIDTH distance away are on a different process
         const auto* NbrPairs = mpiGrid.get_neighbors_of(c, neighborhood);

         // Create list of unique distances
         std::set< int > distancesplus;
         std::set< int > distancesminus;
         std::set<CellID> foundNeighborsP;
         std::set<CellID> foundNeighborsM;
         /** Using sets of cells as well, we should only get one distance per
             (potentially less refined) cell. This should result in safe behaviour
             as long as the neighborhood of a cell does not contain cells with a
             refinement level more than 1 level apart from the cell itself.
         */
         for (const auto& nbrPair : *NbrPairs) {
            if(nbrPair.second[dimension] > 0) {
               if (foundNeighborsP.find(nbrPair.first) == foundNeighborsP.end()) {
                  distancesplus.insert(nbrPair.second[dimension]);
                  foundNeighborsP.insert(nbrPair.first);
               }
            }
            if(nbrPair.second[dimension] < 0) {
               if (foundNeighborsM.find(nbrPair.first) == foundNeighborsM.end()) {
                  distancesminus.insert(-nbrPair.second[dimension]);
                  foundNeighborsM.insert(nbrPair.first);
               }
            }
         }

         foundNeighborsP.clear();
         foundNeighborsM.clear();

         int iSrc = VLASOV_STENCIL_WIDTH - 1;
         // Iterate through positive distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
         for (auto it = distancesplus.begin(); it != distancesplus.end(); ++it) {
            if (ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] == true) iSrc = -1;
            if (iSrc < 0) break; // found enough elements
            // Check all neighbors at distance *it
            for (const auto& nbrPair : *NbrPairs) {
               SpatialCell *ncell = mpiGrid[nbrPair.first];
               if (!ncell) continue;
               int distanceInRefinedCells = nbrPair.second[dimension];
               if (distanceInRefinedCells == *it) {
                  if (foundNeighborsP.find(nbrPair.first) != foundNeighborsP.end()) continue;
                  foundNeighborsP.insert(nbrPair.first);
                  if (!mpiGrid.is_local(nbrPair.first)) {
                     ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] = true;
                     break;
                  }
               }
            } // end loop over neighbors
            iSrc--;
         } // end loop over positive distances

         iSrc = VLASOV_STENCIL_WIDTH - 1;
         // Iterate through negtive distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
         for (auto it = distancesminus.begin(); it != distancesminus.end(); ++it) {
            if (ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] == true) iSrc = -1;
            if (iSrc < 0) break; // found enough elements
            // Check all neighbors at distance *it
            for (const auto& nbrPair : *NbrPairs) {
               SpatialCell *ncell = mpiGrid[nbrPair.first];
               if (!ncell) continue;
               int distanceInRefinedCells = -nbrPair.second[dimension];
               if (distanceInRefinedCells == *it) {
                  if (foundNeighborsM.find(nbrPair.first) != foundNeighborsM.end()) continue;
                  foundNeighborsM.insert(nbrPair.first);
                  if (!mpiGrid.is_local(nbrPair.first)) {
                     ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] = true;
                     break;
                  }
               }
            } // end loop over neighbors
            iSrc--;
         } // end loop over negative distances
      } // end loop over local propagated cells
   } // end loop over dimensions
   return;
}

/* Get cellIDs for spatial cells that are considered target / source cells for a pencil.
 *
 * Source cells are cells that the pencil reads data from to compute polynomial
 * fits that are used for propagation in the vlasov solver. All cells included
 * in the pencil proper + VLASOV_STENCIL_WIDTH cells on both ends are source cells.
 * Invalid cells are replaced by closest good cells.
 * Boundary cells are included.
 *
 * Target cells are cells that the pencil writes data into after translation by
 * the vlasov solver. All cells included in the pencil proper + 1 cells on both ends
 * are target cells.
 *
 * There is only one list of cellIDs for the pencil, where each source cell has also
 * their width stored, and for target cells, the relative contribution is stored. If
 * the cell in question is not a target cell, the contribution is set to zero.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [inout] ids pointer to subsection of vector of ids for actual pencil
 * @param [in] L length of pencil (including stencil cells)
 * @param [in] dimension spatial dimension
 * @param [in] path index of the desired face neighbor when going to a higher refinement level
 * @param [out] source pointer to subsection of vector storing cell widths
 * @param [out] targetRatios pointer to subsection of vector storing relative contribution of pencil to target cells
 */
void computeSpatialSourceCellsForPencil(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                        CellID *ids,
                                        const uint L,
                                        const uint dimension,
                                        std::vector<uint> path,
                                        Realf* sourceDZ,
                                        Realf* targetRatios
                                        ){

   // These neighborhoods now include the AMR addition beyond the regular vlasov stencil
   int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);
   stringstream ss;
   for (uint j = 0; j < L; ++j) {
      ss<< ids[j] << " ";
   }

   // Insert pointers for neighbors of ids.front() and ids.back()
   const auto* frontNbrPairs = mpiGrid.get_neighbors_of(ids[VLASOV_STENCIL_WIDTH], neighborhood);
   const auto* backNbrPairs  = mpiGrid.get_neighbors_of(ids[L-VLASOV_STENCIL_WIDTH-1],  neighborhood);
   // Create list of unique distances in the negative direction from the first cell in pencil
   std::set< int > distances;
   for (const auto& nbrPair : *frontNbrPairs) {
      if(nbrPair.second[dimension] < 0) {
         // gather absolute distance values
         distances.insert(-nbrPair.second[dimension]);
      }
   }

   int iSrc = VLASOV_STENCIL_WIDTH - 1;
   // Iterate through distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
   for (auto it = distances.begin(); it != distances.end(); ++it) {
      if (iSrc < 0) break; // found enough elements

      // Collect all neighbors at distance *it to a vector
      std::vector< CellID > neighbors;
      for (const auto& nbrPair : *frontNbrPairs) {
         int distanceInRefinedCells = -nbrPair.second[dimension];
         if(distanceInRefinedCells == *it) neighbors.push_back(nbrPair.first);
      }
      // Get rid of duplicate neighbor cells at single distance
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());

      // Find source cells (VLASOV_STENCIL_WIDTH at each end)
      int refLvl = mpiGrid.get_refinement_level(ids[VLASOV_STENCIL_WIDTH]);
      int pathPos = 0;
      if (path.size() > refLvl) {
         pathPos = path[refLvl];
      }
      if (neighbors.size() == 1) {
         if (ids[iSrc+1] == neighbors.at(0)) continue; // already found this cell for different distance
         ids[iSrc--] = neighbors.at(0);
      } else if ( pathPos < neighbors.size() ) {
         if (ids[iSrc+1] == neighbors.at(pathPos)) continue; // already found this cell for different distance (should not happen)
         ids[iSrc--] = neighbors.at(pathPos);
      } else {
         ss<<"error too few front neighbors for path! cellid "<<ids[VLASOV_STENCIL_WIDTH]<<
            " Nsize "<<neighbors.size()<<" L "<<L<<" refLvl "<<refLvl<<" iSrc "<<iSrc<<
            " pathsize "<<path.size()<<" path "<<path[refLvl]<<std::endl;
         std::cerr<<ss.str();
      }
   }

   distances.clear();
   // Create list of unique distances in the positive direction from the last cell in pencil
   for (const auto& nbrPair : *backNbrPairs) {
      if(nbrPair.second[dimension] > 0) {
         distances.insert(nbrPair.second[dimension]);
      }
   }

   // Iterate through distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
   // Distances are positive here so smallest distance has smallest value.
   iSrc = L - VLASOV_STENCIL_WIDTH;
   for (auto it = distances.begin(); it != distances.end(); ++it) {
      if (iSrc >= (int)L) break; // Found enough cells

      // Collect all neighbors at distance *it to a vector
      std::vector< CellID > neighbors;
      for (const auto& nbrPair : *backNbrPairs) {
         int distanceInRefinedCells = nbrPair.second[dimension];
         if(distanceInRefinedCells == *it) neighbors.push_back(nbrPair.first);
      }
      // Get rid of duplicate neighbor cells at single distance
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());

      int refLvl = mpiGrid.get_refinement_level(ids[L-VLASOV_STENCIL_WIDTH-1]);
      int pathPos = 0;
      if (path.size() > refLvl) {
         pathPos = path[refLvl];
      }
      if (neighbors.size() == 1) {
         if (ids[iSrc-1] == neighbors.at(0)) continue; // already found this cell for different distance
         ids[iSrc++] = neighbors.at(0);
      } else if ( pathPos < neighbors.size() ) {
         if (ids[iSrc-1] == neighbors.at(pathPos)) continue; // already found this cell for different distance (should not happen)
         ids[iSrc++] = neighbors.at(pathPos);
      } else {
         ss<<"error too few back neighbors for path! cellid "<<ids[L-VLASOV_STENCIL_WIDTH-1]<<
            " Nsize "<<neighbors.size()<<" L "<<L<<" refLvl "<<refLvl<<" iSrc "<<iSrc<<
            " pathsize "<<path.size()<<" path "<<path[refLvl]<<std::endl;
         std::cerr<<ss.str();
      }
   }

   /*loop to negative side and replace all invalid cells with the closest good cell*/
   CellID lastGoodCell = ids[VLASOV_STENCIL_WIDTH];
   for(int i = VLASOV_STENCIL_WIDTH - 1; i >= 0 ;--i){
      bool isGood = false;
      if (ids[i]!=0) {
         if (mpiGrid[ids[i]] != NULL
             || mpiGrid[ids[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE
             || ids[i] != mpiGrid[ids[i]]->parameters[CellParams::CELLID]
             // Last check required as DCCRG may have valid cell pointers which have not been ghost-updated?
            ) {
            if (mpiGrid[ids[i]]->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
               isGood = true;
            }
         }
      }
      if (!isGood) {
         ids[i] = lastGoodCell;
      } else {
         lastGoodCell = ids[i];
      }
   }

   /*loop to positive side and replace all invalid cells with the closest good cell*/
   lastGoodCell = ids[L - VLASOV_STENCIL_WIDTH - 1];
   for(int i = (int)L - VLASOV_STENCIL_WIDTH; i < (int)L; ++i){
      bool isGood = false;
      if (ids[i]!=0) {
         if (mpiGrid[ids[i]] != NULL
             || mpiGrid[ids[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE
             || ids[i] != mpiGrid[ids[i]]->parameters[CellParams::CELLID]
             // Last check required as DCCRG may have valid cell pointers which have not been ghost-updated?
            ) {
            if (mpiGrid[ids[i]]->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
               isGood = true;
            }
         }
      }
      if (!isGood) {
         ids[i] = lastGoodCell;
      } else {
         lastGoodCell = ids[i];
      }
   }

   // Loop over all cells and store widths in translation direction
   for (int i = 0; i < (int)L; ++i) {
      sourceDZ[i] = mpiGrid[ids[i]]->parameters[CellParams::DX+dimension];
   }

   // Loop over all cells and store pencil-to-cell cross-sectional area for target cells
   for (int i = 0; i < (int)L; ++i) {
      if ((i < VLASOV_STENCIL_WIDTH-1) || (i > (int)L-VLASOV_STENCIL_WIDTH)) {
         // Source cell, not a target cell
         targetRatios[i]=0.0;
         continue;
      }
      if (P::vlasovSolverGhostTranslate) {
         if (!check_is_written_to(mpiGrid, ids[i], dimension)) {
            targetRatios[i]=0.0;
            continue;
         }
         if (mpiGrid[ids[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            std::cerr<<"sysboundary targetcell!"<<std::endl;
            abort();
         }
      }

      if (ids[i]) {
         SpatialCell* tc = mpiGrid[ids[i]];
         if (tc && tc->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            // areaRatio is the ratio of the cross-section of the spatial cell to the cross-section of the pencil.
            const int diff = tc->SpatialCell::parameters[CellParams::REFINEMENT_LEVEL] - path.size();
            if(diff>0) {
               // Undefined behaviour! Cell is smaller than pencil (higher reflevel than path size)
               std::cerr<<"Error in path size to cell size: __FILE__:__LINE__"<<std::endl;
               targetRatios[i] = 0.0;
            } else {
               const int ratio = 1 << -diff;
               targetRatios[i] = 1.0 / (ratio*ratio);
            }
         } else { // Don't write to this cell
            targetRatios[i] = 0.0;
         }
      } else { // Don't write to this cell
         std::cerr<<"Found zero id in pencils!"<<std::endl;
         targetRatios[i] = 0.0;
      }
   }
}

/* Select one nearest neighbor of a cell on the + side in a given dimension. If the neighbor
 * has a higher level of refinement, a path variable is needed to make the selection.
 * Returns INVALID_CELLID if the nearest neighbor is not local to this process.
 * (or, if activated, not included in the translation list for local translation)
 *
 * @param grid DCCRG grid object
 * @param id DCCRG cell id
 * @param dimension spatial dimension
 * @param path index of the desired face neighbor
 * @return neighbor DCCRG cell id of the neighbor
 */
CellID selectNeighbor(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid,
                      CellID id, int dimension = 0, uint path = 0) {
   // If face neighbours are at a higher refinement level, only returns the one which
   // which has a neighbor index matching the input path

   //int neighborhood = getNeighborhood(dimension,1);
   //const auto* nbrPairs = grid.get_neighbors_of(id, neighborhood);

   vector < CellID > myNeighbors;
   CellID neighbor = INVALID_CELLID;

   // Iterate through neighbor ids in the positive direction of the chosen dimension,
   // select the neighbor indicated by path, if it is local to this process.
   for (const auto& [neighbor, dir] : grid.find_face_neighbors_of(id)) {
     if (dir == ((int)dimension + 1)) {
	 myNeighbors.push_back(neighbor);
      }
   }
   // TODO Verify: are neighbours always in the same order? Let's sort based
   // on CellID to be sure.
   std::sort(myNeighbors.begin(), myNeighbors.end());

   if( myNeighbors.size() == 0 ) {
      return neighbor; // == INVALID_CELLID
   }

   int neighborIndex = 0;
   if (myNeighbors.size() > 1) {
      neighborIndex = path;
   }

   if (check_is_translated(grid, myNeighbors[neighborIndex], dimension)) {
     neighbor = myNeighbors[neighborIndex];
   }

   return neighbor;
}

/* Recursive function for building one-dimensional pencils to cover local DCCRG cells.
 * Starts from a given seedID and proceeds finding the nearest neighbor in the given dimension
 * and adding it to the pencil until no neighbors are found or an endId is met. When a higher
 * refinement level (ie. multiple nearest neighbors) is met, the pencil splits into four
 * copies to remain at a width of 1 cell. This is done by the function calling itself recursively
 * and passing as inputs the cells added so far. The cell selected by each copy of the function
 * at a split is stored in the path variable, the same path has to be followed if a refinement
 * level is encoutered multiple times.
 *
 * @param [in] grid DCCRG grid object
 * @param [out] pencils Pencil data struct
 * @param [in] seedId DCCRG cell id where we start building the pencil.
 *             The pencil will continue in the + direction in the given dimension until an end condition is met
 * @param [in] dimension Spatial dimension
 * @param [in] path Integer value that determines which neighbor is added to the pencil when a higher refinement level is met
 * @param [in] endIds Prescribed end conditions for the pencil. If any of these cell ids is about to be added to the pencil,
 *             the builder terminates.
 */
void buildPencilsWithNeighbors( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid,
					setOfPencils &pencils, const CellID seedId,
					vector<CellID> ids, const uint dimension,
					vector<uint> path, const vector<CellID> &endIds) {

   const bool debug = false;
   CellID nextNeighbor;
   CellID id = seedId;
   int startingRefLvl = grid.get_refinement_level(id);
   bool periodic = false;
   // If this is a new pencil (instead of being a result of a pencil being split
   if( ids.size() == 0 ) {
      ids.push_back(seedId);
   }
   // If the cell where we start is refined, we need to figure out which path
   // to follow in future refined cells. This is a bit hacky but we have to
   // use the order or the children of the parent cell to figure out which
   // corner we are in.

   std::array<double, 3> coordinates = grid.get_center(seedId);
   int startingPathSize = path.size();

   // Find the "pre-existing" path for new pencils starting at higher reflevels
   if( startingRefLvl > startingPathSize ) {
      CellID myId = seedId;
      for ( int i = path.size(); i < startingRefLvl; ++i) {
         //CellID parentId = grid.mapping.get_parent(myId);
         CellID parentId = grid.get_parent(myId);

         auto myCoords = grid.get_center(myId);
         auto parentCoords = grid.get_center(parentId);
         int ix=0, iy=0;
         switch(dimension) {
         case 0:
            ix = 1;
            iy = 2;
            break;
         case 1:
            ix = 0;
            iy = 2;
            break;
         case 2:
            ix = 0;
            iy = 1;
            break;
         default:
            cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
            abort();
         }
         //int ix = (dimension + 1) % 3; // incorrect for DCCRG
         //int iy = (dimension + 2) % 3;

         int step = -1;

         if        (myCoords[ix] < parentCoords[ix] && myCoords[iy] < parentCoords[iy]) {
            step = 0;
         } else if (myCoords[ix] > parentCoords[ix] && myCoords[iy] < parentCoords[iy]) {
            step = 1;
         } else if (myCoords[ix] < parentCoords[ix] && myCoords[iy] > parentCoords[iy]) {
            step = 2;
         } else if (myCoords[ix] > parentCoords[ix] && myCoords[iy] > parentCoords[iy]) {
            step = 3;
         }

         // path needs to end up in reflevel order, whereas this loop goes in reverse order
         path.insert(path.begin(), step);
         myId = parentId;
      }
   }

   // Now start loop for gathering ids into pencil. Break when next found cell is rejected.
   while (id != INVALID_CELLID) {

      periodic = false;
      bool neighborExists = false;
      int refLvl = 0;

      // Find the refinement level in the neighboring (local) cell. Check all possible neighbors
      // in case some of them are remote.
      for (int tmpPath = 0; tmpPath < 4; ++tmpPath) {
         nextNeighbor = selectNeighbor(grid,id,dimension,tmpPath);
         if(nextNeighbor != INVALID_CELLID) {
            refLvl = max(refLvl,grid.get_refinement_level(nextNeighbor));
            neighborExists = true;
         }
      }

      // If there are no local acceptable neighbors, we can stop. This is not an error.
      if (!neighborExists) {
         break;
      }

      // Do we need to consider refinement?
      if (refLvl > 0) {
         // If we have encountered this refinement level before and stored
         // the path this builder follows, we will just take the same path
         // again.
         if ( static_cast<int>(path.size()) >= refLvl ) {

            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << " ";
               std::cout << std::endl;
            }

            nextNeighbor = selectNeighbor(grid,id,dimension,path[refLvl - 1]);
            if(nextNeighbor != INVALID_CELLID) {
               coordinates = grid.get_center(nextNeighbor);
            }
         } else {
            // We encounter a new refinement level.
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have NOT seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << ' ';
               std::cout << std::endl;
            }

            // Create a path through each neighbor cell
            for ( uint newPath : {0,1,2,3} ) {
               vector < uint > myPath = path;
               // Extend the path to cover the new reflevel
               myPath.push_back(newPath);
               nextNeighbor = selectNeighbor(grid,id,dimension,newPath);

               if ( newPath == 3 ) {
                  // This builder continues with neighbor 3 with an extended path
                  path = myPath;
                  if(nextNeighbor != INVALID_CELLID) {
                     coordinates = grid.get_center(nextNeighbor);
                  }
               } else {
                  // Spawn recursive new builders for neighbors 0,1,2
                  buildPencilsWithNeighbors(grid,pencils,id,ids,dimension,myPath,endIds);
               }
            }
         }
      } // Closes if (refLvl == 0)

      // If we found a neighbor, let's verify if it should be translated
      if(nextNeighbor != INVALID_CELLID) {
         if (debug) {
            std::cout << " Next neighbor is " << nextNeighbor << "." << std::endl;
         }
         // Non-local, non-translated, and ids belonging to other pencils are not included
         if ( std::any_of(endIds.begin(), endIds.end(), [nextNeighbor](uint i){return i == nextNeighbor;}) ||
              !do_translate_cell(grid[nextNeighbor])) {
            nextNeighbor = INVALID_CELLID;
         } else {
            // Yep, this goes in this pencil.
            ids.push_back(nextNeighbor);
         }
      }

      id = nextNeighbor;
   } // Closes while loop - end of pencil reached.

   // Get the x,y - coordinates of the pencil (in the direction perpendicular to the pencil)
   double x,y;
   int ix=0, iy=0;

   switch(dimension) {
   case 0:
      ix = 1;
      iy = 2;
      break;
   case 1:
      ix = 0;
      iy = 2;
      break;
   case 2:
      ix = 0;
      iy = 1;
      break;
   default:
      cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
      abort();
   }
   //ix = (dimension + 1) % 3; // incorrect for DCCRG
   //iy = (dimension + 2) % 3;

   x = coordinates[ix];
   y = coordinates[iy];

   pencils.addPencil(ids,x,y,periodic,path);
   return;
}

/* Determine which cells in the local DCCRG mesh should be starting points for pencils.
 * If a neighbor cell is non-local, across a periodic boundary, or in non-periodic boundary layer 1
 * then we use this cell as a seed for pencils
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] localPropagatedCells List of local cells that get propagated
 * ie. not L2-boundary or DO_NOT_COMPUTE
 * @param [in] dimension Spatial dimension
 * @param [out] seedIds list of cell ids that will be starting points for pencils
 */
void getSeedIds(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                const vector<CellID> &localPropagatedCells,
                const uint dimension,
                vector<CellID> &seedIds) {

   const bool debug = false;
   int myRank;
   if (debug) MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   // These neighborhoods now include the AMR addition beyond the regular vlasov stencil
   int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

#pragma omp parallel for
   for (uint i=0; i<localPropagatedCells.size(); i++) {
      CellID celli = localPropagatedCells[i];

      bool addToSeedIds = P::amrTransShortPencils;
      if (addToSeedIds) {
#pragma omp critical
         seedIds.push_back(celli);
         continue;
      }
      auto myIndices = mpiGrid.mapping.get_indices(celli);
      int myRefLevel;

      /* -----------------------------------------
         | A |   | B |   |_|_|_|_|   |   | C |   |
         |   |   |   |   | | | | |   |   |   |   |
         -----------------------------------------
         For optimal pencil generation, we need seedids at A, B, and C.
         A Is triggered in the first if-clause. Pencils starting from B
         will be split (but won't cause A to split), and pencils from
         C will be able to remain again un-split. These checks need to be done
         only if we aren't already at the maximum refinement level.

      */

      // First check negative face neighbors (A)
      // Returns all neighbors as (id, direction-dimension) pair pointers.
      for (const auto& [neighbor, dir] : mpiGrid.find_face_neighbors_of(celli) ) {
         if ( dir == -((int)dimension + 1) ) {
            // Check that the neighbor is not across a periodic boundary by calculating
            // the distance in indices between this cell and its neighbor.
            auto nbrIndices = mpiGrid.mapping.get_indices(neighbor);

            // If a neighbor is non-local (or not ghost-translated), across a periodic boundary,
            // or in non-periodic boundary layer >1 (non-translated cell)
            // then we use the current cell as a seed for pencils
            if ( (myIndices[dimension] < nbrIndices[dimension]) ||
                 !check_is_translated(mpiGrid, neighbor, dimension) ||
                 !do_translate_cell(mpiGrid[neighbor]) )
            {
               addToSeedIds = true;
               break;
            }
         }
      } // finish check A
      if ( addToSeedIds ) {
#pragma omp critical
         seedIds.push_back(celli);
         continue;
      }
      myRefLevel = mpiGrid.get_refinement_level(celli);
      if (mpiGrid.get_maximum_refinement_level() == myRefLevel) continue;

      /* Proceed with B, checking if the next positive neighbour has the same refinement level as ccell, but the
         second neighbour a higher one. Iterate through positive distances for VLASOV_STENCIL_WIDTH elements
         starting from the smallest distance. */

      // Gather neighbours in neighbourhood stencil
      const auto* nbrPairs  = mpiGrid.get_neighbors_of(celli, neighborhood);
      // Create list of unique neighbour distances in both directions (using ordered sets)
      std::set< int > distancesplus;
      std::set< int > distancesminus;
      for (const auto& nbrPair : *nbrPairs) {
         if(nbrPair.second[dimension] > 0) {
            distancesplus.insert(nbrPair.second[dimension]);
         }
         if(nbrPair.second[dimension] < 0) {
            // gather absolute distance values for correct order
            distancesminus.insert(-nbrPair.second[dimension]);
         }
      }
      int iSrc = VLASOV_STENCIL_WIDTH-1;
      for (auto it = distancesplus.begin(); it != distancesplus.end(); ++it) {
         if (iSrc < 0) break; // found enough elements
         for (const auto& nbrPair : *nbrPairs) {
            int distanceInRefinedCells = nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               // Break search if we are not at the final entry, and have different refinement level
               if (iSrc!=0 && mpiGrid.get_refinement_level(nbrPair.first)!=myRefLevel) {
                  iSrc = -1;
                  break;
               }
               // Flag as seed id if VLASOV_STENCIL_WIDTH positive neighbour is at higher refinement level
               if (iSrc==0 && mpiGrid.get_refinement_level(nbrPair.first)>myRefLevel) {
                  addToSeedIds = true;
                  break;
               }
            }
         }
         iSrc--;
      } // Finish B check

      if ( addToSeedIds ) {
#pragma omp critical
         seedIds.push_back(celli);
         continue;
      }
      /* Proceed with C, checking if the next two negative neighbours have the same refinement level as ccell, but the
         third neighbour a higher one. Iterate through negative distances for VLASOV_STENCIL_WIDTH+1 elements
         starting from the smallest distance. */
      iSrc = VLASOV_STENCIL_WIDTH;
      for (auto it = distancesminus.begin(); it != distancesminus.end(); ++it) {
         if (iSrc < 0) break; // found enough elements
         for (const auto& nbrPair : *nbrPairs) {
            int distanceInRefinedCells = -nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               // Break search if we are not at the final entry, and have different refinement level
               if (iSrc!=0 && mpiGrid.get_refinement_level(nbrPair.first)!=myRefLevel) {
                  iSrc = -1;
                  break;
               }
               // Flag as seed id if VLASOV_STENCIL_WIDTH+1 positive neighbour is at higher refinement level
               if (iSrc==0 && mpiGrid.get_refinement_level(nbrPair.first)>myRefLevel) {
                  addToSeedIds = true;
                  break;
               }
            }
         }
         iSrc--;
      } // Finish C check

      if ( addToSeedIds ) {
#pragma omp critical
         seedIds.push_back(celli);
      }
   }

   if(debug) {
      cout << "Rank " << myRank << ", Seed ids are: ";
      for (const auto seedId : seedIds) {
         cout << seedId << " ";
      }
      cout << endl;
   }
}

/* Check whether the ghost cells around the pencil contain higher refinement than the pencil does.
 * If they do, the pencil must be split to match the finest refined ghost cell.
 *
 * @param mpiGrid DCCRG grid object
 * @param pencils Pencil data struct
 * @param dimension Spatial dimension
 */
void check_ghost_cells(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                       setOfPencils& pencils,
                       uint dimension) {

   const bool debug = false;
   int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

   int myRank;
   if(debug) {
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   }

   std::vector<CellID> pencilIdsToSplit;

#pragma omp parallel for
   for (uint pencili = 0; pencili < pencils.N; ++pencili) {

      // This check isn't in use at the moment, because no pencils are ever flagged periodic..
      if (pencils.periodic[pencili]) {
         continue;
      }

      // This returns a list of only the central cells, excluding the stencil
      auto ids = pencils.getIds(pencili);

      // It is possible that the pencil has already been refined by the pencil building algorithm
      // and is on a higher refinement level than the refinement level of any of the cells it contains
      // due to e.g. process boundaries. Example: x is the pencil, N is non-local cells.
      /* -----------------------------------------
         |   |   |   |   |_|_|_|_| N | N | N | N |
         |xxx|xxx|xxx|xxx|N|N|N|N|   |   |   |   |
         -----------------------------------------
      */
      int maxPencilRefLvl = pencils.path[pencili].size();
      int maxNbrRefLvl = 0;

      const auto* frontNeighbors = mpiGrid.get_neighbors_of(ids.front(),neighborhood);
      const auto* backNeighbors  = mpiGrid.get_neighbors_of(ids.back(),neighborhood);

      // Create list of unique distances in the negative direction from the first cell in pencil
      std::set< int > distances; // is sorted
      for (const auto& nbrPair : *frontNeighbors) {
         if(nbrPair.second[dimension] < 0) {
            // gather absolute distance values
            distances.insert(-nbrPair.second[dimension]);
         }
      }
      int foundcells = 0;
      CellID lastcell = INVALID_CELLID;
      // Iterate through distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
      for (auto it = distances.begin(); it != distances.end(); ++it) {
         for (const auto& nbrPair : *frontNeighbors) {
            if (nbrPair.first==lastcell) continue;
            int distanceInRefinedCells = -nbrPair.second[dimension];
            if (distanceInRefinedCells == *it) {
               maxNbrRefLvl = max(maxNbrRefLvl,mpiGrid.get_refinement_level(nbrPair.first));
               lastcell = nbrPair.first;
               foundcells++;
               continue;
            }
         }
         if (foundcells >= VLASOV_STENCIL_WIDTH) {
            break; // checked enough distances
         }
      }

      // Create list of unique distances in the positive direction from the last cell in pencil
      distances.clear();
      for (const auto& nbrPair : *backNeighbors) {
         if(nbrPair.second[dimension] > 0) {
            distances.insert(nbrPair.second[dimension]);
         }
      }
      foundcells = 0;
      lastcell = INVALID_CELLID;
      for (auto it = distances.begin(); it != distances.end(); ++it) {
         for (const auto& nbrPair : *backNeighbors) {
            if (nbrPair.first==lastcell) continue;
            int distanceInRefinedCells = nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               maxNbrRefLvl = max(maxNbrRefLvl,mpiGrid.get_refinement_level(nbrPair.first));
               lastcell = nbrPair.first;
               foundcells++;
               continue;
            }
         }
         if (foundcells >= VLASOV_STENCIL_WIDTH) {
            break; // checked enough distances
         }
      }

      if (maxNbrRefLvl > maxPencilRefLvl) {
         if(debug) {
            std::cout << "I am rank " << myRank << ". ";
            std::cout << "Found refinement level " << maxNbrRefLvl << " in one of the ghost cells of pencil " << pencili << ". ";
            std::cout << "Highest refinement level in this pencil is " << maxPencilRefLvl;
            std::cout << ". Splitting pencil " << pencili << endl;
         }
         // Let's avoid modifying pencils while we are looping over it. Write down the indices of pencils
         // that need to be split and split them later.
#pragma omp critical
         {
            pencilIdsToSplit.push_back(pencili);
         }
      }
   }

   // No threading here! Splitting requires knowledge of all
   // already existing pencils.
   for (auto pencili: pencilIdsToSplit) {

      Real dx = 0.0;
      Real dy = 0.0;
      auto ids = pencils.getIds(pencili);
      switch(dimension) {
      case 0:
         dx = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DY];
         dy = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DZ];
         break;
      case 1:
         dx = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DX];
         dy = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DZ];
         break;
      case 2:
         dx = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DX];
         dy = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DY];
         break;
      }

// WARNING threading inside this function
      pencils.split(pencili,dx,dy);

   }
}

/* Checks that each local spatial cell appears in pencils at least 1 time.
 *
 * @param mpiGrid DCCRG grid object
 * @param cells Local spatial cells
 * @param pencils Pencil data struct
 */
bool checkPencils(
   const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID> &cells,
   const setOfPencils& pencils
) {
   bool correct = true;
   for (auto id : cells) {
      if (mpiGrid[id]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY )  {
         int myCount = std::count(pencils.ids.begin(), pencils.ids.end(), id);
         if( myCount == 0) {
            std::cerr << "ERROR: Cell ID " << id << " Appears in pencils " << myCount << " times!"<< std::endl;
            correct = false;
         }
      }
   }
   for (uint ipencil = 0; ipencil < pencils.N; ++ipencil) {
      cint nPencilsThroughThisCell = pow(pow(2,pencils.path[ipencil].size()),2);
      auto ids = pencils.getIds(ipencil);
      for (auto id : ids) {
         cint myCount = std::count(pencils.ids.begin(), pencils.ids.end(), id);
         if (myCount > nPencilsThroughThisCell) {
            std::cerr << "ERROR: Cell ID " << id << " Appears in pencils " << myCount << " times!"<< std::endl;
            std::cerr << "       It should not appear more than " << nPencilsThroughThisCell << " times." << std::endl;
            correct = false;
         }
      }
   }
   return correct;
}

/* Debugging function, prints the list of cells in each pencil
 *
 * @param pencils Pencil data struct
 * @param dimension Spatial dimension
 * @param myRank MPI rank
 */
void printPencilsFunc(const setOfPencils& pencils, const uint dimension, const int myRank,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

// Print out ids of pencils (if needed for debugging)
   uint ibeg = 0;
   uint iend = 0;
   stringstream ss;
   ss << "I am rank " << myRank << ", I have " << pencils.N << " pencils along dimension " << dimension << ":\n";
   MPI_Barrier(MPI_COMM_WORLD);
   if(myRank == MASTER_RANK) {
      ss << "(D=DO_NOT_COMPUTE, S=Sysboundary L2, L=Sysboundary L1, N=Non-sysboundary L2, G=Ghost cell)" << std::endl;
      ss << "t, N, mpirank, dimension, length (x, y): indices {path} DZs AreaRatios" << std::endl;
      ss << "----------------------------------------------------------------------" << std::endl;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   for (uint i = 0; i < pencils.N; i++) {
      const uint L = pencils.lengthOfPencils[i];
      iend = ibeg + L;
      ss << P::t << ", ";
      ss << i << ", ";
      ss << myRank << ", ";
      ss << dimension << ", ";
      ss << L << ", ";
      ss << "(" << pencils.x[i] << ", " << pencils.y[i] << "): ";
      for (auto j = pencils.ids.begin() + ibeg; j != pencils.ids.begin() + iend; ++j) {
         ss << *j;
         if (*j && mpiGrid[*j]) {
            SpatialCell* c = mpiGrid[*j];
            if (c->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) ss<<"D";
            if (c->sysBoundaryLayer != 1 && c->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ss<<"S";
            if (c->sysBoundaryLayer == 1 && c->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ss<<"L";
            if (c->sysBoundaryLayer == 2 && c->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ss<<"N";
            if (!mpiGrid.is_local(*j)) ss<<"G";
         }
         ss<< " ";
      }

      ss << "{";
      for (auto step : pencils.path[i]) {
         ss << step << ", ";
      }
      ss << "}";

      ss << "source DZs: ";
      for (auto j = pencils.sourceDZ.begin() + ibeg; j != pencils.sourceDZ.begin() + iend; ++j) {
         ss << *j << " ";
      }

      ss << "target Ratios: ";
      for (auto j = pencils.targetRatios.begin() + ibeg; j != pencils.targetRatios.begin() + iend; ++j) {
         ss << *j << " ";
      }

      ibeg  = iend;
      ss << std::endl;
   }
   std::cout<<std::flush;
   MPI_Barrier(MPI_COMM_WORLD);
   if(myRank == MASTER_RANK) {
      ss << "-----------------------------------------------------------------" << std::endl;
   }
   std::cout<<ss.str();
}

/* Wrapper function for calling seed ID selection and pencil generation, for all dimensions.
 * Includes threading and gathering of pencils into thread-containers.
 *
 * @param [in] mpiGrid DCCRG grid object
 */
void prepareSeedIdsAndPencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      phiprof::Timer timer {"GetSeedIdsAndBuildPencils"};
      // Remove all old pencils now
      for (int dimension=0; dimension<3; dimension++) {
         DimensionPencils[dimension].removeAllPencils();
      }
      for (int dimension=0; dimension<3; dimension++) {
         prepareSeedIdsAndPencils(mpiGrid, dimension);
      }
}

/* Wrapper function for calling seed ID selection and pencil generation, per dimension.
 * Includes threading and gathering of pencils into thread-containers.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] dimension Spatial dimension
 */
void prepareSeedIdsAndPencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const uint dimension) {

   // Optional heavy printouts for debugging
   const bool printPencils = false;
   const bool printSeeds = false;
   int myRank, mpi_size;
   if(printPencils || printSeeds) {
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
   }

   switch (dimension) {
      case 0:
         if(P::xcells_ini == 1) return;
         break;
      case 1:
         if(P::ycells_ini == 1) return;
         break;
      case 2:
         if(P::zcells_ini == 1) return;
         break;
      default:
         std::cerr<<"Error in dimension: __FILE__:__LINE__"<<std::endl;
         abort();
   }
   const vector<CellID>& localCells = getLocalCells();
   vector<CellID> propagatedCells;
   // Figure out which spatial cells are translated,
   // result independent of particle species.
   if (P::vlasovSolverGhostTranslate) {
      // Sets already include check for do_translate_cell
      switch (dimension) {
         case 0:
            propagatedCells.assign(LocalSet_x.begin(), LocalSet_x.end());
            break;
         case 1:
            propagatedCells.assign(LocalSet_y.begin(), LocalSet_y.end());
            break;
         case 2:
            propagatedCells.assign(LocalSet_z.begin(), LocalSet_z.end());
            break;
         default:
            std::cerr<<"Error in dimension: __FILE__:__LINE__"<<std::endl;
            abort();
      }
   } else {
      for (size_t c=0; c<localCells.size(); ++c) {
         if (do_translate_cell(mpiGrid[localCells[c]])) {
            propagatedCells.push_back(localCells[c]);
         }
      }
   }

   phiprof::Timer getSeedIdsTimer {"getSeedIds"};
   vector<CellID> seedIds;
   getSeedIds(mpiGrid, propagatedCells, dimension, seedIds);
   getSeedIdsTimer.stop();
   
   if (printSeeds) {
      for (int rank=0; rank<mpi_size; ++rank) {
         MPI_Barrier(MPI_COMM_WORLD);
         if (rank!=myRank) continue;
         stringstream ss;
         ss<<"Task "<<myRank<<" Seed Ids (D=DO_NOT_COMPUTE, S=Sysboundary L2, L=Sysboundary L1, N=Non-sysboundary L2, G=Ghost cell)"<<std::endl<<std::endl;
         for (uint i = 0; i < seedIds.size(); i++) {
            ss << seedIds.at(i);
            if (seedIds.at(i) && mpiGrid[seedIds.at(i)]) {
               SpatialCell* c = mpiGrid[seedIds.at(i)];
               if (c->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) ss<<"D";
               if (c->sysBoundaryLayer != 1 && c->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ss<<"S";
               if (c->sysBoundaryLayer == 1 && c->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ss<<"L";
               if (c->sysBoundaryLayer == 2 && c->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ss<<"N";
               if (!mpiGrid.is_local(seedIds.at(i))) ss<<"G";
            }
            ss<<" ";
         }
         ss<<std::endl<<std::endl;
         std::cerr<<ss.str();
      }
   }

   phiprof::Timer buildPencilsTimer {"buildPencils"};

#pragma omp parallel
   {
      // Empty vectors for internal use of buildPencilsWithNeighbors. Could be default values but
      // default vectors are complicated. Should overload buildPencilsWithNeighbors like suggested here
      // https://stackoverflow.com/questions/3147274/c-default-argument-for-vectorint
      vector<CellID> ids;
      vector<uint> path;
      // thread-internal pencil set to be accumulated at the end
      setOfPencils thread_pencils;
      // iterators used in the accumulation
      std::vector<CellID>::iterator ibeg, iend;

#pragma omp for schedule(guided,8)
      for (uint i=0; i<seedIds.size(); i++) {
         cuint seedId = seedIds[i];
         // Construct pencils from the seedIds into a set of pencils.
         buildPencilsWithNeighbors(mpiGrid, thread_pencils, seedId, ids, dimension, path, seedIds);
      }

      // accumulate thread results in global set of pencils
#pragma omp critical
      {
         for (uint i=0; i<thread_pencils.N; i++) {
            // Use vector range constructor
            ibeg = thread_pencils.ids.begin() + thread_pencils.idsStart[i];
            iend = ibeg + thread_pencils.lengthOfPencils[i];
            std::vector<CellID> pencilIds(ibeg, iend);
            DimensionPencils[dimension].addPencil(pencilIds,thread_pencils.x[i],thread_pencils.y[i],thread_pencils.periodic[i],thread_pencils.path[i]);
         }
      }
   }

   phiprof::Timer checkGhostCellsTimer {"check_ghost_cells"};
   // Check refinement of two ghost cells on each end of each pencil
   // in case pencil needs to be split.
   // This function contains threading.
   check_ghost_cells(mpiGrid,DimensionPencils[dimension],dimension);
   checkGhostCellsTimer.stop();

   phiprof::Timer findSourceRatiosTimer {"Find_source_cells_ratios_dz"};
   // Compute also the stencil around the pencil (source cells), and
   // Store source cell widths and target cell contribution ratios.
#pragma omp parallel for schedule(guided)
   for (uint i=0; i<DimensionPencils[dimension].N; ++i) {
      const uint L = DimensionPencils[dimension].lengthOfPencils[i];
      CellID *pencilIds = DimensionPencils[dimension].ids.data() + DimensionPencils[dimension].idsStart[i];
      Realf* pencilDZ = DimensionPencils[dimension].sourceDZ.data() + DimensionPencils[dimension].idsStart[i];
      Realf* pencilAreaRatio = DimensionPencils[dimension].targetRatios.data() + DimensionPencils[dimension].idsStart[i];
      computeSpatialSourceCellsForPencil(mpiGrid,pencilIds,L,dimension,DimensionPencils[dimension].path[i],pencilDZ,pencilAreaRatio);
   }
   findSourceRatiosTimer.stop();

   // ****************************************************************************

   // Now gather unordered_set of target cells (used for resetting block data)
   DimensionTargetCells[dimension].clear();
#pragma omp parallel for
   for (uint i=0; i<DimensionPencils[dimension].ids.size(); ++i) {
      const CellID targ = DimensionPencils[dimension].ids[i];
      const Realf ratio = DimensionPencils[dimension].targetRatios[i];      
      if ((targ!=0)&&(ratio>0.0)) {
#pragma omp critical
         {
            DimensionTargetCells[dimension].insert(targ);
         }
      }
   }

   // Warning: checkPencils fails to understand situations where pencils reach across 3 levels of refinement.
   // if(!checkPencils(mpiGrid, localPropagatedCells, pencils)) {
   //    std::cerr<<"abort checkpencils"<<std::endl;
   //    abort();
   // }

   if(printPencils) {
      for (int rank=0; rank<mpi_size; ++rank) {
         MPI_Barrier(MPI_COMM_WORLD);
         if (rank!=myRank) continue;
         printPencilsFunc(DimensionPencils[dimension],dimension,myRank,mpiGrid);
      }
   }

}
