#include <utility>
#include <vector>
#include <array>
#include "vlsv_writer.h"
#include "../spatial_cell.hpp"
#include "../grid.h"


/*! Function for calculating the different populations in distribution for a given spatial cell. Note that the results are currently saved in block_fx array within the SpatialCell.

 \param cell                                A cell from the DCCRG grid
 \param resolution_threshold                A value for determining how large at minimum we want our populations to be. 0.006 seems ok unless there's a reason to believe otherwise.

 */
void populationAlgorithm(
                      spatial_cell::SpatialCell * cell,
                      const Real resolution_threshold=0.1
                         );


/*! Function for writing a data reducer for population merger (writes variables for each population)
 \param local_cells                        Vector of cell ids in the local process domain
 \param mpiGrid                            DCCRG grid which contains pointers to cell data
 \param vlsvWriter                         a VLSV file with a file opan
 * */
/*
bool writePopulationDataReducer(
                               const vector<CellID>& local_cells,
                               dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, 
                               Writer & vlsvWriter
                               );
*/