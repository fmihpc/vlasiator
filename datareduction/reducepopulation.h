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

