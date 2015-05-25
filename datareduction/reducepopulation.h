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
void population_algorithm(
                      spatial_cell::SpatialCell * cell,
                      const Real resolution_threshold=0.006
                         );


/*! Writes the population variables, so variables for different populations. Note that population_algorithm must have been called before this!

 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already
 */
//bool write_population_variables( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, vlsv::Writer & vlsvWriter );

/*! Writes the distribution function for different populations. Note that population_algorithm must have been called before this!

 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already
 \param local_cells             The cells for which we write out the velocity space
 */
//bool write_population_distribution( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const vector<uint64_t> & local_cells, vlsv::Writer & vlsvWriter );
