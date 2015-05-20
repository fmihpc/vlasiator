#include <utility>
#include <vector>
#include <array>
#include "vlsv_writer.h"
#include "../spatial_cell.hpp"
#include "../grid.h"


/*! A function for calculating local and remote velocity cells for each velocity cell within a block. Remote velocity cells are the cells that would be outside a block and local velocity cells are the cells that are inside a block. Note: This is currently only used in the population_algorithm as a way to optimize the code. This way we don't need to calculate the neighbors for each velocity cell every time, and the memory consumption is negligible. After running this one can retrieve a given velocity cell's neighbors by the velocity cell id (0..63) with local_vcell_neighbors[vCellId] and remote neighbors likewise with remote_vcell_neighbors[vCellId]. The remote velocity cells are stored in pairs of blockid and list of neighbors within that block. So, for instance if remote_vcell_neighbors[vCellId] returned this: [<2,[0,4]> <10,[3,4]>] where <,> signals a pair, then the vCellId at block blockId would have remtoe neighbors at block blockId+2 at velocity cells 0 and 4 and at block blockId+10 at velocity cells 3 and 4

 \param local_vcell_neighbors                 Empty vector where local velocity cells get stored.  local_vcell_neighbors[vCellId] gives the local neighbors of a velocity cell (neighbors that are within the same block)
 \param remote_vcell_neighbors                Empty vector where remote velocity cells get stored remote_vcell_neighbors[vCellId] Gives a vector containing remote neighbors of the velocity cell (neighbors that are outside the block) in vector< pair<int16_t, vector<uint16_t> > > format. The pair's first index gives the neighbor block index (check spatial_cell.hpp for more information) and the second index gives the local velocity cells within that neighbor block.

 */
void set_local_and_remote_velocity_cell_neighbors(
       vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> & vmesh,
       std::array<std::vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       std::array< std::vector< std::pair<int16_t, std::vector<uint16_t> > >, VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                 );
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
