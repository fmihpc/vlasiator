#include <utility>
#include <vector>
#include <array>
#include "vlsv_writer.h"
#include "../spatial_cell.hpp"
#include "../grid.h"


void set_local_and_remote_velocity_cell_neighbors(
       std::array<std::vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       std::array< std::vector< std::pair<int16_t, std::vector<uint16_t> > >, VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                 );
Real evaluate_speed(
                spatial_cell::SpatialCell * cell,
                const std::array<std::vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const std::array< std::vector< std::pair<int16_t, std::vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                   );
bool write_population_variables( dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, vlsv::Writer & vlsvWriter );

