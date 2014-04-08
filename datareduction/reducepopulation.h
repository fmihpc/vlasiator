#include <utility>
#include <vector>
#include <array>
#include "../spatial_cell.hpp"


void set_local_and_remote_velocity_cell_neighbors(
       std::array<std::vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       std::array< std::vector< std::pair<uint16_t, std::vector<uint16_t> > >, VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                 );
Real evaluate_speed(
                spatial_cell::SpatialCell * cell,
                const std::array<std::vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const std::array< std::vector< std::pair<uint16_t, std::vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                   );

