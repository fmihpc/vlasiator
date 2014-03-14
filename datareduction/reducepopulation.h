#include <utility>
#include <vector>
#include <array>
#include "../spatial_cell.hpp"


void set_local_and_remote_velocity_cell_neighbors(
       array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       array< vector< pair<uint16_t, vector<uint16_t> > > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                 );
Real evaluate_speed(
                const SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                   );
Real evaluate_speed_parallel(
                const SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                            );
