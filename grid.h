#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include "datareducer.h"

//Init parallel grid
bool initializeGrid(int argn, char **argc,dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);
//Balance load
void balanceLoad(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);

//write out system
bool writeGrid(const dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,DataReducer& dataReducer,const bool& writeRestart);

//Adjust blocks; remove/add based on user-defined limits
bool adjust_all_velocity_blocks(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);

/*
Updates velocity block lists between remote neighbors and prepares local
copies of remote neighbors to receive velocity block data.
*/
void prepare_to_receive_velocity_block_data(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);

#endif
