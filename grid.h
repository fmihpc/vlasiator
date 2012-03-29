#ifndef GRID_H
#define GRID_H

#include "spatial_cell.hpp"
#include <dccrg.hpp>
#include "datareducer.h"

//Init parallel grid
bool initializeGrid(int argn, char **argc,dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);
//Balance load
void balanceLoad(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);

//Adjust which blocks are existing based on the current state. 
bool adjust_velocity_blocks(dccrg::Dccrg<SpatialCell>& mpiGrid);

//write out system
bool writeGrid(const dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,DataReducer& dataReducer,const bool& writeRestart);


#endif
