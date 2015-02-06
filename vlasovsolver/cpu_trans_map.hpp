/*
  This file is part of Vlasiator.
  Copyright 2014-2015 Finnish Meteorological Institute
*/
#ifndef CPU_TRANS_MAP_H
#define CPU_TRANS_MAP_H

#include <vector>

#include "vec.h"
#include "../common.h"
#include "../spatial_cell.hpp"

void clearTargetGrid(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const std::vector<CellID>& cells);
void createTargetGrid(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const std::vector<CellID>& cells);
bool do_translate_cell(spatial_cell::SpatialCell* SC);
void swapTargetSourceGrid(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const std::vector<CellID>& cells);
bool trans_map_1d(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const CellID cellID,const uint dimension,const Realv dt);
void update_remote_mapping_contribution(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint dimension,int direction);
void zeroTargetGrid(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const std::vector<CellID>& cells);

#endif
