/*
This file is part of Vlasiator.
Copyright 2013-2015 Finnish Meteorological Institute
*/

#ifndef CPU_ACC_MAP_H
#define CPU_ACC_MAP_H

#include "../common.h"
#include "../spatial_cell.hpp"
#include "vec.h"

bool map_1d(vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,
            Realv intersection,Realv intersection_di,Realv intersection_dj,Realv intersection_dk,
            uint dimension);

#endif
