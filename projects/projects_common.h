/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

*/

#ifndef PROJECTS_COMMON_H
#define PROJECTS_COMMON_H
#include "../spatial_cell.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

namespace projects {
   enum Neighbours {
      ZM1_YM1_XM1,
      ZM1_YM1_XCC,
      ZM1_YM1_XP1,
      ZM1_YCC_XM1,
      ZM1_YCC_XCC,
      ZM1_YCC_XP1,
      ZM1_YP1_XM1,
      ZM1_YP1_XCC,
      ZM1_YP1_XP1,
      ZCC_YM1_XM1,
      ZCC_YM1_XCC,
      ZCC_YM1_XP1,
      ZCC_YCC_XM1,
      ZCC_YCC_XCC,
      ZCC_YCC_XP1,
      ZCC_YP1_XM1,
      ZCC_YP1_XCC,
      ZCC_YP1_XP1,
      ZP1_YM1_XM1,
      ZP1_YM1_XCC,
      ZP1_YM1_XP1,
      ZP1_YCC_XM1,
      ZP1_YCC_XCC,
      ZP1_YCC_XP1,
      ZP1_YP1_XM1,
      ZP1_YP1_XCC,
      ZP1_YP1_XP1
   };
   
   const uint MISSING_ZNEG = (1 << projects::ZM1_YCC_XCC);
   const uint MISSING_YNEG = (1 << projects::ZCC_YM1_XCC);
   const uint MISSING_XNEG = (1 << projects::ZCC_YCC_XM1);
   const uint MISSING_XPOS = (1 << projects::ZCC_YCC_XP1);
   const uint MISSING_YPOS = (1 << projects::ZCC_YP1_XCC);
   const uint MISSING_ZPOS = (1 << projects::ZP1_YCC_XCC);
   const uint FACE_NBR_BITMASK = (MISSING_ZNEG | MISSING_YNEG | MISSING_XNEG | MISSING_XPOS | MISSING_YPOS | MISSING_ZPOS);
}

// *********************************
// ***** TEMPLATE DECLARATIONS *****
// *********************************

template<typename CELLID,class CONT> bool classifyLevequeGhostCell(const SpatialCell& cell,const CELLID& cellID,const CONT& nbrs);

template<typename CELLID> CELLID getNeighbour(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const CELLID& cellID,const int& i,const int& j,const int& k);


// ********************************
// ***** TEMPLATE DEFINITIONS *****
// ********************************

template<typename CELLID> CELLID getNeighbour(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const CELLID& cellID,const int& i,const int& j,const int& k){
    std::vector<uint64_t> neighbors = mpiGrid.get_neighbors_of_at_offset(cellID, i, j, k);

    //FIXME: support refined grids
    if(neighbors.size() > 0) {
        return neighbors[0];
    } else {
        return INVALID_CELLID;
    }
}


#endif
