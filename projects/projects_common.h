/*
This file is part of Vlasiator.

Copyright 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PROJECTS_COMMON_H
#define PROJECTS_COMMON_H
using namespace spatial_cell;
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


#include <dccrg.hpp>
template<typename CELLID> CELLID getNeighbour(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CELLID& cellID,const int& i,const int& j,const int& k);


// ********************************
// ***** TEMPLATE DEFINITIONS *****
// ********************************

template<typename CELLID> CELLID getNeighbour(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CELLID& cellID,const int& i,const int& j,const int& k){
    std::vector<uint64_t> neighbors = mpiGrid.get_neighbors_of(cellID, i, j, k);

    //FIXME: support refined grids
    if(neighbors.size() > 0) {
        return neighbors[0];
    } else {
        return INVALID_CELLID;
    }
}


#endif
