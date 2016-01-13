/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * 
 */

#ifndef FS_CACHE_H
#define FS_CACHE_H

#include <vector>

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"

namespace fs_cache {
   /*
   enum Element {
        C222, // parameters,derivatives
        C122, // parameters,derivatives
        C322, // parameters,derivatives
        C212, // parameters,derivatives
        C232, // parameters,derivatives
        C221, // parameters,derivatives
        C223, // parameters,derivatives
        C233, // parameters
        C323, // parameters
        C332, // parameters
        C211, // parameters,derivatives
        C311, // parameters,derivatives
        C312, // parameters,derivatives
        C321, // parameters,derivatives
        N_ELEMENTS
   };*/
         
   template<typename INT> INT calculateNbrID(const INT& I,const INT& J,const INT& K) {
      return K*9 + J*3 + I;
   }
    
   struct CellCache {
      CellCache();
      CellCache(Real** ptrs,Real** derivs,Real* derivsBVOL);
      
      CellID cellID;
      uint existingCellsFlags;
      uint sysBoundaryFlag;
      //Real* parameters[N_ELEMENTS];
      //Real* derivatives[N_ELEMENTS];
      //Real* derivativesBVOL;

      spatial_cell::SpatialCell* cells[27];
   };
   
   struct CacheContainer {
      static long int cacheCalculatedStep;
      static std::vector<fs_cache::CellCache> localCellsCache;        /**< Cache for all local cells.*/
      static std::vector<uint16_t> boundaryCellsWithLocalNeighbours;
      static std::vector<uint16_t> boundaryCellsWithRemoteNeighbours;
      static std::vector<uint16_t> cellsWithLocalNeighbours;          /**< Local IDs of get_local_cells_not_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID) cells.
                                                                       * Cells with sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE are not included.
                                                                       * Stored values are used to index CacheContainer::localCellsCache.*/
      
      static std::vector<uint16_t> cellsWithRemoteNeighbours;         /**< Local IDs of get_local_cells_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID) cells.
                                                                       * Cells with sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE are not included.
                                                                       * Stored values are used to index CacheContainer::localCellsCache.*/
      static std::vector<uint16_t> local_NOT_DO_NOT_COMPUTE;          /**< Exclude DO_NOT_COMPUTE cells.*/
      static std::vector<uint16_t> local_NOT_SYSBOUND_DO_NOT_COMPUTE; /**< Exclude DO_NOT_COMPUTE and system boundary cells.*/
      
      static void clear();
   };
   
   void calculateCache(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells
   );
   
   CacheContainer& getCache();

} // namespace fs_cache
   
#endif
