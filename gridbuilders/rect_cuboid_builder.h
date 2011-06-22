#ifndef RECT_CUBOID_BUILDER_H
#define RECT_CUBOID_BUILDER_H

#include <vector>
#include <map>
#include "mpibuilder.h"

// This GridBuilder creates a rectangular cuboid 
// (see, e.g. Wikipedia). In short, a box-shaped grid.

class RectCuboidBuilder: public MPIBuilder {
 public:
   RectCuboidBuilder();
   ~RectCuboidBuilder();
   
   bool finalize();

   bool calculatesAnalyticInitialState();
   bool getCellBlockData(const VirtualCell::ID& cellID,cuint& N_blocks,Real* blocks,Real* blockParams,uint* nbrsVel);
   bool getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs);
   bool getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* nbrTypes);
   bool getCellNumberOfBlocks(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,uint* N_blocks);
   bool getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(VirtualCell::ID& N_cells);
   bool initialize(MPI_Comm comm,const int& MASTER_RANK);
   
 protected:
   bool initialized;  /**< If true, RectCuboidBuilder initialized successfully.*/
   
   Real dx;           /**< Size of spatial cells in x-direction.*/
   Real dy;           /**< Size of spatial cells in y-direction.*/
   Real dz;           /**< Size of spatial cells in z-direction.*/
   Real xmin;         /**< Minimum value of spatial grid x-coordinate.*/
   Real xmax;         /**< Maximum value of spatial grid x-coordinate.*/
   Real ymin;         /**< Minimum value of spatial grid y-coordinate.*/
   Real ymax;         /**< Maximum value of spatial grid y-coordinate.*/
   Real zmin;         /**< Minimum value of spatial grid z-coordinate.*/
   Real zmax;         /**< Maximum value of spatial grid z-coordinate.*/
   uint xsize;        /**< Initial number of cells in x-direction.*/
   uint ysize;        /**< Initial number of cells in y-direction.*/
   uint zsize;        /**< Initial number of cells in z-direction.*/
   
   uint vx_blocks;    /**< Initial number of velocity blocks in vx-direction.*/
   uint vy_blocks;    /**< Initial number of velocity blocks in vy-direction.*/
   uint vz_blocks;    /**< Initial number of velocity blocks in vz-direction.*/
   Real vx_min;       /**< Minimum value of velocity block vx coordinate.*/
   Real vx_max;       /**< Maximum value of velocity block vx coordinate.*/
   Real vy_min;       /**< Minimum value of velocity block vy coordinate.*/
   Real vy_max;       /**< Maximum value of velocity block vy coordinate.*/
   Real vz_min;       /**< Minimum value of velocity block vz coordinate.*/
   Real vz_max;       /**< Maximum value of velocity block vz coordinate.*/

   bool periodicInX;  /**< If true, grid is periodic in x-direction.*/
   bool periodicInY;  /**< If true, grid is periodic in y-direction.*/
   bool periodicInZ;  /**< If true, grid is periodic in z-direction.*/
   
   std::map<std::string,std::string> options;

   VirtualCell::ID calculateNeighbourID(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k,
					const int& i_nbr,const int& j_nbr,const int& k_nbr);
   uint calculateNeighbours(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k,
			    VirtualCell::ID& x_neg,VirtualCell::ID& x_pos,VirtualCell::ID& y_neg,
			    VirtualCell::ID& y_pos,VirtualCell::ID& z_neg,VirtualCell::ID& z_pos);
   uint countNeighbours(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k);
   VirtualCell::ID spatCellIndex(const VirtualCell::ID& i,const VirtualCell::ID& j,const VirtualCell::ID& k);
   uint velBlockIndex(cuint& iv,cuint& jv,cuint& kv);
};

#endif
