#include "dccrg.hpp"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include "../../definitions.h"
#include "../../parameters.h"

using namespace std;

struct grid_data {

   int value = 0;

   std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
   {
      return std::make_tuple(this, 0, MPI_BYTE);
   }
    
};

struct setOfPencils {

   uint N; // Number of pencils in the set
   uint sumOfLengths;
   std::vector<uint> lengthOfPencils; // Lengths of pencils
   std::vector<CellID> ids; // List of cells
   std::vector<Real> x,y; // x,y - position 

   setOfPencils() {
      N = 0;
      sumOfLengths = 0;
   }

   void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn) {

      N += 1;
      sumOfLengths += idsIn.size();
      lengthOfPencils.push_back(idsIn.size());
      ids.insert(ids.end(),idsIn.begin(),idsIn.end());
      x.push_back(xIn);
      y.push_back(yIn);
  
   }

   std::vector<CellID> getIds(uint pencilId) {

      if (pencilId > N) {
         vector<CellID> foo;
         return foo;
      }

      CellID ibeg = 0;
      for (uint i = 0; i < pencilId; i++) {
         ibeg += lengthOfPencils[i];
      }
      CellID iend = ibeg + lengthOfPencils[pencilId];

      vector<CellID> idsOut;
    
      for (uint i = ibeg; i <= iend; i++) {
         idsOut.push_back(ids[i]);
      }

      return idsOut;
   }

};


CellID selectNeighbor(dccrg::Dccrg<grid_data> &grid, CellID id, int dimension = 0, uint path = 0) {

   vector < CellID > myNeighbors;
   // Collect neighbor ids in the positive direction of the chosen dimension.
   // Note that dimension indexing starts from 1 (of course it does)
   for (const auto& [neighbor, dir] : grid.get_face_neighbors_of(id)) {
      if (cell.second == dimension + 1)
         myNeighbors.push_back(cell.first);
   }

   CellID neighbor;
  
   switch( myNeighbors.size() ) {
      // Since refinement can only increase by 1 level the only possibilities
      // Should be 0 neighbors, 1 neighbor or 4 neighbors.
   case 0 : {
      // did not find neighbors
      neighbor = INVALID_CELLID;
      break;
   }
   case 1 : {
      neighbor = myNeighbors[0];
      break;
   }
   case 4 : {
      neighbor = myNeighbors[path];
      break;
   }
   default: {
      // something is wrong
      neighbor = INVALID_CELLID;
      throw "Invalid neighbor count!";
      break;
   }
   }

   return neighbor;
  
}

setOfPencils buildPencilsWithNeighbors( dccrg::Dccrg<grid_data> &grid, 
					setOfPencils &pencils, CellID startingId,
					vector<CellID> ids, uint dimension, 
					vector<uint> path) {

   const bool debug = false;
   CellID nextNeighbor;
   uint id = startingId;
   uint startingRefLvl = grid.get_refinement_level(id);

   if( ids.size() == 0 )
      ids.push_back(startingId);

   // If the cell where we start is refined, we need to figure out which path
   // to follow in future refined cells. This is a bit hacky but we have to
   // use the order or the children of the parent cell to figure out which
   // corner we are in.
   // Maybe you could use physical coordinates here?
   if( startingRefLvl > path.size() ) {
      for ( uint i = path.size(); i < startingRefLvl; i++) {
         auto parent = grid.mapping.get_parent(id);

         std::array<uint64_t, 8> childrenarr = mpiGrid.mapping.get_all_children(parent);
         vector<CellID> children(childrenarr.begin(), childrenarr.end());
         auto it = std::find(children.begin(),children.end(),id);
         auto index = std::distance(children.begin(),it);
         auto index2 = index;
      
         switch( dimension ) {
         case 0: {
            index2 = index / 2;
            break;
         }
         case 1: {
            index2 = index - index / 2;
            break;
         }
         case 2: {
            index2 = index % 4;
            break;
         }
         }
         path.insert(path.begin(),index2);
         id = parent;
      }
   }

   id = startingId;
  
   while (id > 0) {

      // Find the refinement level in the neighboring cell. Any neighbor will do
      // since refinement level can only increase by 1 between neighbors.
      nextNeighbor = selectNeighbor(grid,id,dimension);

      // If there are no neighbors, we can stop.
      if (nextNeighbor == 0)
         break;
    
      uint refLvl = grid.get_refinement_level(nextNeighbor);

      if (refLvl > 0) {
    
         // If we have encountered this refinement level before and stored
         // the path this builder follows, we will just take the same path
         // again.
         if ( path.size() >= refLvl ) {
      
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << " ";
               std::cout << std::endl;
            }
	
            nextNeighbor = selectNeighbor(grid,id,dimension,path[refLvl-1]);      
	
         } else {
	
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have NOT seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << ' ';
               std::cout << std::endl;
            }
	
            // New refinement level, create a path through each neighbor cell
            for ( uint i : {0,1,2,3} ) {
	  
               vector < uint > myPath = path;
               myPath.push_back(i);
	  
               nextNeighbor = selectNeighbor(grid,id,dimension,myPath.back());
	  
               if ( i == 3 ) {
	    
                  // This builder continues with neighbor 3
                  ids.push_back(nextNeighbor);
                  path = myPath;
	    
               } else {
	    
                  // Spawn new builders for neighbors 0,1,2
                  buildPencilsWithNeighbors(grid,pencils,id,ids,dimension,myPath);
	    
               }
	  
            }
	
         }

      } else {
         if(debug) {
            std::cout << "I am cell " << id << ". ";
            std::cout << " I am on refinement level 0." << std::endl;
         }
      }// Closes if (refLvl == 0)

      // If we found a neighbor, add it to the list of ids for this pencil.
      if(nextNeighbor != INVALID_CELLID) {
         if (debug) {
            std::cout << " Next neighbor is " << nextNeighbor << "." << std::endl;
         }
         ids.push_back(nextNeighbor);
      }

      // Move to the next cell.
      id = nextNeighbor;
    
   } // Closes while loop

   // Get the x,y - coordinates of the pencil (in the direction perpendicular to the pencil)
   const auto coordinates = grid.get_center(ids[0]);
   double x,y;
   uint ix,iy,iz;
   switch(dimension) {
   case 0: 
      ix = 1;
      iy = 2;
      iz = 0;
      break;
   
   case 1: 
      ix = 2;
      iy = 0;
      iz = 1;
      break;
   
   case 2: 
      ix = 0;
      iy = 1;
      iz = 2;
      break;
   
   default: 
      ix = 0;
      iy = 1;
      iz = 2;
      break;
   
   }
   
   x = coordinates[ix];
   y = coordinates[iy];
   // z = vector<Real>;
  
   // for( auto id: ids ) {
   //   coordinates = grid.get_center(id);
   //   z.push_back(coordinates[iz])
   // }
  
   pencils.addPencil(ids,x,y);
   return pencils;
  
}


void printVector(vector<CellID> v) {

   for (auto k = v.begin(); k != v.end(); ++k)
      std::cout << *k << ' ';
   std::cout <<  "\n";

}

int main(int argc, char* argv[]) {

   if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      // cerr << "Coudln't initialize MPI." << endl;
      abort();
   }
  
   MPI_Comm comm = MPI_COMM_WORLD;
  
   int rank = 0, comm_size = 0;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &comm_size);
  
   dccrg::Dccrg<grid_data> grid;

   // paremeters
   const uint xDim = 9;
   const uint yDim = 3;
   const uint zDim = 1;
   const std::array<uint64_t, 3> grid_size = {{xDim,yDim,zDim}};
   const int dimension = 0;
   const bool doRefine = true;
   const std::array<uint,4> refinementIds = {{10,14,64,72}};

   grid.initialize(grid_size, comm, "RANDOM", 1);
   
   typedef dccrg::Types<3>::neighborhood_item_t neigh_t;
   std::vector<neigh_t> neighborhood_x;
   std::vector<neigh_t> neighborhood_y;

   int neighborhood_width = 2;
   for (int d = -neighborhood_width; d <= neighborhood_width; d++) {
     if (d != 0) {
        neighborhood_x.push_back({{d, 0, 0}});
        neighborhood_y.push_back({{0, d, 0}});
     }
   }
   grid.add_neighborhood(1, neighborhood_x);
   grid.add_neighborhood(2, neighborhood_y);      

   grid.balance_load();

   if(doRefine) {
      for(uint i = 0; i < refinementIds.size(); i++) {
         if(refinementIds[i] > 0) {
            grid.refine_completely(refinementIds[i]);
            grid.stop_refining();
         }
      }
   }
   
   grid.balance_load();

   auto cells = grid.cells;
   sort(cells.begin(), cells.end());

   vector<CellID> ids;
   vector<CellID> startingIds;  
  
   for (const auto& cell: cells) {
      // std::cout << "Data of cell " << cell.id << " is stored at " << cell.data << std::endl;
      // Collect a list of cell ids.
      ids.push_back(cell.id);

      // Collect a list of cell ids that do not have a neighbor in the negative direction
      vector<CellID> negativeNeighbors;
      for (auto neighbor : grid.get_face_neighbors_of(cell.id)) {      
      
         if (neighbor.second == - (dimension + 1))
            negativeNeighbors.push_back(neighbor.first);
      }
      if (negativeNeighbors.size() == 0)
         startingIds.push_back(cell.id);    
   }

   std::cout << "Starting cell ids for pencils are ";
   printVector(startingIds);

   sort(ids.begin(),ids.end());     
  
   vector<CellID> idsInitial;
   vector<uint> path;  
   setOfPencils pencils;

   for (auto id : startingIds) {
      pencils = buildPencilsWithNeighbors(grid,pencils,id,idsInitial,dimension,path);
   }

   uint ibeg = 0;
   uint iend = 0;
  
   std::cout << "I have created " << pencils.N << " pencils along dimension " << dimension << ":\n";
   std::cout << "(x, y): indices " << std::endl;
   std::cout << "-----------------------------------------------------------------" << std::endl;
   for (uint i = 0; i < pencils.N; i++) {
      iend += pencils.lengthOfPencils[i];
      std::cout << "(" << pencils.x[i] << ", " << pencils.y[i] << "): ";
      for (auto j = pencils.ids.begin() + ibeg; j != pencils.ids.begin() + iend; ++j) {
         std::cout << *j << " ";
      }
      ibeg  = iend;
      std::cout << std::endl;
   }


   CellID id = 3;
   const vector<CellID>* neighbors = grid.get_neighbors_of(id, 1);
   if (neighbors != NULL) {
      std::cout << "Neighbors of cell " << id << std::endl;
      for (auto neighbor : *neighbors) {
         std::cout << neighbor << std::endl;
      }
   }


   
   std::ofstream outfile;
  
   grid.write_vtk_file("test.vtk");

   outfile.open("test.vtk", std::ofstream::app);
   // write each cells id
   outfile << "CELL_DATA " << cells.size() << std::endl;
   outfile << "SCALARS id int 1" << std::endl;
   outfile << "LOOKUP_TABLE default" << std::endl;
   for (const auto& cell: cells) {
      outfile << cell.id << std::endl;
   }
   outfile.close();
		
   MPI_Finalize();

   return 0;
    
}
