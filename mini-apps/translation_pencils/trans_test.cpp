#include <iostream>
#include <fstream>
#include <vector>
//#include "dccrg.hpp"
#include "../../grid.h"
#include "mpi.h"
#include "../../definitions.h"
#include "../../parameters.h"
#include "../../vlasovsolver/cpu_trans_map.hpp"

using namespace std;

// struct grid_data {

//   int value = 0;

//   std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
//   {
//     return std::make_tuple(this, 0, MPI_BYTE);
//   }
    
// };

int main(int argc, char* argv[]) {

   if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      // cerr << "Coudln't initialize MPI." << endl;
      abort();
   }

   MPI_Comm comm = MPI_COMM_WORLD;

   int rank = 0, comm_size = 0;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &comm_size);

   const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry> grid;

   const uint dimension = 0;
   const uint xDim = 9;
   const uint yDim = 9;
   const uint zDim = 1;
   const std::array<uint64_t, 3> grid_size = {{xDim,yDim,zDim}};

   int argn;
   char **argc;
   
   initializeGrid(argn,argc,grid,sysBoundaries,project);
   //grid.initialize(grid_size, comm, "RANDOM", 1);

   grid.balance_load();

   bool doRefine = false;
   const std::array<uint,4> refinementIds = {{1,2,3,4}};
   if(doRefine) {
      for(uint i = 0; i < refinementIds.size(); i++) {
         if(refinementIds[i] > 0) {
            grid.refine_completely(refinementIds[i]);
            grid.stop_refining();
         }
      }
   }

   grid.balance_load();

   setOfPencils pencils;
   vector<CellID> seedIds;
   vector<CellID> localPropagatedCells;
   vector<CellID> ids;
   vector<uint> path;
   
   for (CellID i = 0; i < xDim * yDim * zDim; i++) localPropagatedCells.push_back( i + 1 );
   get_seed_ids(grid, localPropagatedCells, dimension, seedIds);
   for (const auto seedId : seedIds) {
      // Construct pencils from the seedIds into a set of pencils.
      pencils = buildPencilsWithNeighbors(grid, pencils, seedId, ids, dimension, path);      
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

}
