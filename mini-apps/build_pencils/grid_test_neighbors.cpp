#include "dccrg.hpp"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include "../../definitions.h"
#include "cpu_sort_ids.hpp"

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
  std::vector<uint> lengthOfPencils; // Lengths of pencils
  std::vector<CellID> ids; // List of cells
  std::vector<Real> x,y; // x,y - position (Maybe replace with uint width?)

  setOfPencils() {
    N = 0;
  }

  void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn) {

    N += 1;
    lengthOfPencils.push_back(idsIn.size());
    ids.insert(ids.end(),idsIn.begin(),idsIn.end());
    x.push_back(xIn);
    y.push_back(yIn);
  
  }
    
};


CellID selectNeighbor(dccrg::Dccrg<grid_data> grid, CellID id, int dimension = 0, int path = 0) {

  const auto neighbors = grid.get_face_neighbors_of(id);

  vector < CellID > myNeighbors;
  
  for (const auto cell : neighbors) {
    if (cell.second == dimension + 1)
      myNeighbors.push_back(cell.first);
  }

  CellID retVal;
  
  switch( myNeighbors.size() ) {
  case 0 : {
    // did not find neighbors
    retVal = 0;
  }
    break;
  case 1 : {
    retVal = myNeighbors[0];
  }
    break;
  case 4 : {
    retVal = myNeighbors[path];
  }
    break;
  default: {
    // something is wrong
    retVal = 0;
  }
    break;
  }

  return retVal;
  
}

setOfPencils buildPencilsWithNeighbors( dccrg::Dccrg<grid_data> grid, 
					setOfPencils &pencils, CellID startingId,
					vector<CellID> ids, uint dimension, 
					vector<int> path) {

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
  if( startingRefLvl > 0 ) {
    for ( uint i = 0; i < startingRefLvl; i++) {
      auto parent = grid.get_parent(id);
      auto children = grid.get_all_children(parent);
      auto it = std::find(children.begin(),children.end(),id);
      auto index = std::distance(children.begin(),it);
      auto index2 = index;
      switch( dimension ) {
      case 0: {
	index2 = index / 2;
      }
	break;
      case 1: {
	index2 = index - index / 2;
      }
	break;
      case 2: {
	index2 = index % 4;
      }
	break;
      }
      path.insert(path.begin(),index2);
      id = parent;
    }
  }

  id = startingId;
  
  while (id > 0) {

    // Find the refinement level in the neighboring cell
    nextNeighbor = selectNeighbor(grid,id,dimension);

    if (nextNeighbor == 0)
      break;
    
    uint refLvl = grid.get_refinement_level(nextNeighbor);

    //std::cout << "I am cell " << id << ". Next neighbor is " << nextNeighbor << " at refinement level " << refLvl << std::endl;
    
    if (refLvl == 0 ) {

      // If the neighbor is unrefined, add it to the pencil
      ids.push_back(nextNeighbor);
      
    } else {

      // Check if we have encountered this refinement level before and stored
      // the path this builder follows.
      if ( path.size() >= refLvl ) {
	
	// std::cout << " I have seen refinement level " << refLvl << " before. Path is ";
	// for (auto k = path.begin(); k != path.end(); ++k)
	//   std::cout << *k << ' ';
	// std::cout << std::endl;

	nextNeighbor = selectNeighbor(grid,id,dimension,path[refLvl-1]);
	if (nextNeighbor == 0)
	  break;
	
	ids.push_back(nextNeighbor);
	
      } else {

	// std::cout << " I have NOT seen refinement level " << refLvl << " before. Path is ";
	// for (auto k = path.begin(); k != path.end(); ++k)
	//   std::cout << *k << ' ';
	// std::cout << std::endl;

	
	// New refinement level, create a path through each neighbor cell
	for ( uint i : {0,1,2,3} ) {
	  
	  vector < int > myPath = path;
	  myPath.push_back(i);
	  
	  nextNeighbor = selectNeighbor(grid,id,dimension,myPath.back());
	  if (nextNeighbor == 0)
	    break;
	  
	  if ( i == 3 ) {

	    // This builder continues along neighbor 3
	    ids.push_back(nextNeighbor);
	    path = myPath;
	    
	  } else {

	    // Spawn new builders for neighbors 0,1,2
	    buildPencilsWithNeighbors(grid,pencils,id,ids,dimension,myPath);
	    
	  }
	  
	}
	
      }

    } // Closes if (refLvl == 0)

    id = nextNeighbor;
    
  } // Closes while loop
  
  pencils.addPencil(ids,0.0,0.0);
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

  const uint xDim = 9;
  const uint yDim = 3;
  const uint zDim = 1;
  const std::array<uint64_t, 3> grid_size = {{xDim,yDim,zDim}};
  
  grid.initialize(grid_size, comm, "RANDOM", 1);

  grid.balance_load();

  bool doRefine = true;
  const std::array<uint,4> refinementIds = {{10,14,64,72}};
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
  
  int dimension = 0;
  
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
  vector<int> path;  
  setOfPencils pencils;

  for (auto id : startingIds) {
    pencils = buildPencilsWithNeighbors(grid,pencils,id,idsInitial,dimension,path);
  }

  uint ibeg = 0;
  uint iend = 0;
  
  std::cout << "I have created " << pencils.N << " pencils along dimension " << dimension << ":\n";
  for (uint i = 0; i < pencils.N; i++) {
    iend += pencils.lengthOfPencils[i];
    for (auto j = pencils.ids.begin() + ibeg; j != pencils.ids.begin() + iend; ++j) {
      std::cout << *j << " ";
    }
    ibeg  = iend;
    std::cout << "\n";
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
