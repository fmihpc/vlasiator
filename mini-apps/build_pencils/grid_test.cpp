#include "dccrg.hpp"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include "../../definitions.h"
#include <list>
#include "cpu_sort_ids.hpp"
#include <map>

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

  void addPencil(vector<CellID> idsIn, Real xIn, Real yIn) {

    N += 1;
    lengthOfPencils.push_back(idsIn.size());
    ids.insert(ids.end(),idsIn.begin(),idsIn.end());
    x.push_back(xIn);
    y.push_back(yIn);
  
  }
    
};


  
void insertVectorIntoVector(vector<CellID> &v1,vector<CellID> v2,uint i) {

  vector<CellID> tmp(v1.begin(),v1.begin() + i);
  tmp.insert(tmp.end(),v2.begin(),v2.end());
  tmp.insert(tmp.end(),v1.begin() + i, v1.end());
  v1.clear();
  v1 = tmp;
  v2.clear();
  
}

vector<CellID> getMyChildren(vector<CellID> &children,
		    uint dimension, bool up, bool left ) {

  bool down  = !up;
  bool right = !left;
  
  uint i1 = 999;
  uint i2 = 999;
  
  switch(dimension) {
  case 0 :
    if (up && left) {
      i1 = 0;
      i2 = 1;
      break;
    }
    if (down && left) {
      i1 = 2;
      i2 = 3;
      break;
    }
    if (up && right) {
      i1 = 4;
      i2 = 5;
      break;
    }
    if (down && right) {
      i1 = 6;
      i2 = 7;
      break;
    }
  case 1 :
    if (up && left) {
      i1 = 0;
      i2 = 2;
      break;
    }
    if (down && left) {
      i1 = 1;
      i2 = 3;
      break;
    }
    if (up && right) {
      i1 = 4;
      i2 = 6;
      break;
    }
    if (down && right) {
      i1 = 5;
      i2 = 7;
      break;
    }
  case 2 :
    if (up && left) {
      i1 = 0;
      i2 = 4;
      break;
    }
    if (down && left) {
      i1 = 1;
      i2 = 5;
      break;
    }
    if (up && right) {
      i1 = 2;
      i2 = 6;
      break;
    }
    if (down && right) {
      i1 = 3;
      i2 = 7;
      break;
    }
  default:    
    break;

  }
  
  vector<CellID> myChildren {children[i1],children[i2]};
  return myChildren;
  
}

void printVector(vector<CellID> v) {

  for (auto k = v.begin(); k != v.end(); ++k)
    std::cout << *k << ' ';
  std::cout <<  "\n";

}

setOfPencils buildPencils( dccrg::Dccrg<grid_data> grid,
			   setOfPencils &pencils, vector<CellID> idsOut,
			   vector<CellID> idsIn, int dimension, 
			   vector<pair<bool,bool>> path) {

  // Not necessary since c++ passes a copy by default.
  // Copy the input ids to a working set of ids
  // vector<int> ids( idsIn );
  // Copy the already computed pencil to the output list
  // vector<int> idsOut( idsInPencil );
  
  uint i = 0;
  uint length = idsIn.size();
  
  // Walk along the input pencil. Initially length is equal to the length of the
  // Unrefined pencil. When refined cells are encountered, the length is increased
  // accordingly to go through the entire pencil.
  while (i < length) {

    uint i1 = i + 1;
    uint id = idsIn[i];

    
    std::array<uint64_t, 8> children = mpiGrid.mapping.get_all_children(id);
    bool hasChildren = ( grid.mapping.get_parent(children[0]) == id );
    
    // Check if the current cell contains refined cells
    if (hasChildren) {
      
      // Check if we have encountered this refinement level before and stored
      // the path this builder followed
      if (path.size() > grid.get_refinement_level(id)) {

	// Get children using the stored path	
	vector<CellID> myChildren = getMyChildren(children,dimension,
						  path[grid.get_refinement_level(id)].first,
						  path[grid.get_refinement_level(id)].second);
	
	// Add the children to the working set at index i1

	insertVectorIntoVector(idsIn,myChildren,i1);
	length += myChildren.size();	
      
      } else {

	// Spawn new builders to construct pencils at the new refinement level
	
	for (bool left : { true, false }) {	
	  for (bool up : { true, false }) {
	    
	    // Store the path this builder has chosen
	    vector < pair <bool,bool>> myPath = path;
	    myPath.push_back(pair<bool, bool>(up,left));
	    
	    // Get children along my path.
	    vector<CellID> myChildren = getMyChildren(children,dimension,up,left);
	    // Get the ids that have not been processed yet.
	    vector<CellID> remainingIds(idsIn.begin() + i1, idsIn.end());
	    
	    // The current builder continues along the bottom-right path.
	    // Other paths will spawn a new builder.
	    if (!up && !left) {
	      
	      // Add the children to the working set. Next iteration of the
	      // main loop (over idsIn) will start on the first child

	      // Add the children to the working set at index i1
	      insertVectorIntoVector(idsIn,myChildren,i1);
	      length += myChildren.size();
	      path = myPath;
	      
	    } else {
	      
	      // Create a new working set by adding the remainder of the old
	      // working set to the end of the current children list

	      myChildren.insert(myChildren.end(),remainingIds.begin(),remainingIds.end());

	      buildPencils(grid,pencils,idsOut,myChildren,dimension,myPath);
	      
	    };
	    
	  };	  
	};      
      };
    
    } else {

      // Add unrefined cells to the pencil directly

      idsOut.push_back(id);
      
    }; // closes if(isRefined)

    // Move to the next cell
    i++;
    
  }; // closes loop over ids

  pencils.addPencil(idsOut,0.0,0.0);
  return pencils;
  
} // closes function

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
  
  std::cout << "Grid size at 0 refinement is " << xDim << " x " << yDim << " x " << zDim << std::endl;
  for (const auto& cell: cells) {
    // std::cout << "Data of cell " << cell.id << " is stored at " << cell.data << std::endl;
    // Collect a list of cell ids.
    ids.push_back(cell.id);
    
    // Add parent cells of refined cells to the list of cell ids.
    CellID parent = grid.mapping.get_parent(cell.id);
    if (parent > 0 &&
	!(std::find(ids.begin(), ids.end(), parent) != ids.end())) {
      ids.push_back(parent);
      std::cout << "Cell " << parent << " at refinement level " << grid.get_refinement_level(parent) << " has been refined into ";
      for (const auto& child: grid.mapping.get_all_children(parent)) {
	std::cout << child << " ";
      }
      std::cout << "\n";
    }
  }

  sort(ids.begin(),ids.end());
  
  uint ibeg = 0;
  uint iend = 0;

  uint dimension = 0;
  vector<uint> dims;
  
  switch( dimension ) {
  case 0 : {
    dims = {zDim,yDim,xDim};
  }
    break;
  case 1 : {
    dims = {zDim,xDim,yDim};
  }
    break;
  case 2 : {
    dims = {yDim,xDim,zDim};
  }
    break;
  default : {
    dims = {0,0,0};
  }
  };

  map <CellID,CellID> mapping;
  sortIds< CellID, dccrg::Grid_Length::type >(dimension, grid_size, ids, mapping);
  
  list < vector < CellID >> unrefinedPencils;
  std::cout << "The unrefined pencils along dimension " << dimension << " are:\n";
  for (uint i = 0; i < dims[0]; i++) {
    for (uint j = 0; j < dims[1]; j++) {
      vector <CellID> unrefinedIds;
      ibeg = 1 + i * dims[2] * dims[1] + j * dims[2];
      iend = 1 + i * dims[2] * dims[1] + (j + 1) * dims[2];
      for (uint k = ibeg; k < iend; k++) {
	std::cout << mapping[k] << " ";
	unrefinedIds.push_back(mapping[k]);
	//unrefinedIds.push_back(ids[k]);
      }
      unrefinedPencils.push_back(unrefinedIds);
      std::cout << "\n";
    }
  }

  ibeg = 0;
  iend = 0;
  
  setOfPencils pencilInitial;
  vector<CellID> idsInitial;
  vector<pair<bool,bool>> path;
  
  setOfPencils pencils;
  for ( auto &unrefinedIds : unrefinedPencils ) {
    pencils = buildPencils(grid, pencilInitial, idsInitial, unrefinedIds, dimension, path);
  }

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
