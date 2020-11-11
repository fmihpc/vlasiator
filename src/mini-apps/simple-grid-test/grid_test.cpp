#include "dccrg.hpp"
#include "mpi.h"
#include <iostream>
#include "fstream"

struct grid_data {

  int value = 0;

  std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
  {
    return std::make_tuple(this, 0, MPI_BYTE);
  }
    
};

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

  const int xDim = 16;
  const int yDim = 16;
  const int zDim = 16;
  const std::array<uint64_t, 3> grid_size = {{xDim,yDim,zDim}};
  
  grid.initialize(grid_size, comm, "RANDOM", 1);

  grid.balance_load();

  bool doRefine = true;
  const std::array<uint,4> refinementIds = {{4,22,29,4104}};
  if(doRefine) {
    for(uint i = 0; i < refinementIds.size(); i++) {
      if(refinementIds[i] > 0) {
	grid.refine_completely(refinementIds[i]);
	grid.stop_refining();
      }
    }
  }

  grid.balance_load();

  //  std::vector<uint64_t> cells = grid.get_cells()

  auto cells = grid.cells;
  sort(cells.begin(), cells.end());

  std::cout << "Grid size at 0 refinement is " << xDim << " x " << yDim << " x " << zDim << std::endl;
  for (const auto& cell: cells) {
    std::cout << "Data of cell " << cell.id << " is stored at " << cell.data << std::endl;
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
