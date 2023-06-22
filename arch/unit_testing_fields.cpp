
/* This is a unit testing suite for the Vlasiator unified CPU-GPU for loop interface.
 * Compile for device execution with
 *  `USE_CUDA=1 VLASIATOR_ARCH=mahti_cuda make -j12 unit_testing_fields`
 * and for host execution with
 *  `USE_CUDA=0 VLASIATOR_ARCH=mahti_cuda make -j12 unit_testing_fields`
 */

/* Included standard headers */
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <tuple> 
#include <vector>
#include <fsgrid.hpp>

#define ARCH_MAIN 1

/* Include the tested architecture-specific header */
#include "arch_device_api.h"
#include "../mpiconversion.h"
#include "../common.h"

/* Host execution of min() and max() require using std namespace */
using namespace std;

/* Auxiliary function for result evaluation and printing */
void result_eval(std::tuple<bool, double, double> res, const uint test_id){
  std::string success = std::get<0>(res) == true ? "PASSED" : "FAILED";
  printf("Test %d %s - Arch: %9.2f µs, Host: %9.2f µs\n", test_id, success.c_str(), std::get<1>(res), std::get<2>(res));
}

/* The test functions are all named as `test`, and only differentiated 
 * by their ascending template id number `I`. This allows executing tests
 * nicely using an array of function pointers, and not calling 
 * each test by a separate name. New tests can be added by inserting 
 * a new function named `test` with the `std::enable_if<I == value, ...` 
 * construct with the next unused value for `I`. 
 */


// Define a test function that compares the performance of two different implementations of a loop that sets the value of a field in a 3D grid.
template<uint I>
typename std::enable_if<I == 0, std::tuple<bool, double, double>>::type test(){

  // Define the size of the 3D grid
  int32_t gridDims[3] = {100, 100, 100};
  uint gridDims3 = gridDims[0] * gridDims[1] * gridDims[2];

  // Initialize MPI and grid coupling information
  MPI_Comm comm = MPI_COMM_WORLD;
  FsGridCouplingInformation gridCoupling;

  // Set the periodicity of the grid
  std::array<bool,3> periodicity{true, true, true};

  // Create the 3D grid with a field of type fsgrids::technical
  FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> technicalGrid(gridDims, comm, periodicity,gridCoupling); 
  FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> technicalGrid2(gridDims, comm, periodicity,gridCoupling); 

  // Create a buffer object that provides a convenient interface for accessing the grid data on the device
  arch::buf<fsgrids::technical> dataBuffer(&technicalGrid.get(0), gridDims3 * sizeof(fsgrids::technical));  
  arch::buf<fsgrids::technical> dataBuffer2(&technicalGrid2.get(0), gridDims3 * sizeof(fsgrids::technical));  

  // Execute the loop in parallel on the device using CUDA
  clock_t arch_start = clock();
  arch::parallel_for({(uint)gridDims3}, ARCH_LOOP_LAMBDA(int i) {
    dataBuffer[i].fsGridRank=2; 
  }); 
  double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);

  // sync device to host
  dataBuffer.syncHostData();

  // Execute the loop on the host
  clock_t host_start = clock();
  for (uint k = 0; k < gridDims3; ++k){
    dataBuffer2[k].fsGridRank=2;
  }
  double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC); 

  // Check whether the underlying data arrays on host have the same data
  bool success = true;
  if (dataBuffer2[10].fsGridRank != 2) { // TODO: check pointers!!!
    cout << "dataBuffer2 value is not 2!" << endl; 
    success = false;
  }
  if (technicalGrid2.get(10).fsGridRank != 2) { // TODO: check pointers!!!
    cout << "technicalgrid2 value is not 2!" << endl; 
    success = false;
  }
  if (technicalGrid.get(10).fsGridRank != 2) { // TODO: check pointers!!!
    cout << "technicalgrid value is not 2!" << endl; 
    success = false;
  }

  // Return a tuple containing the success status and the execution times for the device and host implementations
  return std::make_tuple(success, arch_time, host_time);
}
// This test function compares the performance of two 
// different implementations of a loop that sets the value of a field in a 3D grid.
template<uint I>
typename std::enable_if<I == 1, std::tuple<bool, double, double>>::type test(){

  // Define the size of the 3D grid
  // int32_t gridDims[3] = {100, 100, 100};
  int32_t fsGridDimensions[3] = {100, 100, 100};

  // Initialize MPI and grid coupling information
  MPI_Comm comm = MPI_COMM_WORLD;
  FsGridCouplingInformation gridCoupling;

  // Set the periodicity of the grid
  std::array<bool,3> periodicity{true, true, true};

  // Create the 3D grid with a field of type std::array<Real, fsgrids::bfield::N_BFIELD>
  FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> perBGrid(fsGridDimensions, comm, periodicity,gridCoupling); 
  FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> perBGrid2(fsGridDimensions, comm, periodicity,gridCoupling); 
  int32_t *gridDims = perBGrid.getStorageSize();

  // Create a buffer object that provides a convenient interface for accessing the grid data on the device
  arch::buf<FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> perBGridBuf(&perBGrid);
  arch::buf<FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> perBGridBuf2(&perBGrid2);

  // Execute the loop in parallel on the device using CUDA
  clock_t arch_start = clock();
  arch::parallel_for({(uint)gridDims[0], (uint)gridDims[1], (uint)gridDims[2]}, ARCH_LOOP_LAMBDA(int i, int j, int k) {
    perBGridBuf.get(i, j, k).at(fsgrids::bfield::PERBX) = 2;
  });  
  double arch_time = (double)((clock() - arch_start) * 1e6 / CLOCKS_PER_SEC);
  perBGridBuf.syncHostData();

  // Execute the loop on the host
  clock_t host_start = clock();
  for (uint k = 0; k < gridDims[2]; ++k){
    for (uint j = 0; j < gridDims[1]; ++j){
      for (uint i = 0; i < gridDims[0]; ++i) {
        perBGridBuf2.get(i, j, k).at(fsgrids::bfield::PERBX) = 2;
      }
    } 
  }
  double host_time = (double)((clock() - host_start) * 1e6 / CLOCKS_PER_SEC); 

  // Check whether the test was successful
  bool success = true;
  if (perBGridBuf.get(10,10,10)[fsgrids::bfield::PERBX] != 2)
    success = false; 
  if (perBGridBuf2.get(10,10,10)[fsgrids::bfield::PERBX] != 2)
    success = false; 
  return std::make_tuple(success, arch_time, host_time); 
}


/* Instantiate each test function by recursively calling the
 * driver function in a descending order beginning from `N - 1`
 */
template<uint N, uint I>
struct test_instatiator {
  static void driver(std::tuple<bool, double, double>(*fptr_test[N])()){
    fptr_test[I - 1] = &test<I - 1>;
    test_instatiator<N, I - 1>::driver(fptr_test);
  }
};

/* Specialization for the instantiation end condition `I = 0` */
template<uint N>
struct test_instatiator<N, 0> {
  static void driver(std::tuple<bool, double, double>(*fptr_test[N])()){}
};

/* The main function */
int main(int argn,char* args[]) {

  // init MPI, required for grids initialization
  int required=MPI_THREAD_FUNNELED;
  int provided;
  MPI_Init_thread(&argn,&args,required,&provided);
    
  /* Specify the number of tests and set function pointers */
  constexpr uint n_tests = 2;
  std::tuple<bool, double, double>(*fptr_test[n_tests])();
  test_instatiator<n_tests, n_tests>::driver(fptr_test);

  /* Indicate for what backend option the test suite is compiled */
  #ifdef USE_CUDA
    printf("Run tests for Arch = CUDA (USE_CUDA defined)\n");
  #else
    printf("Run tests for Arch = HOST (USE_CUDA not defined)\n");
  #endif

  
  /* Evaluate all test cases using the array of function pointers */
  for(uint i = 0; i < n_tests; i++) {
    printf("Running Test %d\n", i);
    fflush(stdout);
    result_eval(fptr_test[i](), i);
  }

  /* Finalize MPI */
  MPI_Finalize(); 

  return 0;
}
