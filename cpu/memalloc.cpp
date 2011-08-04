#include <cstdlib>
#include <iostream>

#include "../memalloc.h"

using namespace std;

/* Allocate a float array of requested size. The functionality of this function depends on CUDA
 * preprocessor flag. If memalloc.cpp is compiled with -DCUDA flag, the memory is allocated as 
 * page-locked (pinned) as required by GPU Vlasov solver. If -DCUDA is not defined, the array 
 * is allocated with new operator. With -DCUDA flag this function will abort the simulation 
 * if the requested array could not be allocated.
 * @param ptr Pointer to the array.
 * @param elements Requested number of elements in the array.
 */
void allocateArray(float** ptr,const size_t& elements) {
   try {
	   *ptr = new float[elements];
   }
   catch (exception& e) {
      cerr << __FILE__ << ":" << __LINE__
         << "Couldn't allocate memory: " << e.what()
         << endl;
      abort();
   }

}

/** Free an array which has been previously allocated with allocateArray. The functionality 
 * of this function depends on CUDA preprocessor flag. If memalloc.cpp is compiled with 
 * -DCUDA flag, the memory is freed with GPU driver function. Otherwise the memory is 
 * freed with delete operator.
 * @param ptr Pointer to the array to be deallocated.
 */
void freeArray(float* ptr) {
   delete ptr;
}





