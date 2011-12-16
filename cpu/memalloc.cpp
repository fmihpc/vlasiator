/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>

#include "../definitions.h"
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
void allocateArray(Real** ptr,const size_t& elements) {
   try {
	   *ptr = new Real[elements];
   }
   catch (exception& e) {
      cerr << __FILE__ << ":" << __LINE__
         << " Couldn't allocate memory: " << e.what()
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
void freeArray(Real* ptr) {
   delete ptr;
}





