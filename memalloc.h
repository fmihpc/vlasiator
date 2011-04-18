#ifndef MEMALLOC_H
#define MEMALLOC_H

#include <cstdlib>
#include <iostream>

void allocateArray(float** ptr,const size_t& size);
void freeArray(float* ptr);

#endif
