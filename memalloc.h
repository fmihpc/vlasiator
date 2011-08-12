#ifndef MEMALLOC_H
#define MEMALLOC_H

#include <cstdlib>
#include <iostream>
#include "definitions.h"

void allocateArray(Real** ptr,const size_t& size);
void freeArray(Real* ptr);

#endif
