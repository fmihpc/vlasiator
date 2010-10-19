#ifndef CUDAFUNCS_H
#define CUDAFUNCS_H

#include <iostream>

#include "definitions.h"

void* gpuCreateArray(const std::string& name,const size_t& bytes);
void gpuDeleteArray(const std::string& name,void* ptr);
void gpuCopyArray(const std::string& name,const size_t& bytes,void* cpuPtr,void* gpuPtr,const bool& cpuToGpu);

bool deviceCreateArray(real*& arrptr,const size_t& bytes);
bool deviceCreateArray(uint*& arrptr,const size_t& bytes);
bool deviceDeleteArray(real*& arrptr);
bool deviceDeleteArray(uint*& arrptr);

#endif
