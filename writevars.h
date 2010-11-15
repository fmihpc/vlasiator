#ifndef WRITEVARS_H
#define WRITEVARS_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include "definitions.h"

// Variable saving needs to be completely rewritten

void openOutputFile(const std::string& fileName);
void closeOutputFile(const std::string& fileName);

void writeVelocityBlockGrid(const std::string& gridName,cuint& BLOCKS,Real* blockParams);
void writeVelocityBlockScalar(const std::string& varName,const std::string& gridName,cuint& BLOCKS,Real* array);

#endif
