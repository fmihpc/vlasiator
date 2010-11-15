#include <cstdlib>
#include <iostream>
#include <sstream>

#include "logger.h"
#include "writevars.h"
#include "silowriter.h"

using namespace std;

extern Logger logger;

void openOutputFile(const std::string& fileName) {
   openOutputFile(fileName.c_str(),"test");
}

void closeOutputFile(const std::string& fileName) {
   closeOutputFile();
}

void writeVelocityBlockGrid(const string& gridName,cuint& BLOCKS,Real* blockParams) {
   //cout << "Writing " << BLOCKS << " blocks to silo" << endl;
   writeVelocityBlockGrid3D(gridName,BLOCKS,blockParams);
}

void writeVelocityBlockScalar(const string& varName,const string& gridName,cuint& BLOCKS,Real* array) {
   //cout << "Writing data from address " << array << endl;
   writeVelocityBlockGridScalar3D(varName,gridName,BLOCKS,array);
}

