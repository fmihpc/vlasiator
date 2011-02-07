#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "asciibuilder.h"
#include "asciireaders.h"

using namespace std;

static fstream filein;

AsciiBuilder::AsciiBuilder(): GridBuilder() {initialized = false;}

AsciiBuilder::~AsciiBuilder() {filein.close();}

bool AsciiBuilder::finalize() {filein.close(); return true;}

bool AsciiBuilder::getNextCell(uint maxNbrs,lluint& cellID,Real* coords,Real* sizes,lluint* nbrs,uchar* nbrTypes) { 
   if (initialized == false) return false;

   // Read next cell:
   vector<string> line;
   if (readline(filein,line) == false) return false;
   if (line.size() < 1 + 2*N_dimensions) return false;
   
   cellID = atoi(line[0].c_str());
   unsigned int offset = 1;
   for (int i=0; i<N_dimensions; ++i) {coords[i] = atof(line[i+offset].c_str());}
   offset += N_dimensions;
   for (int i=0; i<N_dimensions; ++i) {sizes[i] = atof(line[i+offset].c_str());}
   offset += N_dimensions;

   unsigned int i=0;
   while (offset < line.size()) {
      const unsigned long long int nbrID = atoi(line[offset].substr(0,line[offset].find(':')).c_str());
      const unsigned int nbrType = atoi(line[offset].substr(line[offset].find(':')+1,line[offset].size()-line[offset].find(':')-1).c_str());
      nbrs[i] = nbrID;
      nbrTypes[i] = nbrType;
      ++offset;
      ++i;
   }
   return true;
}

bool AsciiBuilder::getParameter(const std::string& parameterName,std::string& value) {
   return false;
}

bool AsciiBuilder::getTotalNumberOfCells(lluint& N_cells) {
   N_cells = this->N_cells;
   return initialized;
}

bool AsciiBuilder::initialize() {
   // Open input file. Get the filename from somewhere.
   filein.open("vlasov.grid");
   if (filein.good() == false) return false;
   initialized = true;
   cerr << "File opened" << endl;
   vector<string> line;
   // Parse grid dimensionality:
   if (readline(filein,line) == false) {initialized = false; return false;}
   if (line.size() != 2) {initialized = false; return false;}
   if (line[0] != "DIMENSIONS") {initialized = false; return false;}
   N_dimensions = atoi(line[1].c_str());
   line.clear();
   cerr << "Dimensions = " << N_dimensions << endl;
   // Parse the number of cells:
   if (readline(filein,line) == false) {initialized = false; return false;}
   if (line.size() != 2) {initialized = false; return false;}
   if (line[0] != "CELLS") {initialized = false; return false;}
   N_cells = atoi(line[1].c_str());
   cerr << "Cells = " << N_cells << endl;
   
   if (initialized == true) cerr << "init OK" << endl;
   else cerr << "init FAIL" << endl;
   return initialized;
}

// Register AsciiBuilder:

GridBuilder* newAsciiBuilder() {return new AsciiBuilder;}

class Dummy {
 public:
   Dummy() {
      GridBuilderFactory::registerBuilder(newAsciiBuilder);
   }
   ~Dummy() { }
};

static Dummy dummy;
