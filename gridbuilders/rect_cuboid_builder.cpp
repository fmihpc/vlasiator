#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "../definitions.h"
#include "../parameters.h"
#include "../utility.h"
#include "rect_cuboid_builder.h"

using namespace std;

RectCuboidBuilder::RectCuboidBuilder(): GridBuilder(),initialized(false) { }

RectCuboidBuilder::~RectCuboidBuilder() { }

bool RectCuboidBuilder::finalize() {return true;}

bool RectCuboidBuilder::getNextCell(uint maxNbrs,lluint& cellID,Real* coords,Real* sizes,lluint* nbrs,uchar* nbrTypes) {
   
}

bool RectCuboidBuilder::getParameter(const std::string& parameterName,std::string& value) {
   typedef Parameters P;
   value = "";
   map<string,string>::const_iterator it = parameters.find(parameterName);
   if (it == parameters.end()) return false;
   value = it->second;
   return true;
}

bool RectCuboidBuilder::getTotalNumberOfCells(lluint& N_cells) {
   
   return false;
}

bool RectCuboidBuilder::initialize() {
   typedef Parameters P;
   bool initialized = true;
   parameters["xmin"]   = floatToString(P::xmin);
   parameters["xmax"]   = floatToString(P::xmax);
   parameters["ymin"]   = floatToString(P::ymin);
   parameters["ymax"]   = floatToString(P::ymax);
   parameters["zmin"]   = floatToString(P::zmin);
   parameters["zmax"]   = floatToString(P::zmax);
   parameters["dx"]     = floatToString(P::dx_ini);
   parameters["dy"]     = floatToString(P::dy_ini);
   parameters["dz"]     = floatToString(P::dz_ini);
   parameters["xcells"] = floatToString(P::xcells_ini);
   parameters["ycells"] = floatToString(P::ycells_ini);
   parameters["zcells"] = floatToString(P::zcells_ini);
   return initialized;
}

// Register RectCuboidBuilder:
GridBuilder* newBuilder() {return new RectCuboidBuilder;}

class Dummy {
 public:
   Dummy() {
      GridBuilderFactory::registerBuilder(newBuilder);
   }
   ~Dummy() { }
};

static Dummy dummy;

