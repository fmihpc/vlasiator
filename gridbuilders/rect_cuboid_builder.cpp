#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include "../definitions.h"
#include "../parameters.h"
#include "rect_cuboid_builder.h"

using namespace std;

RectCuboidBuilder::RectCuboidBuilder(): GridBuilder(),initialized(false) { }

RectCuboidBuilder::~RectCuboidBuilder() { }

uint RectCuboidBuilder::calcIndex(cuint& i,cuint& j,cuint& k) {return k*ysize*xsize + j*xsize + i;}

bool RectCuboidBuilder::finalize() {return true;}

bool RectCuboidBuilder::getNextCell(lluint& cellID,Real* coords,Real* sizes,std::vector<std::pair<lluint,uchar> >& nbrs) {
   // Create a counter for counting the number of returned cells.
   static long unsigned int counter = 0;
   if (counter == xsize*ysize*zsize) {
      return false;
   }
   
   cellID = counter;
   cuint k = counter / (xsize*ysize);
   cuint j = (counter - k*xsize*ysize)/xsize;
   cuint i = counter - k*xsize*ysize - j*xsize;
   
   coords[0] = xmin + i*dx;
   coords[1] = ymin + j*dy;
   coords[2] = zmin + k*dz;
   sizes[0] = dx;
   sizes[1] = dy;
   sizes[2] = dz;
  
   // Push existing neighbours to nbrs:
   nbrs.clear();
   if (i != 0)       nbrs.push_back(make_pair(calcIndex(i-1,j,k), 0));
   if (i != xsize-1) nbrs.push_back(make_pair(calcIndex(i+1,j,k), 8));
   if (j != 0)       nbrs.push_back(make_pair(calcIndex(i,j-1,k),16));
   if (j != ysize-1) nbrs.push_back(make_pair(calcIndex(i,j+1,k),24));
   if (k != 0)       nbrs.push_back(make_pair(calcIndex(i,j,k-1),32));
   if (k != zsize-1) nbrs.push_back(make_pair(calcIndex(i,j,k+1),40));

   ++counter;
   return true;
}

bool RectCuboidBuilder::getParameter(const std::string& parameterName,std::string& value) {
   return false;
}

bool RectCuboidBuilder::getTotalNumberOfCells(lluint& N_cells) {
   N_cells = xsize*ysize*zsize;
   return true;
}

bool RectCuboidBuilder::initialize() {
   typedef Parameters P;
   bool initialized = true;

   // Define the required parameters and query their values from Parameters:
   P::add("gridbuilder.x_min","minimum value of the x-coordinate",xmin,numeric_limits<Real>::max());
   P::add("gridbuilder.x_max","minimum value of the x-coordinate",xmax,numeric_limits<Real>::max());
   P::add("gridbuilder.y_min","minimum value of the y-coordinate",ymin,numeric_limits<Real>::max());
   P::add("gridbuilder.y_max","minimum value of the y-coordinate",ymax,numeric_limits<Real>::max());
   P::add("gridbuilder.z_min","minimum value of the z-coordinate",zmin,numeric_limits<Real>::max());
   P::add("gridbuilder.z_max","minimum value of the z-coordinate",zmax,numeric_limits<Real>::max());
   P::add("gridbuilder.x_length","number of cells in x-direction in initial grid",xsize,numeric_limits<uint>::max());
   P::add("gridbuilder.y_length","number of cells in y-direction in initial grid",ysize,numeric_limits<uint>::max());
   P::add("gridbuilder.z_length","number of cells in z-direction in initial grid",zsize,numeric_limits<uint>::max());
   P::parse();
   
   // Check that we received sane values:
   if (xmin == numeric_limits<Real>::max()) initialized = false;
   if (xmax == numeric_limits<Real>::max()) initialized = false;
   if (ymin == numeric_limits<Real>::max()) initialized = false;
   if (ymax == numeric_limits<Real>::max()) initialized = false;
   if (zmin == numeric_limits<Real>::max()) initialized = false;
   if (zmax == numeric_limits<Real>::max()) initialized = false;
   if (xmax < xmin || (ymax < ymin || zmax < zmin)) initialized = false;
   if (xsize == numeric_limits<uint>::max()) initialized = false;
   if (ysize == numeric_limits<uint>::max()) initialized = false;
   if (zsize == numeric_limits<uint>::max()) initialized = false;
   
   dx = (xmax-xmin)/xsize;
   dy = (ymax-ymin)/ysize;
   dz = (zmax-zmin)/zsize;
   
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

