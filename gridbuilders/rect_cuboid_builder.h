#ifndef RECT_CUBOID_BUILDER_H
#define RECT_CUBOID_BUILDER_H

#include <vector>
#include <map>
#include "../gridbuilder.h"

// This GridBuilder creates a rectangular cuboid 
// (see, e.g. Wikipedia). 

class RectCuboidBuilder: public GridBuilder {
 public:
   RectCuboidBuilder();
   ~RectCuboidBuilder();
   
   bool finalize();
   bool getNextCell(lluint& cellID,Real* coords,Real* sizes,std::vector<std::pair<lluint,uchar> >& nbrs);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(lluint& N_cells);
   bool initialize();
   
 protected:
   bool initialized;
   Real dx;
   Real dy;
   Real dz;
   Real xmin;
   Real xmax;
   Real ymin;
   Real ymax;
   Real zmin;
   Real zmax;
   uint xsize;
   uint ysize;
   uint zsize;
   
   uint calcIndex(cuint& i,cuint& j,cuint& k);
};

#endif
