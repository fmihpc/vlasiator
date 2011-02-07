#ifndef RECT_CUBOID_BUILDER_H
#define RECT_CUBOID_BUILDER_H

#include <map>
#include "../gridbuilder.h"

// This GridBuilder creates a rectangular cuboid 
// (see, e.g. Wikipedia). 

class RectCuboidBuilder: public GridBuilder {
 public:
   RectCuboidBuilder();
   ~RectCuboidBuilder();
   
   bool finalize();
   bool getNextCell(uint maxNbrs,lluint& cellID,Real* coords,Real* sizes,lluint* nbrs,uchar* nbrTypes);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(lluint& N_cells);
   bool initialize();
   
 protected:
   bool initialized;
   std::map<std::string,std::string> parameters;
};

#endif
