#ifndef ASCIIREADER_H
#define ASCIIREADER_H

#include <vector>
#include "../gridbuilder.h"

class AsciiBuilder: public GridBuilder {
 public:
   AsciiBuilder();
   ~AsciiBuilder();
   
   bool finalize();
   bool getNextCell(lluint& cellID,Real* coords,Real* sizes,std::vector<std::pair<lluint,uchar> >& nbrs);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(lluint& N_cells);
   bool initialize();
   
 protected:
   bool initialized;
   uint N_cells;
   uint N_dimensions;
};

#endif
