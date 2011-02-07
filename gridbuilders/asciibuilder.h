#ifndef ASCIIREADER_H
#define ASCIIREADER_H

#include "../gridbuilder.h"

class AsciiBuilder: public GridBuilder {
 public:
   AsciiBuilder();
   ~AsciiBuilder();
   
   bool finalize();
   bool getNextCell(uint maxNbrs,lluint& cellID,Real* coords,Real* sizes,lluint* nbrs,uchar* nbrTypes);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(lluint& N_cells);
   bool initialize();
   
 protected:
   bool initialized;
   uint N_cells;
   uint N_dimensions;
};

#endif
