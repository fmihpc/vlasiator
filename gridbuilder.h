#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include "definitions.h"

class GridBuilder {
 public:
   GridBuilder();
   virtual ~GridBuilder();

   virtual bool finalize() =0;
   virtual bool getNextCell(uint maxNbrs,lluint& cellID,Real* coords,Real* sizes,lluint* nbrs,uchar* nbrTypes) = 0;
   virtual bool getParameter(const std::string& parameterName,std::string& value) =0;
   virtual bool getTotalNumberOfCells(lluint& N_cells) =0;
   virtual bool initialize() =0;
   
 protected:

};

typedef GridBuilder* (*GridBuilderCreator)();

class GridBuilderFactory {
 public:
   GridBuilderFactory();
   ~GridBuilderFactory();
   
   static GridBuilder* createBuilder();
   static bool deleteBuilder(GridBuilder*& gb);
   static bool registerBuilder(GridBuilderCreator gbc);
   
 private:
   static GridBuilderCreator gbc;
};

#endif