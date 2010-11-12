#ifndef GRID_H
#define GRID_H

#include "definitions.h"

class Grid {
 public:
   Grid();
   ~Grid();

   uint getFreeMemory(cuint& BLOCKS);
   bool addReference(cuint& INDEX);
   bool removeReference(cuint& INDEX);

   uint* getNbrsVel() const {return nbrsVel;}
   Real* getBlockParams() const {return blockParams;}
   Real* getAvgs() const {return avgs;}
   Real* getFx() const {return fx;}
   Real* getFy() const {return fy;}
   Real* getFz() const {return fz;}
   Real* getD1x() const {return d1x;}
   Real* getD1y() const {return d1y;}
   Real* getD1z() const {return d1z;}
   Real* getD2x() const {return d2x;}
   Real* getD2y() const {return d2y;}
   Real* getD2z() const {return d2z;}
   
   void printReferences();
   
 private:
   Grid(const Grid& g);
   
   uint* nbrsVel; /**< Pointer to array which is used to store velocity block neighbour lists in CPU memory.*/
   Real* blockParams; /**< Pointer to array which is used to store velocity block parameters in CPU memory.*/
   Real* avgs; /**< Pointer to array which is used to store velocity blocks in CPU memory.*/
   Real* fx; /**< Pointer to array containing fluxes to x-direction.*/
   Real* fy; /**< Pointer to array containing fluxes to y-direction.*/
   Real* fz; /**< Pointer to array containing fluxes to z-direction.*/
   Real* d1x; /**< Pointer to array containing 1st derivatives to x-direction.*/
   Real* d1y; /**< Pointer to array containing 1st derivatives to y-direction.*/
   Real* d1z; /**< Pointer to array containing 1st derivatives to z-direction.*/
   Real* d2x; /**< Pointer to array containing 2nd derivatives to x-direction.*/
   Real* d2y; /**< Pointer to array containing 2nd derivatives to y-direction.*/
   Real* d2z; /**< Pointer to array containing 2nd derivatives to z-direction.*/
};

#endif
