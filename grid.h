#ifndef GRID_H
#define GRID_H

#include <vector>
#include "definitions.h"
#include "parameters.h"
#include "cell_spatial.h"

class Grid {
 public:
   Grid();
   ~Grid();

   /** Query if the Grid has been properly initialized.
    * @return If true, the Grid is ready for use. If false, an error has occurred and
    * the simulation should be aborted.
    */
   bool isInitialized() const {return initialized;}
   
   uint getSpatialCell(const uint& velBlocks);

   /** Get a pointer to the array containing velocity grid blocks.
    * @return A pointer to blockArray.
    */
   real* getBlockArray() const {return blockArray;}
   /** Get a pointer to the array containing parameters of the velocity grid blocks.
    * @return A pointer to blockParams.
    */
   real* getBlockParams() const {return blockParams;}
   /** Get a pointer to the array containing parameters of the spatial cells.
    * @return A pointer to cellParams.
    */
   real* getCellParams() const {return cellParams;}
   /** Get a pointer to the array containing neighbour lists for spatial cells.
    * @return A pointer to nbrsSpa.
    */
   uint* getNbrsSpa() const {return nbrsSpa;}
   /** Get a pointer to the array containing the neighbour lists for the velocity blocks.
    * @return A pointer to nbrsVel.
    */
   uint* getNbrsVel() const {return nbrsVel;}

   /** Get the number of spatial cells currently stored in the Grid.
    * @return Number of stored cells.
    */
   uint size() const {return nextSpaCell;}
   
   /** Return a pointer to the requested spatial cell.
    * @param index The index of the cell in blockArray.
    * @return A pointer to the requested cell, or a NULL pointer if the cell does not exist.
    */
   SpatialCell* operator[](cuint& index) {return &(cells[index]);}
   
 private:
   /** Private copy constructor to prevent the creation of new copies of Grid.
    * @param g The Grid to be copied.
    */
   Grid(const Grid& g);
   
   bool initialized; /**< If true, Grid has been initialized correctly and is ready for use.*/
   uint nextVelBlock; /**< The index of the next free velocity grid block in blockArray.*/
   uint nextSpaCell; /**< The index of the next free spatial cell in array cells.*/
   
   uint* nbrsSpa; /**< Pointer to array which is used to store spatial cell neighbour lists in CPU memory.*/
   uint* nbrsVel; /**< Pointer to array which is used to store velocity block neighbour lists in CPU memory.*/
   real* blockParams; /**< Pointer to array which is used to store velocity block parameters in CPU memory.*/
   real* cellParams; /**< Pointer to array which is used to store spatial cell parameters in CPU memory.*/
   real* blockArray; /**< Pointer to array which is used to store velocity blocks in CPU memory.*/

   SpatialCell* cells; /**< Pointer to array which is used to store the spatial cells.*/   
};

#endif
