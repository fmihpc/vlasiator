#ifndef DEVICEGRID_H
#define DEVICEGRID_H

#include <vector>
#include "common.h"
#include "definitions.h"
#include "cell_spatial.h"
#include "parameters.h"

/** This class manages the velocity blocks in the device memory.
 */
class DeviceGrid {
 public:
   /** Create a new DeviceGrid and allocate memory for velocity blocks and related book-keeping data in 
    * device memory.
    */
   DeviceGrid();
   
   //DeviceGrid(const std::vector<uint>& realElemSizes,const std::vector<uint>& uintElemSizes);
   
   /** Deallocates the memory that has been previously allocated by a call to DeviceGrid().*/
   ~DeviceGrid();
   
   /** Query if the DeviceGrid is properly initialized. 
    * @return If true, all memory was allocated successfully and DeviceGrid is ready for use.
    * If false, an error has occurred and the simulation should be aborted.
    */
   bool isInitialized() const {return initialized;}
   
   void initSpatialCell(SpatialCell& cell);

   real* getAvgnbrx() const {return avgnbrx;}
   real* getAvgnbry() const {return avgnbry;}
   real* getAvgnbrz() const {return avgnbrz;}
   real* getBlockArray() const {return blockArray;}
   real* getBlockParams() const {return blockParams;}
   real* getCellParams() const {return cellParams;}
   real* getD1x() const {return d1x;}
   real* getD1y() const {return d1y;}
   real* getD1z() const {return d1z;}
   real* getD2x() const {return d2x;}
   real* getD2y() const {return d2y;}
   real* getD2z() const {return d2z;}
   real* getD1xnbr() const {return d1xnbr;}
   real* getD1ynbr() const {return d1ynbr;}
   real* getD1znbr() const {return d1znbr;}
   real* getD2xnbr() const {return d2xnbr;}
   real* getD2ynbr() const {return d2ynbr;}
   real* getD2znbr() const {return d2znbr;}
   real* getFx() const {return fx;}
   real* getFy() const {return fy;}
   real* getFz() const {return fz;}
   real* getFxnbr() const {return fxnbr;}
   real* getFynbr() const {return fynbr;}
   real* getFznbr() const {return fznbr;}
   uint* getNbrsSpa() const {return nbrsVel;}
   uint* getNbrsVel() const {return nbrsVel;}

   uint sizeBlockArray() const {return MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real);}
   uint sizeBlockParams() const {return MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(real);}
   uint sizeD1x() const {return MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real);}
   uint sizeD1y() const {return MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real);}
   uint sizeD1z() const {return MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real);}
   uint sizeCellParams() const {return MAX_SPA_CELLS*SIZE_CELLPARAMS*sizeof(real);}
   uint sizeNbrsVel() const {return MAX_VEL_BLOCKS*SIZE_NBRS_VEL*sizeof(uint);}
   
 private:
   bool initialized; /**< If true, DeviceGrid was initialized successfully.*/
   uint nextInner; /**< A pointer to a position in blockArray where the next inner velocity block can be written. */
   uint nextBoundary; /**< A pointer to a position in blockArray where the next boundary velocity block can be written. */
   
   real* cellParams; /**< A pointer to array containing spatial cell parameters in device memory.*/
   
   real* blockArray; /**< A pointer to the array containing velocity blocks in device memory. */
   real* blockParams; /**< A pointer to an array containing velocity block parameter data in device memory.*/
   uint* nbrsVel; /**< A pointer to an array containing neighbour lists in velocity space.*/

   real* d1x; /**< Pointer to array containing 1st derivatives to x-direction for each velocity grid block.*/
   real* d1y; /**< Pointer to array containing 1st derivatives to y-direction for each velocity grid block.*/
   real* d1z; /**< Pointer to array containing 1st derivatives to z-direction for each velocity grid block.*/
   real* d2x; /**< Pointer to array containing 2nd derivatives to x-direction for each velocity grid block.*/
   real* d2y; /**< Pointer to array containing 2nd derivatives to y-direction for each velocity grid block.*/
   real* d2z; /**< Pointer to array containing 2nd derivatives to z-direction for each velocity grid block.*/
   real* fx; /**< Pointer to array containing the fluxes to x-direction for each velocity grid block.*/
   real* fy; /**< Pointer to array containing the fluxes to y-direction for each velocity grid block.*/
   real* fz; /**< Pointer to array containing the fluxes to z-direction for each velocity grid block.*/
   real* avgnbrx; /**< Pointer to array containing the ghost cells for volume averages in x-direction.*/
   real* avgnbry; /**< Pointer to array containing the ghost cells for volume averages in y-direction.*/
   real* avgnbrz; /**< Pointer to array containing the ghost cells for volume averages in z-direction.*/
   real* d1xnbr; /**< Pointer to array containing ghost cells for 1st derivatives in x-direction. */
   real* d1ynbr; /**< Pointer to array containing ghost cells for 1st derivatives in y-direction. */
   real* d1znbr; /**< Pointer to array containing ghost cells for 1st derivatives in z-direction. */
   real* d2xnbr; /**< Pointer to array containing ghost cells for 2nd derivatives in x-direction. */
   real* d2ynbr; /**< Pointer to array containing ghost cells for 2nd derivatives in y-direction. */
   real* d2znbr; /**< Pointer to array containing ghost cells for 2nd derivatives in z-direction. */
   real* fxnbr; /**< Pointer to array containing ghost cells for fluxes in x-direction.*/
   real* fynbr; /**< Pointer to array containing ghost cells for fluxes in y-direction.*/
   real* fznbr; /**< Pointer to array containing ghost cells for fluxes in z-direction.*/
   
   bool allocateArray(const std::string& name,const uint& bytes,real*& arrptr);
   bool allocateArray(const std::string& name,const uint& bytes,uint*& arrptr);
};

#endif
