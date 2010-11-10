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
   
   //DeviceGrid(const std::vector<uint>& RealElemSizes,const std::vector<uint>& uintElemSizes);
   
   /** Deallocates the memory that has been previously allocated by a call to DeviceGrid().*/
   ~DeviceGrid();
   
   /** Query if the DeviceGrid is properly initialized. 
    * @return If true, all memory was allocated successfully and DeviceGrid is ready for use.
    * If false, an error has occurred and the simulation should be aborted.
    */
   bool isInitialized() const {return initialized;}
   
   void initSpatialCell(SpatialCell& cell);

   Real* getAvgnbrx() const {return avgnbrx;}
   Real* getAvgnbry() const {return avgnbry;}
   Real* getAvgnbrz() const {return avgnbrz;}
   Real* getBlockArray() const {return blockArray;}
   Real* getBlockParams() const {return blockParams;}
   Real* getCellParams() const {return cellParams;}
   Real* getD1x() const {return d1x;}
   Real* getD1y() const {return d1y;}
   Real* getD1z() const {return d1z;}
   Real* getD2x() const {return d2x;}
   Real* getD2y() const {return d2y;}
   Real* getD2z() const {return d2z;}
   Real* getD1xnbr() const {return d1xnbr;}
   Real* getD1ynbr() const {return d1ynbr;}
   Real* getD1znbr() const {return d1znbr;}
   Real* getD2xnbr() const {return d2xnbr;}
   Real* getD2ynbr() const {return d2ynbr;}
   Real* getD2znbr() const {return d2znbr;}
   Real* getFx() const {return fx;}
   Real* getFy() const {return fy;}
   Real* getFz() const {return fz;}
   Real* getFxnbr() const {return fxnbr;}
   Real* getFynbr() const {return fynbr;}
   Real* getFznbr() const {return fznbr;}
   uint* getNbrsSpa() const {return nbrsVel;}
   uint* getNbrsVel() const {return nbrsVel;}

   uint sizeBlockArray() const {return MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real);}
   uint sizeBlockParams() const {return MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(Real);}
   uint sizeD1x() const {return MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real);}
   uint sizeD1y() const {return MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real);}
   uint sizeD1z() const {return MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real);}
   uint sizeCellParams() const {return MAX_SPA_CELLS*SIZE_CELLPARAMS*sizeof(Real);}
   uint sizeNbrsVel() const {return MAX_VEL_BLOCKS*SIZE_NBRS_VEL*sizeof(uint);}
   
 private:
   bool initialized; /**< If true, DeviceGrid was initialized successfully.*/
   uint nextInner; /**< A pointer to a position in blockArray where the next inner velocity block can be written. */
   uint nextBoundary; /**< A pointer to a position in blockArray where the next boundary velocity block can be written. */
   
   Real* cellParams; /**< A pointer to array containing spatial cell parameters in device memory.*/
   
   Real* blockArray; /**< A pointer to the array containing velocity blocks in device memory. */
   Real* blockParams; /**< A pointer to an array containing velocity block parameter data in device memory.*/
   uint* nbrsVel; /**< A pointer to an array containing neighbour lists in velocity space.*/

   Real* d1x; /**< Pointer to array containing 1st derivatives to x-direction for each velocity grid block.*/
   Real* d1y; /**< Pointer to array containing 1st derivatives to y-direction for each velocity grid block.*/
   Real* d1z; /**< Pointer to array containing 1st derivatives to z-direction for each velocity grid block.*/
   Real* d2x; /**< Pointer to array containing 2nd derivatives to x-direction for each velocity grid block.*/
   Real* d2y; /**< Pointer to array containing 2nd derivatives to y-direction for each velocity grid block.*/
   Real* d2z; /**< Pointer to array containing 2nd derivatives to z-direction for each velocity grid block.*/
   Real* fx; /**< Pointer to array containing the fluxes to x-direction for each velocity grid block.*/
   Real* fy; /**< Pointer to array containing the fluxes to y-direction for each velocity grid block.*/
   Real* fz; /**< Pointer to array containing the fluxes to z-direction for each velocity grid block.*/
   Real* avgnbrx; /**< Pointer to array containing the ghost cells for volume averages in x-direction.*/
   Real* avgnbry; /**< Pointer to array containing the ghost cells for volume averages in y-direction.*/
   Real* avgnbrz; /**< Pointer to array containing the ghost cells for volume averages in z-direction.*/
   Real* d1xnbr; /**< Pointer to array containing ghost cells for 1st derivatives in x-direction. */
   Real* d1ynbr; /**< Pointer to array containing ghost cells for 1st derivatives in y-direction. */
   Real* d1znbr; /**< Pointer to array containing ghost cells for 1st derivatives in z-direction. */
   Real* d2xnbr; /**< Pointer to array containing ghost cells for 2nd derivatives in x-direction. */
   Real* d2ynbr; /**< Pointer to array containing ghost cells for 2nd derivatives in y-direction. */
   Real* d2znbr; /**< Pointer to array containing ghost cells for 2nd derivatives in z-direction. */
   Real* fxnbr; /**< Pointer to array containing ghost cells for fluxes in x-direction.*/
   Real* fynbr; /**< Pointer to array containing ghost cells for fluxes in y-direction.*/
   Real* fznbr; /**< Pointer to array containing ghost cells for fluxes in z-direction.*/
   
   bool allocateArray(const std::string& name,const uint& bytes,Real*& arrptr);
   bool allocateArray(const std::string& name,const uint& bytes,uint*& arrptr);
};

#endif
