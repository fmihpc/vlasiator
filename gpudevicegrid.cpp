#include <cstdlib>
#include <iostream>
#include <limits>
#include <cuda_runtime.h>

#include "devicegrid.h"
#include "logger.h"

using namespace std;

extern Logger logger;

DeviceGrid::DeviceGrid() {
   bool success = true;
   cudaError_t error;
   
   logger << "(DEVICEGRID): Initializing." << endl;
   if (allocateArray("nbrsVel", MAX_VEL_BLOCKS*SIZE_NBRS_VEL*sizeof(uint), nbrsVel) == false) success = false;
   if (allocateArray("blockArray", MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real), blockArray) == false) success = false;
   if (allocateArray("blockParams", MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(real), blockParams) == false) success = false;
   if (allocateArray("cellParams", MAX_SPA_CELLS*SIZE_CELLPARAMS*sizeof(real), cellParams) == false) success = false;
   if (allocateArray("avgnbrx", MAX_VEL_BLOCKS*SIZE_BOUND*sizeof(real), avgnbrx) == false) success = false;
   if (allocateArray("avgnbry", MAX_VEL_BLOCKS*SIZE_BOUND*sizeof(real), avgnbry) == false) success = false;
   if (allocateArray("avgnbrz", MAX_VEL_BLOCKS*SIZE_BOUND*sizeof(real), avgnbrz) == false) success = false;
   if (allocateArray("fx     ", MAX_VEL_BLOCKS*SIZE_FLUXS*sizeof(real), fx) == false) success = false;
   if (allocateArray("fy     ", MAX_VEL_BLOCKS*SIZE_FLUXS*sizeof(real), fy) == false) success = false;
   if (allocateArray("fz     ", MAX_VEL_BLOCKS*SIZE_FLUXS*sizeof(real), fz) == false) success = false;
   if (allocateArray("fxnbr  ", MAX_VEL_BLOCKS*SIZE_BFLUX*sizeof(real), fxnbr) == false) success = false;
   if (allocateArray("fynbr  ", MAX_VEL_BLOCKS*SIZE_BFLUX*sizeof(real), fynbr) == false) success = false;
   if (allocateArray("fznbr  ", MAX_VEL_BLOCKS*SIZE_BFLUX*sizeof(real), fznbr) == false) success = false;
   if (allocateArray("d1x    ", MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real), d1x) == false) success = false;
   if (allocateArray("d1y    ", MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real), d1y) == false) success = false;
   if (allocateArray("d1z    ", MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real), d1z) == false) success = false;
   if (allocateArray("d2x    ", MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real), d2x) == false) success = false;
   if (allocateArray("d2y    ", MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real), d2y) == false) success = false;
   if (allocateArray("d2z    ", MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(real), d2z) == false) success = false;
   if (allocateArray("d1xnbr ", MAX_VEL_BLOCKS*SIZE_BDERI*sizeof(real), d1xnbr) == false) success = false;
   if (allocateArray("d1ynbr ", MAX_VEL_BLOCKS*SIZE_BDERI*sizeof(real), d1ynbr) == false) success = false;
   if (allocateArray("d1znbr ", MAX_VEL_BLOCKS*SIZE_BDERI*sizeof(real), d1znbr) == false) success = false;
   if (allocateArray("d2xnbr ", MAX_VEL_BLOCKS*SIZE_BDERI*sizeof(real), d2xnbr) == false) success = false;
   if (allocateArray("d2ynbr ", MAX_VEL_BLOCKS*SIZE_BDERI*sizeof(real), d2ynbr) == false) success = false;
   if (allocateArray("d2znbr ", MAX_VEL_BLOCKS*SIZE_BDERI*sizeof(real), d2znbr) == false) success = false;
   
   if (success == false) {
      nextInner = 0;
      nextBoundary = 0;
   } else {
      nextInner = 0;
      nextBoundary = MAX_VEL_BLOCKS-1;
   }
}

DeviceGrid::~DeviceGrid() {
   logger << "(DEVICEGRID): Shutting down." << endl;
   
   logger << "\t Deallocating 'nbrsVel': " << cudaGetErrorString( cudaFree(nbrsVel) ) << endl;
   logger << "\t Deallocating 'cellParams': " << cudaGetErrorString( cudaFree(cellParams) ) << endl;
   logger << "\t Deallocating 'blockParams': " << cudaGetErrorString( cudaFree(blockParams) ) << endl;
   logger << "\t Deallocating 'blockArray': " << cudaGetErrorString( cudaFree(blockArray) ) << endl;
   logger << "\t Deallocating 'avgnbrx': " << cudaGetErrorString( cudaFree(avgnbrx) ) << endl;
   logger << "\t Deallocating 'avgnbry': " << cudaGetErrorString( cudaFree(avgnbry) ) << endl;
   logger << "\t Deallocating 'avgnbrz': " << cudaGetErrorString( cudaFree(avgnbrz) ) << endl;
   logger << "\t Deallocating 'fx': " << cudaGetErrorString( cudaFree(fx) ) << endl;
   logger << "\t Deallocating 'fy': " << cudaGetErrorString( cudaFree(fy) ) << endl;
   logger << "\t Deallocating 'fz': " << cudaGetErrorString( cudaFree(fz) ) << endl;
   logger << "\t Deallocating 'fxnbr': " << cudaGetErrorString( cudaFree(fxnbr) ) << endl;
   logger << "\t Deallocating 'fynbr': " << cudaGetErrorString( cudaFree(fynbr) ) << endl;
   logger << "\t Deallocating 'fznbr': " << cudaGetErrorString( cudaFree(fznbr) ) << endl;
   logger << "\t Deallocating 'd1x': " << cudaGetErrorString( cudaFree(d1x) ) << endl;
   logger << "\t Deallocating 'd1y': " << cudaGetErrorString( cudaFree(d1y) ) << endl;
   logger << "\t Deallocating 'd1z': " << cudaGetErrorString( cudaFree(d1z) ) << endl;
   logger << "\t Deallocating 'd2x': " << cudaGetErrorString( cudaFree(d2x) ) << endl;
   logger << "\t Deallocating 'd2y': " << cudaGetErrorString( cudaFree(d2y) ) << endl;
   logger << "\t Deallocating 'd2z': " << cudaGetErrorString( cudaFree(d2z) ) << endl;
   logger << "\t Deallocating 'd1xnbr': " << cudaGetErrorString( cudaFree(d1xnbr) ) << endl;
   logger << "\t Deallocating 'd1ynbr': " << cudaGetErrorString( cudaFree(d1ynbr) ) << endl;
   logger << "\t Deallocating 'd1znbr': " << cudaGetErrorString( cudaFree(d1znbr) ) << endl;
   logger << "\t Deallocating 'd2xnbr': " << cudaGetErrorString( cudaFree(d2xnbr) ) << endl;
   logger << "\t Deallocating 'd2ynbr': " << cudaGetErrorString( cudaFree(d2ynbr) ) << endl;
   logger << "\t Deallocating 'd2znbr': " << cudaGetErrorString( cudaFree(d2znbr) ) << endl;
}

bool DeviceGrid::allocateArray(const std::string& name,const uint& bytes,real*& arrptr) {
   logger << "\t Reserving " << bytes << " bytes of GPU memory for array '" << name << "': ";
   cudaError_t error = cudaMalloc(&arrptr,bytes);
   logger << cudaGetErrorString(error) << endl;
   if (error == cudaSuccess) return true;
   return false;
}

bool DeviceGrid::allocateArray(const std::string& name,const uint& bytes,uint*& arrptr) {
   logger << "\t Reserving " << bytes << " bytes of GPU memory for array '" << name << "': ";
   cudaError_t error = cudaMalloc(&arrptr,bytes);
   logger << cudaGetErrorString(error) << endl;
   if (error == cudaSuccess) return true;
   return false;
}

/** Initialize a spatial cell so that its GPU array pointers point to correct 
 * locations in GPU arrays.
 * @param cell The cell whose pointers should be initialized.
 */
void DeviceGrid::initSpatialCell(SpatialCell& cell) {
   cell.gpu_cellParams  = cellParams  + cell.cellIndex*SIZE_CELLPARAMS;
   cell.gpu_avgs        = blockArray  + cell.velBlockIndex*SIZE_VELBLOCK;
   cell.gpu_nbrsVel     = nbrsVel     + cell.velBlockIndex*SIZE_NBRS_VEL;
   cell.gpu_blockParams = blockParams + cell.velBlockIndex*SIZE_BLOCKPARAMS;

   //cerr << "SpatialCell #" << cell.cellIndex << " gpu_avgs = " << cell.gpu_avgs << "\t blockArray = " << blockArray << endl;   
   //cerr << "GPU memory offsets:" << endl;
   //cerr << "\t gpu_avgs = " << cell.gpu_avgs - blockArray << endl;
   //cerr << "\t gpu_nbrsVel = " << cell.gpu_nbrsVel - nbrsVel << endl;
   //cerr << "\t gpu_blockParams = " << cell.gpu_blockParams - blockParams << endl;
}
