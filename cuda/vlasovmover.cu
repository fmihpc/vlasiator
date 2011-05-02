#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <cuda_runtime.h>
#include <list>
#include <set>

#include "vlasovmover.h"
#include "devicegrid.h"
#include "cuda_acc_leveque.cu"
#include "cuda_trans_leveque.cu"
#include "streampool.h"
#include "eventcontainer.h"
#include "priorityqueue.h"
#include <blockarray.h>
#include <parameters.h>
#include <memalloc.h>
#include <mpilogger.h>

extern MPILogger mpilogger;

using namespace std;

namespace Main {
   #ifdef PARGRID
      vector<ID::type> cells;
   #else
      #include <stdint.h>
      vector<uint64_t> cells;
   #endif
   
   std::vector<const SpatialCell*> nbrPtrs(6,NULL);
   SpatialCell* cellPtr;
}

struct CellState {
   bool ghostFlag;
   bool deadFlag;
   bool remoteFlag;
   uint copyNbr;
   
   CellState(): ghostFlag(false),deadFlag(false),remoteFlag(false) { }
};

map<uint,CellState> cellStates;

const int MAX_COMPUTE_SUBMITS = 5;

DeviceGrid gpu_avgs;
DeviceGrid gpu_dfdt;
DeviceGrid gpu_flux;
DeviceGrid gpu_nbrsSpa;
DeviceGrid gpu_nbrsVel;
DeviceGrid gpu_blockParams;
DeviceGrid gpu_cellParams;

// Following are for reducing velocity moments:
DeviceGrid gpuRhoBuffers;
DeviceGrid gpuRhovxBuffers;
DeviceGrid gpuRhovyBuffers;
DeviceGrid gpuRhovzBuffers;
cuint VEL_MOMENT_BUFFERS = 100;

DeviceGrid gpuFluxBuffers; // Buffers for receiving flux updates from remote processes

int* gpu_dfdt_flags;

EventContainer computationEvents;
EventContainer transfersToCpuEvents;
EventContainer transfersToGpuEvents;
StreamPool streamPool;

PriorityQueue<uint> readyGpuCells;

struct RemoteUpdateInfo {
   uint computed;
   uint required;
   
   RemoteUpdateInfo(): computed(0),required(0) { }
};

map<uint,RemoteUpdateInfo> remoteUpdates;   

struct StencilAcceleration {
   list<uint> outgoing; /**< Remote neighbours who to send avgs.*/
   set<int> outgoingHosts;
};
map<uint,StencilAcceleration> stencilAcceleration;
set<uint> avgsReceives;

struct StencilSpatialFluxes {
   uint received;            /**< Number of received remote neighbour avgs so far.*/
   list<uint> incoming;      /**< Remote neighbours whom to receive avgs.*/
   list<uint> outgoing;      /**< Remote neighbours who to send dfdt.*/

   set<uint> receiveBuffers; /**< Receive buffers allocated to this cell (if any). */
   
   StencilSpatialFluxes(): received(0) { }
};
map<uint,Real*> spatialFluxSends;
map<uint,StencilSpatialFluxes> stencilSpatialFluxes;
multimap<uint,uint> remoteToLocalAvgs;

struct StencilSpatialProp {
   uint received;                    /**< Number of received neighbour dfdt received so far.*/         
   set<int> incoming;                /**< List of hosts whom to receive data.*/
   
   StencilSpatialProp(): received(0) { }
   ~StencilSpatialProp() { }
};
map<pair<uint,int>,uint> spatialFluxReceives;     // pair<local ID,host>,bufferID
map<uint,StencilSpatialProp> stencilSpatialProp;
vector<Real*> fluxReceiveBuffers;


/** Calculate the block index corresponding to the given 3-dimensional indices.
 * This is a temporary solution.
 */
uint calcBlockIndex(cuint& bix,cuint& biy,cuint& biz,cuint& N_blocks) {
   if (bix >= Parameters::vxblocks_ini) return N_blocks;
   if (biy >= Parameters::vyblocks_ini) return N_blocks;
   if (biz >= Parameters::vzblocks_ini) return N_blocks;
   return biz*Parameters::vyblocks_ini*Parameters::vxblocks_ini + biy*Parameters::vxblocks_ini + bix;
}

uint calcCellIndex(const uint& cellIndex,const int& i,const int& j,const int& k,ParGrid<SpatialCell>& mpiGrid,uint& exists) {
   typedef Parameters P;
   int ii = i;
   int jj = j;
   int kk = k;
   
   // If neighbouring cell does not exist, the index of this cell is returned instead. 
   // In practice this behaviour should only occur in boundary cells, and it forces the 
   // solver to fetch (non-existing) neighbour avgs from the boundary cell.
   exists = 1;
   if (ii < 0) {ii = 0; exists = 0;}
   if (jj < 0) {jj = 0; exists = 0;}
   if (kk < 0) {kk = 0; exists = 0;}
   if (ii > P::xcells_ini-1) {ii = P::xcells_ini-1; exists = 0;}
   if (jj > P::ycells_ini-1) {jj = P::ycells_ini-1; exists = 0;}
   if (kk > P::zcells_ini-1) {kk = P::zcells_ini-1; exists = 0;}
   uint index = kk*P::ycells_ini*P::xcells_ini + jj*P::xcells_ini + ii;
   
   // The index calculated above can be valid, but it does not mean that the 
   // cell exists. Here we check from the grid that the cell with computer index 
   // exists on this process:
   if (mpiGrid[index] == NULL) {
      exists = 0;
      return cellIndex;
   }
   return index;
}

bool bindTexture(texture<uint,1,cudaReadModeElementType>& texRef,uint* arrptr,cluint& BYTES,size_t& offset) {
   bool success = true;
   #ifndef NDEBUGCUDA
      cout << "bindTexture(uint) called with arrptr " << arrptr << " BYTES " << BYTES << endl;
      if (BYTES/sizeof(uint) > 134217728) cerr << "\t 1D uint array too large!" << endl;
   #endif
   
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindUnsigned);
   if (cudaBindTexture(NULL,&texRef,arrptr,&channelDesc,BYTES) != cudaSuccess) {
      cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' while binding textures!" << endl;
      success = false;
   } else {
      #ifndef NDEBUGCUDA
         cout << "\t 1D textures bound, offset = " << offset << endl;
      #endif
   }
   return success;
}

bool bindTexture(texture<Real,1,cudaReadModeElementType>& texRef,Real* arrptr,cluint& BYTES,size_t& offset) {
   bool success = true;
   #ifndef NDEBUGCUDA
      cout << "bindTexture(Real) called with arrptr " << arrptr << " BYTES " << BYTES << endl;
      if (BYTES/sizeof(Real) > 134217728) cerr << "\t 1D Real array too large!" << endl;
   #endif
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
   if (cudaBindTexture(&offset,&texRef,arrptr,&channelDesc,BYTES) != cudaSuccess) {
      cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while binding textures!" << endl;
      success = false;
   } else {
      #ifndef NDEBUGCUDA
         cout << "\t 1D textures bound, offset = " << offset << endl;
      #endif
   }
   return success;
}

bool bindTexture2D(texture<Real,2,cudaReadModeElementType>& texRef,Real* arrptr,cluint& height,size_t& offset) {
   bool success = true;
   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32,0,0,0,cudaChannelFormatKindFloat);
   if (cudaBindTexture2D(&offset,&texRef,arrptr,&channelDesc,CUDA_WIDTH,height,CUDA_WIDTH*sizeof(Real)) != cudaSuccess) {
      cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while binding textures!" << endl;
      success = false;
   } else {
      #ifndef NDEBUGCUDA
         cout << "\t 2D textures bound, width = " << CUDA_WIDTH << " height = " << height << " pitch = " << CUDA_WIDTH*sizeof(Real) << endl;
      #endif
   }
   return success;
}

void appendExistingBlocks(cuint& N_blocks,Real* blockParams,cuint& offset,map<uint,uint>& searchTree) {
   cerr << "Appending existing blocks, offset = " << offset << endl;
   // Go through all velocity blocks on this spatial cell, calculate their index and 
   // insert existing block indices to search tree:
   uint blockIndex;
   for (uint block=0; block<N_blocks; ++block) {
      BlockArray::calcIndex(blockParams + block*SIZE_BLOCKPARAMS,blockIndex);
      BlockArray::appendToSearchTree(blockIndex,offset + block,searchTree);
   }
}

// Create a spatial neighbour list entry for each velocity block in the given spatial cell.
// The last entry in the neighbour list of each block is a status word, in which each existing 
// neighbouring block has its bit flipped to unit value. Function spatialFluxes needs this information 
// when it is copying the calculated df/dt values to neighbours - changes calculated for 
// non-existing neighbours are simply thrown away, i.e. memory copies are not made.
void populateSpatNbrList(cuint& cellIndex,uint* tempNbrs,cuint& N_blocks,const DeviceGrid& gpu_avgs,ParGrid<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   
   uint tmpIndex;
   const int k = cellIndex / (P::ycells_ini*P::xcells_ini);
   tmpIndex = cellIndex - k*P::ycells_ini*P::xcells_ini;
   const int j = tmpIndex / P::xcells_ini;
   tmpIndex -= j*P::xcells_ini;
   const int i = tmpIndex;
   
   // Create a spatial neighbour list entry for each velocity block:
   for (uint bk=0; bk<P::vzblocks_ini; ++bk) for (uint bj=0; bj<P::vyblocks_ini; ++bj) for (uint bi=0; bi<P::vxblocks_ini; ++bi) {
      cuint offset      = calcBlockIndex(bi,bj,bk,N_blocks)*SIZE_NBRS_SPA;
      cuint blockOffset = calcBlockIndex(bi,bj,bk,N_blocks);
      int index=0;
      uint exists;
      tempNbrs[offset+30] = 0; // Status word - all existing neighbours have their bits flipped to unit value
      for (int kk=-1; kk<2; ++kk) for (int jj=-1; jj<2; ++jj) for (int ii=-1; ii<2; ++ii) {	 
	 cuint cellOffset = gpu_avgs.getOffset(calcCellIndex(cellIndex,i+ii,j+jj,k  ,mpiGrid,exists))/WID3;
	 if (cellOffset == numeric_limits<uint>::max()) {
	    cerr << "ERROR: Could not find cell #" << calcCellIndex(cellIndex,i+ii,j+jj,k  ,mpiGrid,exists) << " offset!" << endl;
	 }
	 tempNbrs[offset+index] = blockOffset + cellOffset;
	 tempNbrs[offset+30] = (tempNbrs[offset+30] | (exists << index));
	 ++index;
      }
      tempNbrs[offset+27] = blockOffset + gpu_avgs.getOffset(calcCellIndex(cellIndex,i-2,j  ,k  ,mpiGrid,exists))/WID3;  // --x neighbour
      tempNbrs[offset+30] = (tempNbrs[offset+30] | (exists << 27));
      tempNbrs[offset+28] = blockOffset + gpu_avgs.getOffset(calcCellIndex(cellIndex,i  ,j-2,k  ,mpiGrid,exists))/WID3;  // --y neighbour
      tempNbrs[offset+30] = (tempNbrs[offset+30] | (exists << 28));
      tempNbrs[offset+29] = blockOffset + gpu_avgs.getOffset(cellIndex)/WID3;                                     // --z, Periodic
      tempNbrs[offset+30] = (tempNbrs[offset+30] | (1 << 29));
      tempNbrs[offset+30] = (tempNbrs[offset+30] | (1 << 13));
   }
}

void populateVelNbrList(cuint& N_blocks,uint* tempNbrs,Real* blockParams) {
   uint ghostBlocksX = 0;
   uint ghostBlocksY = 0;
   uint ghostBlocksZ = 0;
   for (uint bk=0; bk<Parameters::vzblocks_ini; ++bk) for (uint bj=0; bj<Parameters::vyblocks_ini; ++bj) for (uint bi=0; bi<Parameters::vxblocks_ini; ++bi) {
      cuint offset = calcBlockIndex(bi,bj,bk,N_blocks)*28;
      
      tempNbrs[offset+ 0] = calcBlockIndex(bi-1,bj-1,bk-1,N_blocks);
      tempNbrs[offset+ 1] = calcBlockIndex(bi  ,bj-1,bk-1,N_blocks);
      tempNbrs[offset+ 2] = calcBlockIndex(bi+1,bj-1,bk-1,N_blocks);
      tempNbrs[offset+ 3] = calcBlockIndex(bi-1,bj  ,bk-1,N_blocks);
      tempNbrs[offset+ 4] = calcBlockIndex(bi  ,bj  ,bk-1,N_blocks);
      tempNbrs[offset+ 5] = calcBlockIndex(bi+1,bj  ,bk-1,N_blocks);
      tempNbrs[offset+ 6] = calcBlockIndex(bi-1,bj+1,bk-1,N_blocks);
      tempNbrs[offset+ 7] = calcBlockIndex(bi  ,bj+1,bk-1,N_blocks);
      tempNbrs[offset+ 8] = calcBlockIndex(bi+1,bj+1,bk-1,N_blocks);
         
      tempNbrs[offset+ 9] = calcBlockIndex(bi-1,bj-1,bk  ,N_blocks);
      tempNbrs[offset+10] = calcBlockIndex(bi  ,bj-1,bk  ,N_blocks);
      tempNbrs[offset+11] = calcBlockIndex(bi+1,bj-1,bk  ,N_blocks);
      tempNbrs[offset+12] = calcBlockIndex(bi-1,bj  ,bk  ,N_blocks);
      tempNbrs[offset+13] = calcBlockIndex(bi  ,bj  ,bk  ,N_blocks);
      tempNbrs[offset+14] = calcBlockIndex(bi+1,bj  ,bk  ,N_blocks);
      tempNbrs[offset+15] = calcBlockIndex(bi-1,bj+1,bk  ,N_blocks);
      tempNbrs[offset+16] = calcBlockIndex(bi  ,bj+1,bk  ,N_blocks);
      tempNbrs[offset+17] = calcBlockIndex(bi+1,bj+1,bk  ,N_blocks);
      
      tempNbrs[offset+18] = calcBlockIndex(bi-1,bj-1,bk+1,N_blocks);
      tempNbrs[offset+19] = calcBlockIndex(bi  ,bj-1,bk+1,N_blocks);
      tempNbrs[offset+20] = calcBlockIndex(bi+1,bj-1,bk+1,N_blocks);
      tempNbrs[offset+21] = calcBlockIndex(bi-1,bj  ,bk+1,N_blocks);
      tempNbrs[offset+22] = calcBlockIndex(bi  ,bj  ,bk+1,N_blocks);
      tempNbrs[offset+23] = calcBlockIndex(bi+1,bj  ,bk+1,N_blocks);
      tempNbrs[offset+24] = calcBlockIndex(bi-1,bj+1,bk+1,N_blocks);
      tempNbrs[offset+25] = calcBlockIndex(bi  ,bj+1,bk+1,N_blocks);
      tempNbrs[offset+26] = calcBlockIndex(bi+1,bj+1,bk+1,N_blocks);
      
      // Index of +vx face neighbour is at tempNbrs[14]
      //          +vy                      tempNbrs[16]
      //          +vz                      tempNbrs[22]
      if (tempNbrs[offset+14] == N_blocks) ++ghostBlocksX;
      if (tempNbrs[offset+16] == N_blocks) ++ghostBlocksY;
      if (tempNbrs[offset+22] == N_blocks) ++ghostBlocksZ;
   }
}

bool finalizeMover() {
   bool success = true;
   gpu_avgs.finalize();
   gpu_dfdt.finalize();
   gpu_flux.finalize();
   gpu_nbrsSpa.finalize();
   gpu_nbrsVel.finalize();
   gpu_blockParams.finalize();
   gpu_cellParams.finalize();

   gpuRhoBuffers.finalize();
   gpuRhovxBuffers.finalize();
   gpuRhovyBuffers.finalize();
   gpuRhovzBuffers.finalize();
   
   gpuFluxBuffers.finalize();
   
   // Free page-locked host memory
   for (map<uint,Real*>::iterator it = spatialFluxSends.begin(); it != spatialFluxSends.end(); ++it) 
     freeArray(it->second);

   for (size_t i=0; i<fluxReceiveBuffers.size(); ++i) {
      freeArray(fluxReceiveBuffers[i]);
   }
   
   return success;
}

bool initCudaDevices(const int& mpirank) {
   bool success = true;
   int deviceNumber = -1000;
   int cudaDevices;
   cudaGetDeviceCount(&cudaDevices);
   
   mpilogger << "Number of CUDA devices is #" << cudaDevices << endl;
   mpilogger << "Trying to set device #" << mpirank << " for this process." << endl;
   if (cudaSetDevice(mpirank) != cudaSuccess) {
      mpilogger << "\t ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred when trying to change device!" << endl;
      cerr << "CRITICAL ERROR occurred when trying to set device!" << endl;
      success = false;
   }
   cudaGetDevice(&deviceNumber);
   mpilogger << "\t CUDA reports that device #" << deviceNumber << " is now in use." << endl << write;

   cerr << "Proc #" << mpirank << " is using device #" << deviceNumber << endl;
   
   if (deviceNumber != mpirank) success = false;   
   return success;
}

bool initializeMover(ParGrid<SpatialCell>& mpiGrid) {
   bool success = true;
   
   size_t arrayHeight = MAX_VEL_BLOCKS*WID3 / CUDA_WIDTH;
   if (MAX_VEL_BLOCKS*WID3 % CUDA_WIDTH > 0) ++arrayHeight;
   
   // Initialize device arrays. Arrays avgs, dfdt and flux may be so large that 1D textures 
   // cannot be used to fetch values, thus some extra padding is necessary for 2D textures:
   if (gpu_avgs.initialize(CUDA_WIDTH*arrayHeight,sizeof(Real)) == false) success = false;
   if (gpu_dfdt.initialize(CUDA_WIDTH*arrayHeight,sizeof(Real)) == false) success = false;
   if (gpu_flux.initialize(CUDA_WIDTH*arrayHeight,sizeof(Real)) == false) success = false;
   if (gpu_nbrsSpa.initialize(MAX_VEL_BLOCKS*SIZE_NBRS_SPA,sizeof(uint)) == false) success = false;
   if (gpu_nbrsVel.initialize(MAX_VEL_BLOCKS*28,sizeof(uint)) == false) success = false; // HACK
   if (gpu_blockParams.initialize(MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS,sizeof(Real)) == false) success = false;
   if (gpu_cellParams.initialize(MAX_SPAT_CELLS*SIZE_CELLPARAMS,sizeof(Real)) == false) success = false;

   if (gpuRhoBuffers.initialize(VEL_MOMENT_BUFFERS*1024,sizeof(Real)) == false) success = false;
   if (gpuRhovxBuffers.initialize(VEL_MOMENT_BUFFERS*1024,sizeof(Real)) == false) success = false;
   if (gpuRhovyBuffers.initialize(VEL_MOMENT_BUFFERS*1024,sizeof(Real)) == false) success = false;
   if (gpuRhovzBuffers.initialize(VEL_MOMENT_BUFFERS*1024,sizeof(Real)) == false) success = false;

   if (success == false) {
      cerr << "CRITICAL ERROR: Failed to allocate memory, exiting!" << endl;
      exit(1);
   }
   
   // Initialize velocity moment buffers to zero value:
   Real* tmpMoments = new Real[VEL_MOMENT_BUFFERS*1024];
   for (size_t i=0; i<VEL_MOMENT_BUFFERS*1024; ++i) tmpMoments[i] = 0.0;
   cudaMemcpy(gpuRhoBuffers.getArray<Real>()  ,tmpMoments,VEL_MOMENT_BUFFERS*1024*sizeof(Real),cudaMemcpyHostToDevice);
   cudaMemcpy(gpuRhovxBuffers.getArray<Real>(),tmpMoments,VEL_MOMENT_BUFFERS*1024*sizeof(Real),cudaMemcpyHostToDevice);
   cudaMemcpy(gpuRhovyBuffers.getArray<Real>(),tmpMoments,VEL_MOMENT_BUFFERS*1024*sizeof(Real),cudaMemcpyHostToDevice);
   cudaMemcpy(gpuRhovzBuffers.getArray<Real>(),tmpMoments,VEL_MOMENT_BUFFERS*1024*sizeof(Real),cudaMemcpyHostToDevice);
   delete tmpMoments;
   tmpMoments = NULL;
   
   // Bind textures to device arrays:
   size_t texOffset = 0;    
   bindTexture2D(texRef_avgs2D,gpu_avgs.getArray<Real>(),arrayHeight,texOffset);
   bindTexture2D(texRef_flux2D,gpu_flux.getArray<Real>(),arrayHeight,texOffset);
   bindTexture2D(texRef_dfdt2D,gpu_dfdt.getArray<Real>(),arrayHeight,texOffset);
   bindTexture(texRef_nbrsSpa    ,gpu_nbrsSpa.getArray<uint>()    ,MAX_VEL_BLOCKS*SIZE_NBRS_SPA   *sizeof(uint),texOffset);
   bindTexture(texRef_nbrsVel    ,gpu_nbrsVel.getArray<uint>()    ,MAX_VEL_BLOCKS*28              *sizeof(uint),texOffset); // HACK
   bindTexture(texRef_blockParams,gpu_blockParams.getArray<Real>(),MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(Real),texOffset);
   bindTexture(texRef_cellParams ,gpu_cellParams.getArray<Real>() ,MAX_SPAT_CELLS*SIZE_CELLPARAMS *sizeof(Real),texOffset);
   
   // Allocate a mutex for each velocity block on the device, and initialize their 
   // values to zero:
   if (cudaMalloc(&gpu_dfdt_flags,MAX_VEL_BLOCKS*sizeof(uint)) != cudaSuccess) 
     cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while allocating gpu_dfdt_flags!" << endl;
   uint* tmp = new uint[MAX_VEL_BLOCKS];
   for (uint i=0; i<MAX_VEL_BLOCKS; ++i) tmp[i] = 0;
   if (cudaMemcpy(gpu_dfdt_flags,tmp,MAX_VEL_BLOCKS*sizeof(uint),cudaMemcpyHostToDevice) != cudaSuccess) 
     cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying zeroes to gpu_dfdt_flags!" << endl;
   delete tmp;
   tmp = NULL;

   // Allocate arrays for velocity blocks on device:
   Real cpuBoundaryBlock[WID3];
   for (uint i=0; i<WID3; ++i) cpuBoundaryBlock[i] = 0.0;
   
   cuint N_ghostBlocks = 1000;
   cuint ghostAvgsOffset = gpu_avgs.reserveArray(numeric_limits<uint>::max()-1,N_ghostBlocks*WID3);
   cuint ghostFluxOffset = gpu_flux.reserveArray(numeric_limits<uint>::max()-1,N_ghostBlocks*WID3);

   mpiGrid.getAllCells(Main::cells);
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellIndex = Main::cells[cell];
      Main::cellPtr = mpiGrid[cellIndex];

      cuint offsetAvgs        = gpu_avgs.reserveArray(cellIndex,(Main::cellPtr->N_blocks+1)*WID3);
      cuint offsetDerivs      = gpu_dfdt.reserveArray(cellIndex,Main::cellPtr->N_blocks*WID3);
      cuint offsetFlux        = gpu_flux.reserveArray(cellIndex,(Main::cellPtr->N_blocks+1)*WID3);
      cuint offsetNbrs        = gpu_nbrsVel.reserveArray(cellIndex,Main::cellPtr->N_blocks*28); // HACK
      cuint offsetBlockParams = gpu_blockParams.reserveArray(cellIndex,Main::cellPtr->N_blocks*SIZE_BLOCKPARAMS);
      cuint offsetCellParams  = gpu_cellParams.reserveArray(cellIndex,SIZE_CELLPARAMS);

      if (offsetAvgs == numeric_limits<uint>::max()) {cerr << "offsetAvgs = " << offsetAvgs << " cell #" << cellIndex << endl; exit(1);}
      if (offsetNbrs == numeric_limits<uint>::max()) {cerr << "offsetNbrs = " << offsetNbrs << " cell #" << cellIndex << endl; exit(1);}
      
      Real* gpuAvgs        = gpu_avgs.getArray<Real>(offsetAvgs);
      Real* gpuFlux        = gpu_flux.getArray<Real>(offsetFlux);
      Real* gpuDerivs      = gpu_dfdt.getArray<Real>(offsetDerivs);
      uint* gpuNbrs        = gpu_nbrsVel.getArray<uint>(offsetNbrs);
      Real* gpuBlockParams = gpu_blockParams.getArray<Real>(offsetBlockParams);
      Real* gpuCellParams  = gpu_cellParams.getArray<Real>(offsetCellParams);
      
      uint* tempNbrs = new uint[28*Main::cellPtr->N_blocks];
      populateVelNbrList(Main::cellPtr->N_blocks,tempNbrs,Main::cellPtr->cpu_blockParams); // HACK
      
      cudaMemcpy(gpuAvgs+Main::cellPtr->N_blocks*WID3,cpuBoundaryBlock,WID3*sizeof(Real),cudaMemcpyHostToDevice); // TEST
      cudaMemcpy(gpuFlux+Main::cellPtr->N_blocks*WID3,cpuBoundaryBlock,WID3*sizeof(Real),cudaMemcpyHostToDevice);
      
      if (cudaMemcpy(gpuAvgs       ,Main::cellPtr->cpu_avgs       ,Main::cellPtr->N_blocks*WID3*sizeof(Real)   ,cudaMemcpyHostToDevice) != cudaSuccess) 
	cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying cpu_avgs of cell #" << cellIndex << endl;
      if (cudaMemcpy(gpuNbrs       ,tempNbrs                      ,Main::cellPtr->N_blocks*28*sizeof(uint)              ,cudaMemcpyHostToDevice) != cudaSuccess)
	cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying velNbrs of cell #" << cellIndex << endl;
      if (cudaMemcpy(gpuBlockParams,Main::cellPtr->cpu_blockParams,Main::cellPtr->N_blocks*SIZE_BLOCKPARAMS*sizeof(Real),cudaMemcpyHostToDevice) != cudaSuccess)
	cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying cpu_blockParams of cell #" << cellIndex << endl;
      if (cudaMemcpy(gpuCellParams ,Main::cellPtr->cpu_cellParams ,SIZE_CELLPARAMS*sizeof(Real)                         ,cudaMemcpyHostToDevice) != cudaSuccess)
	cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying cpu_cellParams of cell #" << cellIndex << endl;
      delete tempNbrs;
      tempNbrs = NULL;
   }
   
   // Populate spatial neighbour lists:
   mpiGrid.getCells(Main::cells);
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellIndex = Main::cells[cell];
      Main::cellPtr = mpiGrid[cellIndex];
      
      cuint offsetNbrs = gpu_nbrsSpa.reserveArray(cellIndex,(Main::cellPtr->N_blocks)*SIZE_NBRS_SPA);
      uint* gpuNbrsSpa    = gpu_nbrsSpa.getArray<uint>(offsetNbrs);
      
      uint* tempNbrs = new uint[SIZE_NBRS_SPA*Main::cellPtr->N_blocks];
      populateSpatNbrList(cellIndex,tempNbrs,Main::cellPtr->N_blocks,gpu_avgs,mpiGrid);
      if (cudaMemcpy(gpuNbrsSpa,tempNbrs,Main::cellPtr->N_blocks*SIZE_NBRS_SPA*sizeof(uint),cudaMemcpyHostToDevice) != cudaSuccess) 
	cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying spatNbrs of cell #" << cellIndex << endl;
      
      delete tempNbrs;
      tempNbrs = NULL;
   }

   // Ghost cell classification. Cell is a ghost cell if even one 
   // face neighbour does not exist.
   uint ghostCellCount = 0;
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellIndex = Main::cells[cell];
      const int j = cellIndex / Parameters::xcells_ini;
      const int i = cellIndex - j*Parameters::xcells_ini;
      
      int existingNeighbours = 0;
      uint exists;
      
      // +/- x-neighbours:
      calcCellIndex(cellIndex,i-1,j  ,0,mpiGrid,exists);
      if (exists > 0) ++existingNeighbours;
      calcCellIndex(cellIndex,i+1,j  ,0,mpiGrid,exists);
      if (exists > 0) ++existingNeighbours;
      
      // +/- y-neighbours:
      calcCellIndex(cellIndex,i  ,j-1,0,mpiGrid,exists);
      if (exists > 0) ++existingNeighbours;
      calcCellIndex(cellIndex,i  ,j+1,0,mpiGrid,exists);
      if (exists > 0) ++existingNeighbours;
      
      // z-periodic:
      
      cellStates[cellIndex];
      if (existingNeighbours != 4) {
	 cellStates[cellIndex].ghostFlag = true;
	 ++ghostCellCount;
      } else cellStates[cellIndex].ghostFlag = false;
   }
   // TEST: Insert remote neighbours to cellStates as ghosts:
   mpiGrid.getRemoteCells(Main::cells);
   for (uint i=0; i<Main::cells.size(); ++i) {
      cellStates[Main::cells[i]].ghostFlag = true;
      cellStates[Main::cells[i]].remoteFlag = true;
   }   
   cout << "Total number of ghost cells: " << ghostCellCount << endl;

   // Dead cell classification. Cell is dead if it has face neighbours, 
   // but they are all non-remote ghost cells.
   uint deadCellCount = 0;
   mpiGrid.getCells(Main::cells);
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellIndex = Main::cells[cell];
      if (cellStates[cellIndex].ghostFlag == false) continue;
      
      const int j = cellIndex / Parameters::xcells_ini;
      const int i = cellIndex - j*Parameters::xcells_ini;
      
      int existingNeighbours = 0;
      uint existsLeft,existsRight,existsBottom,existsTop;
      
      cuint nbrIndexLeft = calcCellIndex(cellIndex,i-1,j  ,0,mpiGrid,existsLeft);
      if (existsLeft > 0) {
	 if (cellStates.find(nbrIndexLeft) == cellStates.end()) 
	   cerr << "Could not find " << nbrIndexLeft << " from cellStates!" << endl;
	 else {
	    if (cellStates[nbrIndexLeft].remoteFlag == true) ++existingNeighbours; 
	    if (cellStates[nbrIndexLeft].ghostFlag == false) ++existingNeighbours;
	 }
      }
      cuint nbrIndexRight = calcCellIndex(cellIndex,i+1,j  ,0,mpiGrid,existsRight);
      if (existsRight > 0) {
	 if (cellStates.find(nbrIndexRight) == cellStates.end())
	   cerr << "Could not find " << nbrIndexRight << " from cellStates!" << endl;
	 else {
	    if (cellStates[nbrIndexRight].remoteFlag == true) ++existingNeighbours;
	    if (cellStates[nbrIndexRight].ghostFlag == false) ++existingNeighbours;
	 }
      }
      cuint nbrIndexBottom = calcCellIndex(cellIndex,i  ,j-1,0,mpiGrid,existsBottom);
      if (existsBottom > 0) {
	 if (cellStates.find(nbrIndexBottom) == cellStates.end())
	   cerr << "Could not find " << nbrIndexBottom << " from cellStates!" << endl;
	 else {
	    if (cellStates[nbrIndexBottom].remoteFlag == true) ++existingNeighbours;
	    if (cellStates[nbrIndexBottom].ghostFlag != true) ++existingNeighbours;
	 }
      }
      cuint nbrIndexTop = calcCellIndex(cellIndex,i  ,j+1,0,mpiGrid,existsTop);
      if (existsTop > 0) {
	 if (cellStates.find(nbrIndexTop) == cellStates.end())
	   cerr << "Could not find " << nbrIndexTop << " from cellStates!" << endl;
	 else {
	    if (cellStates[nbrIndexTop].remoteFlag == true) ++existingNeighbours;
	    if (cellStates[nbrIndexTop].ghostFlag != true) ++existingNeighbours;
	 }
      }
      // Mark cell as dead if all its existing neighbours are ghost cells:
      if (existingNeighbours == 0) {
	 cellStates[cellIndex].deadFlag = true;
	 ++deadCellCount;
      } else cellStates[cellIndex].deadFlag = false;
      
      // Determine how the boundary condition for this cell is calculated. For cells 
      // on the external boundary do nothing (leave initial state), for internal boundaries 
      // copy from an existing face neighbour:
      bool noCopy = false;
      if (i == Parameters::xcells_ini-1) noCopy = true;
      if (j == 0 || j == Parameters::ycells_ini-1) noCopy = true;
      if (noCopy == true) {
	 cellStates[cellIndex].copyNbr = numeric_limits<uint>::max();
	 continue;
      }
      
      // Cell is an internal boundary cell. Figure out the neighbour whom to copy volume averages:
      if (existsRight > 0 && cellStates[nbrIndexRight].ghostFlag != true) cellStates[cellIndex].copyNbr = nbrIndexRight;
      else if (existsBottom > 0 && cellStates[nbrIndexBottom].ghostFlag != true) cellStates[cellIndex].copyNbr = nbrIndexBottom;
      else if (existsTop > 0 && cellStates[nbrIndexTop].ghostFlag != true) cellStates[cellIndex].copyNbr = nbrIndexTop;
      else cellStates[cellIndex].copyNbr = numeric_limits<uint>::max();
      
      cuint j_nbr = cellStates[cellIndex].copyNbr/Parameters::xcells_ini;
      cuint i_nbr = cellStates[cellIndex].copyNbr - j_nbr*Parameters::xcells_ini;
   }
   cout << "Total number of dead cells: " << deadCellCount << endl;

   // Init streams:
   streamPool.initialize(100,15);

   avgsReceives.clear();
   stencilAcceleration.clear();
   stencilSpatialFluxes.clear();
   stencilSpatialProp.clear();
   mpiGrid.getCells(Main::cells);
   for (size_t i=0; i<Main::cells.size(); ++i) {
      // Fill stencils:
      stencilAcceleration[Main::cells[i]];
      stencilSpatialFluxes[Main::cells[i]];
      stencilSpatialProp[Main::cells[i]];
            
      int nbrHost;
      for (uchar nbrTypeID=9; nbrTypeID<18; ++nbrTypeID) {
	 if (nbrTypeID == 13) continue; // This cell
	 cuint nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],nbrTypeID);
	 if (nbrCellID == numeric_limits<uint>::max()) continue; // Neighbour does not exist or this process has it
	 //cerr << "Proc #" << mpiGrid.rank() << " local #" << Main::cells[i] << " nbrType #" << (int)nbrTypeID << " index = " << nbrCellID << endl;

	 stencilSpatialFluxes[Main::cells[i]].outgoing.push_back(nbrCellID);
	 spatialFluxSends.insert(make_pair(nbrCellID,(Real*)NULL));
	 
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilSpatialProp[Main::cells[i]].incoming.insert(nbrHost);
	 spatialFluxReceives.insert(make_pair(make_pair(Main::cells[i],nbrHost),0));
	 
	 ++remoteUpdates[nbrCellID].required;
      }
      
      // Face neighbours:
      uint nbrCellID;
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],12);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 avgsReceives.insert(nbrCellID);
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoingHosts.insert(nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoing.push_back(nbrCellID);
	 stencilSpatialFluxes[Main::cells[i]].incoming.push_back(nbrCellID);
	 remoteToLocalAvgs.insert(make_pair(nbrCellID,Main::cells[i]));
      }
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],14);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 avgsReceives.insert(nbrCellID);
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoingHosts.insert(nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoing.push_back(nbrCellID);
	 stencilSpatialFluxes[Main::cells[i]].incoming.push_back(nbrCellID);
	 remoteToLocalAvgs.insert(make_pair(nbrCellID,Main::cells[i]));
      }
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],10);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 avgsReceives.insert(nbrCellID);
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoingHosts.insert(nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoing.push_back(nbrCellID);
	 stencilSpatialFluxes[Main::cells[i]].incoming.push_back(nbrCellID);
	 remoteToLocalAvgs.insert(make_pair(nbrCellID,Main::cells[i]));
      }
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],16);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 avgsReceives.insert(nbrCellID);
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoingHosts.insert(nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoing.push_back(nbrCellID);
	 stencilSpatialFluxes[Main::cells[i]].incoming.push_back(nbrCellID);
	 remoteToLocalAvgs.insert(make_pair(nbrCellID,Main::cells[i]));
      }
      
      // x+2,y+2 neighbours to outgoing avgs:
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],28);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoingHosts.insert(nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoing.push_back(nbrCellID);
      }
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],30);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 mpiGrid.getHost(nbrCellID,nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoingHosts.insert(nbrHost);
	 stencilAcceleration[Main::cells[i]].outgoing.push_back(nbrCellID);
      }
      
      // x-2,y-2 neighbours to incoming avgs:
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],27);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 avgsReceives.insert(nbrCellID);
	 stencilSpatialFluxes[Main::cells[i]].incoming.push_back(nbrCellID);
	 remoteToLocalAvgs.insert(make_pair(nbrCellID,Main::cells[i]));
      }
      nbrCellID = mpiGrid.getRemoteNeighbour(Main::cells[i],29);
      if (nbrCellID != numeric_limits<uint>::max()) {
	 avgsReceives.insert(nbrCellID);
	 stencilSpatialFluxes[Main::cells[i]].incoming.push_back(nbrCellID);
	 remoteToLocalAvgs.insert(make_pair(nbrCellID,Main::cells[i]));
      }
      
   }

   // Allocate memory for flux receive buffers:
   const size_t fluxReceiveBufferSize = mpiGrid[Main::cells[0]]->N_blocks*WID3;
   fluxReceiveBuffers.resize(spatialFluxReceives.size());
   uint bufferCounter = 0;
   for (map<pair<uint,int>,uint>::iterator it=spatialFluxReceives.begin(); it!=spatialFluxReceives.end(); ++it) {
      it->second = bufferCounter;
      stencilSpatialFluxes[it->first.first].receiveBuffers.insert(bufferCounter);      
      ++bufferCounter;
   }
   for (size_t i=0; i<fluxReceiveBuffers.size(); ++i) {
      allocateArray(&(fluxReceiveBuffers[i]),fluxReceiveBufferSize);
   }
   
   // Each cell that receives a remote flux update needs a buffer on the GPU as well. Count the number of required buffers:
   uint N_requiredBuffers = 0;
   for (map<uint,StencilSpatialFluxes>::const_iterator it=stencilSpatialFluxes.begin(); it!=stencilSpatialFluxes.end(); ++it) {
      if (it->second.receiveBuffers.size() > 0) ++N_requiredBuffers;
   }
   mpilogger << "Need GPU flux transfer buffers for " << N_requiredBuffers << " cells" << endl << write;
   if (gpuFluxBuffers.initialize(N_requiredBuffers*fluxReceiveBufferSize,sizeof(Real)) == false) {
      mpilogger << "CRITICAL ERROR: Could not allocate enough device memory for GPU flux receive buffers!" << endl << write;
      exit(1);
   }
   
   for (map<uint,StencilSpatialFluxes>::const_iterator it=stencilSpatialFluxes.begin(); it!=stencilSpatialFluxes.end(); ++it) {
      if (it->second.receiveBuffers.size() == 0) continue;
      if (gpuFluxBuffers.reserveArray(it->first,fluxReceiveBufferSize) == numeric_limits<uint>::max()) {
	 mpilogger << "CRITICAL ERROR: Could not allocate flux transfer buffer on GPU for cell #" << it->first << endl << write;
	 exit(1);
      }
   }

   // Allocate memory for flux send buffers:
   const size_t fluxSendBufferSize = mpiGrid[Main::cells[0]]->N_blocks*WID3;
   for (map<uint,Real*>::iterator it=spatialFluxSends.begin(); it!=spatialFluxSends.end(); ++it) {
      allocateArray(&(it->second),fluxSendBufferSize);
   }
   return success;
}

void calculateAcceleration(ParGrid<SpatialCell>& mpiGrid) {   
   dim3 gridSize;
   dim3 blockSize(WID,WID,WID);

   cuint dynamicMemory = 0;
   uint cellsToCompute = 0;
   uint computedCells = 0;
   uint cellIndex;
   uint streamID;
   uint priority;
   uint64_t cellID;
   cudaStream_t stream;
   int counter = 0;
   // Initialization. Push all local cells to readyGpuCells. 
   // Number of cell's remote neighbours is used as its priority.
   mpiGrid.getCells(Main::cells);
   vector<uint>::iterator cellIterator = Main::cells.begin();   
   for (cellIterator = Main::cells.begin(); cellIterator != Main::cells.end(); ++cellIterator) {
      if (cellStates[*cellIterator].ghostFlag == true && stencilAcceleration[*cellIterator].outgoing.size() == 0) continue;
      readyGpuCells.insert(*cellIterator,mpiGrid.getNumberOfRemoteNeighbours(*cellIterator));
      ++cellsToCompute;
   }

   cuint avgsByteSize = mpiGrid[Main::cells[0]]->N_blocks*WID3*sizeof(Real);
   
   // Post receives for avgs (needed in spatial fluxes):
   mpiGrid.startSingleMode();
   for (map<uint,StencilSpatialFluxes>::iterator it=stencilSpatialFluxes.begin(); it!=stencilSpatialFluxes.end(); ++it) it->second.received = 0;
   for (set<uint>::const_iterator it=avgsReceives.begin(); it!=avgsReceives.end(); ++it) {
      mpiGrid.singleReceive(*it,*it,avgsByteSize,reinterpret_cast<char*>(mpiGrid[*it]->cpu_avgs),*it);
   }
   
   bool allTasksCompleted = false;
   while (allTasksCompleted == false) {
      allTasksCompleted = true;
      
      // Check for completed GPU -> CPU transfers:
      while (transfersToCpuEvents.getCompleted(streamID,cellID) == true) {
	 if (streamPool.insert(streamID) == false)
	   cerr << "ERROR occurred while inserting idle streams to streampool!" << endl;
	 
	 // A cell with remote neighbour(s) has been transferred to CPU. Send it 
	 // to remote neighbours:
	 for (set<int>::const_iterator destHost = stencilAcceleration[cellID].outgoingHosts.begin(); destHost != stencilAcceleration[cellID].outgoingHosts.end(); ++destHost) {
	    mpiGrid.singleSend(*destHost,cellID,avgsByteSize,reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs),cellID);
	 }
      }
      if (transfersToCpuEvents.empty() == false) allTasksCompleted = false;
      
      // Check for completed computations:
      while (computationEvents.getCompleted(streamID,cellID) == true) {
	 // Test if the cell has remote neighbours and needs to be 
	 // transferred to CPU and sent to other process(es):
	 if (stencilAcceleration[cellID].outgoingHosts.size() > 0) {
	    Main::cellPtr = mpiGrid[cellID];
	    cuint offsetAvgs = gpu_avgs.getOffset(cellID);
	    Real* gpuAvgs    = gpu_avgs.getArray<Real>(offsetAvgs);
	    cudaMemcpyAsync(Main::cellPtr->cpu_avgs,gpuAvgs,Main::cellPtr->N_blocks*WID3*sizeof(Real),cudaMemcpyDeviceToHost,streamPool.get(streamID));
	    transfersToCpuEvents.insert(streamID,streamPool.get(streamID),cellID);
	    allTasksCompleted = false;
	 } else {
	    streamPool.insert(streamID);
	 }
	 ++computedCells;
      }
      if (computedCells != cellsToCompute) allTasksCompleted = false;
      
      // Check if streams are available for computing cells:
      counter = 0;
      while (readyGpuCells.empty() == false) {
	 if (streamPool.get(streamID,stream,false) == true) {
	    readyGpuCells.pop(cellIndex,priority);
	    Main::cellPtr = mpiGrid[cellIndex];
	    
	    // Set stream to calculate the next uncalculated cell:
	    cuint offsetAvgs        = gpu_avgs.getOffset(cellIndex);
	    cuint offsetDerivs      = gpu_dfdt.getOffset(cellIndex);
	    cuint offsetFlux        = gpu_flux.getOffset(cellIndex);
	    cuint offsetNbrs        = gpu_nbrsVel.getOffset(cellIndex);
	    cuint offsetBlockParams = gpu_blockParams.getOffset(cellIndex);
	    cuint offsetCellParams  = gpu_cellParams.getOffset(cellIndex);
	    
	    Real* gpuAvgs        = gpu_avgs.getArray<Real>(offsetAvgs);
	    Real* gpuFlux        = gpu_flux.getArray<Real>(offsetFlux);
	    Real* gpuDerivs      = gpu_dfdt.getArray<Real>(offsetDerivs);
	    uint* gpuNbrs        = gpu_nbrsVel.getArray<uint>(offsetNbrs);
	    Real* gpuBlockParams = gpu_blockParams.getArray<Real>(offsetBlockParams);
	    Real* gpuCellParams  = gpu_cellParams.getArray<Real>(offsetCellParams);
	    
	    gridSize.x = Main::cellPtr->N_blocks;
	    gridSize.y = 1;

	    // Calculate contribution to df/dt coming from vx-faces:
	    cuda_acc_xfluxes<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_xfluxes3DcorrY<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_xfluxes3DcorrZ<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_vx_changes<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuDerivs,offsetDerivs,offsetFlux,offsetNbrs,offsetBlockParams,Parameters::dt);
	    // Calculate contribution to df/dt coming from vy-faces:
	    cuda_acc_yfluxes<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_yfluxes3DcorrX<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_yfluxes3DcorrZ<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_vy_changes<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuDerivs,offsetDerivs,offsetFlux,offsetNbrs,offsetBlockParams,Parameters::dt);
	    // Calculate contribution to df/dt coming from vz-faces:
	    cuda_acc_zfluxes<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_zfluxes3DcorrX<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_zfluxes3DcorrY<<<gridSize,blockSize,dynamicMemory,stream>>>(offsetAvgs,offsetFlux,offsetBlockParams,offsetCellParams,offsetNbrs,gpuFlux,Parameters::dt,Parameters::q_per_m);
	    cuda_acc_vz_changes<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuDerivs,offsetDerivs,offsetFlux,offsetNbrs,offsetBlockParams,Parameters::dt);
	    // Propagate volume averages in velocity space:
	    cuda_acc_propagate<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuAvgs,gpuDerivs,offsetAvgs,offsetDerivs);

	    // Insert computation event to computationEvents:
	    computationEvents.insert(streamID,stream,cellIndex);
	    ++cellIterator;
	 } else break;
	 
	 ++counter;
	 if (counter == MAX_COMPUTE_SUBMITS) break;
	 allTasksCompleted = false;
      }
   } // while (allTasksCompleted == false)
}

void calculateSpatialDerivatives(ParGrid<SpatialCell>& mpiGrid) { }

void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid) { 
   mpiGrid.getCells(Main::cells);

   // Initialization. Only done once per simulation:
   static bool initTransfer = false;
   if (initTransfer == false) {
      // Post first receives for avgs here (needed in spatial fluxes):
      cuint avgsByteSize = mpiGrid[Main::cells[0]]->N_blocks*WID3*sizeof(Real);
      if (mpiGrid.processes() > 1) mpiGrid.startSingleMode();
      for (map<uint,StencilSpatialFluxes>::iterator it=stencilSpatialFluxes.begin(); it!=stencilSpatialFluxes.end(); ++it) it->second.received = 0;
      for (set<uint>::const_iterator it=avgsReceives.begin(); it!=avgsReceives.end(); ++it) {
	 mpiGrid.singleReceive(*it,*it,avgsByteSize,reinterpret_cast<char*>(mpiGrid[*it]->cpu_avgs),*it);
      }
      // Post first sends for avgs:
      for (map<uint,StencilAcceleration>::const_iterator it=stencilAcceleration.begin(); it!=stencilAcceleration.end(); ++it) {
	 for (set<int>::const_iterator dest=it->second.outgoingHosts.begin(); dest!=it->second.outgoingHosts.end(); ++dest) {
	    mpiGrid.singleSend(*dest,it->first,avgsByteSize,reinterpret_cast<char*>(mpiGrid[it->first]->cpu_avgs),it->first);
	 }
      }
      initTransfer = true;
   }

   int counter = 0;
   cuint dynamicMemory = 0;
   dim3 gridSize;
   dim3 blockSize(WID,WID,WID);
   uint cellsToCompute = 0;
   uint computedCells = 0;
   uint cellIndex;
   uint streamID;
   uint priority;
   uint64_t cellID;
   cudaStream_t stream;
   list<uint> remoteToLocalAvgsMap;
   map<uint,RemoteUpdateInfo>::iterator updateIterator;

   cuint fluxByteSize = mpiGrid[Main::cells[0]]->N_blocks*WID3*sizeof(Real);
   
   int myrank = mpiGrid.rank();

   // Post receives for spatial fluxes. These are needed in spatial propagation:
   mpiGrid.startSingleMode2();
   for (map<uint,StencilSpatialProp>::iterator it=stencilSpatialProp.begin(); it!=stencilSpatialProp.end(); ++it) it->second.received = 0;
   
   for (map<pair<uint,int>,uint>::iterator it=spatialFluxReceives.begin(); it!=spatialFluxReceives.end(); ++it) {
      mpiGrid.singleReceive2(it->first.second,it->first.first,fluxByteSize,reinterpret_cast<char*>(fluxReceiveBuffers[it->second]),it->first.first);
   }
   
   // Reset spatial fluxes to zero on all cells (incl. ghosts):
   mpiGrid.getAllCells(Main::cells);
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellIndex = Main::cells[cell];
      Main::cellPtr = mpiGrid[Main::cells[cell]];
      
      cuint offsetFlux = gpu_flux.getOffset(cellIndex);
      Real* gpuFlux    = gpu_flux.getArray<Real>(offsetFlux);

      if (offsetFlux == numeric_limits<uint>::max()) {
	 cerr << "CRITICAL ERROR could not get offsetFlux for cell #" << cellIndex << endl;
      }
      
      gridSize.x = Main::cellPtr->N_blocks;
      gridSize.y = 1;
      cuda_clearFluxes<<<gridSize,blockSize>>>(gpuFlux);
   }
   
   // Initialization. Push all inner cells (=cells with no remote neighbours) to readyCells.
   // Boundary cells (=cells with remote neighbours) are inserted when their neighbour
   // data arrives. Number of cell's remote neighbours is used as its priority.
   for (map<uint,StencilSpatialFluxes>::const_iterator cell = stencilSpatialFluxes.begin(); cell != stencilSpatialFluxes.end(); ++cell) {
      if (cell->second.incoming.size() > 0) continue;
      readyGpuCells.insert(cell->first,cell->second.outgoing.size());
   }
   for (updateIterator=remoteUpdates.begin(); updateIterator!=remoteUpdates.end(); ++updateIterator) updateIterator->second.computed = 0;

   // All local cell df/dt contributions need to be computed:
   cellsToCompute = mpiGrid.getNumberOfLocalCells();

   list<uint> transfersToCpu;
   list<uint> transfersToGpu;
   
   bool allTasksCompleted = false;
   while (allTasksCompleted == false) {
      allTasksCompleted = true;
      // Check for completed GPU -> CPU transfers:
      while (transfersToCpuEvents.getCompleted(streamID,cellID) == true) {
	 // A cell with remote neighbour(s) has been transferred to CPU. Send it to remote host:
       	 int nbrHost;
	 mpiGrid.getHost(cellID,nbrHost);
	 Real* buffer = spatialFluxSends[cellID];
	 
	 mpiGrid.singleSend2(nbrHost,cellID,fluxByteSize,reinterpret_cast<char*>(buffer),cellID);
	 
	 if (streamPool.insert(streamID) == false)
	   cerr << "ERROR occurred while inserting idle streams to streampool!" << endl;
      }
      if (transfersToCpuEvents.empty() == false) allTasksCompleted = false;

      // Check for completed computations:
      while (computationEvents.getCompleted(streamID,cellID) == true) {
	 if (streamPool.insert(streamID) == false)
	   cerr << "ERROR occurred while inserting idle streams to streampool!" << endl;

	 // Increase counter on all remote neighbours that might need the just computed fluxes.
	 // If all local modifications have been made, the data can be transferred to CPU and 
	 // sent over MPI:
	 for (list<uint>::const_iterator it=stencilSpatialFluxes[cellID].outgoing.begin(); it!=stencilSpatialFluxes[cellID].outgoing.end(); ++it) {
	    updateIterator = remoteUpdates.find(*it);
	    ++updateIterator->second.computed;
	    if (updateIterator->second.computed == updateIterator->second.required) {
	       transfersToCpu.push_back(*it);
	    }
	 }
	 ++computedCells;
      }
      if (computedCells != cellsToCompute) allTasksCompleted = false;
	    
      while (transfersToCpu.empty() == false) {
	 if (streamPool.get(streamID,stream) == false) break;
	 cellIndex = transfersToCpu.front();
	 transfersToCpu.pop_front();
	 
	 cuint offsetFlux = gpu_flux.getOffset(cellIndex);
	 Real* gpuFlux    = gpu_flux.getArray<Real>(offsetFlux);
	 Real* cpuFlux    = spatialFluxSends[cellIndex];

	 if (cudaMemcpyAsync(cpuFlux,gpuFlux,fluxByteSize,cudaMemcpyDeviceToHost,stream) != cudaSuccess) 
	   cerr << "Proc #" << mpiGrid.rank() << " failed to transfer remote flux #" << cellIndex << " GPU -> CPU" << endl;
	 
	 transfersToCpuEvents.insert(streamID,stream,cellIndex);
	 allTasksCompleted = false;
      }
      if (transfersToCpu.empty() == false) allTasksCompleted = false;
      
      // Check if streams are available for computing cells:
      counter = 0;
      while (readyGpuCells.empty() == false) {
	 if (streamPool.get(streamID,stream,false) == true) {
	    readyGpuCells.pop(cellIndex,priority);
	    Main::cellPtr = mpiGrid[cellIndex];
	    
	    cuint offsetAvgs    = gpu_avgs.getOffset(cellIndex);
	    Real* gpuAvgs          = gpu_avgs.getArray<Real>(offsetAvgs);
	    
	    // Ghost cells: copy volume averages from a predetermined neighbour. 
	    // Remote cells are left unmodified.
	    if (cellStates[cellIndex].ghostFlag == true && cellStates[cellIndex].remoteFlag == false) {
	       cuint copyNbr = cellStates[cellIndex].copyNbr;
	       if (copyNbr != numeric_limits<uint>::max()) {
		  cuint nbrAvgsOffset = gpu_avgs.getOffset(copyNbr);
		  Real* nbrAvgs       = gpu_avgs.getArray<Real>(nbrAvgsOffset);
		  
		  if (cudaMemcpyAsync(gpuAvgs,nbrAvgs,Main::cellPtr->N_blocks*WID3*sizeof(Real),cudaMemcpyDeviceToDevice,stream) != cudaSuccess)
		    cerr << "ERROR " << cudaGetErrorString(cudaGetLastError()) << " while copying ghost avgs" << endl;
	       }
	    }
	    
	    cuint offsetNbrs        = gpu_nbrsSpa.getOffset(cellIndex);
	    cuint offsetBlockParams = gpu_blockParams.getOffset(cellIndex);
	    cuint offsetCellParams  = gpu_cellParams.getOffset(cellIndex);
	    Real* gpuFlux        = gpu_flux.getArray<Real>();
	          
	    gridSize.x = Main::cellPtr->N_blocks;
	    gridSize.y = 1;	                
	    cuda_trans_fluxes<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuAvgs,gpuFlux,offsetNbrs,offsetCellParams,offsetBlockParams,gpu_dfdt_flags,HALF*Parameters::dt);
	    
	    // Insert computation event to computationEvents:
	    computationEvents.insert(streamID,stream,cellIndex);
	    allTasksCompleted = false;
	 } else break;
	 ++counter;
	 if (counter == MAX_COMPUTE_SUBMITS) break;
      }
      if (computationEvents.empty() == false) allTasksCompleted = false;
      
      // Wait for some remote neighbour data to arrive. 
      // NOTE: This should really be testSomeReceives.
      mpiGrid.singleModeWaitSome();
      
      // Check if enough remote neighbour data has arrived to make 
      // at least one boundary cell ready. NOTE: This is a separate call from 
      // testing/waiting for remote neighbour data because a cell might be pushed 
      // back to ParGrid if enough stream(s) were not available:
      while (mpiGrid.getReadyCell(cellIndex) == true) {
	 transfersToGpu.push_back(cellIndex);
	 allTasksCompleted = false;
      }

      while (transfersToGpu.empty() == false) {
	 // Get a stream for transferring data to GPU.
	 if (streamPool.get(streamID,stream) == false) break;

	 cellIndex = transfersToGpu.front();
	 transfersToGpu.pop_front();
	 
	 cuint offsetAvgs = gpu_avgs.getOffset(cellIndex);
	 Real* gpuAvgs    = gpu_avgs.getArray<Real>(offsetAvgs);
	 Main::cellPtr    = mpiGrid[cellIndex];
	 
	 cudaMemcpyAsync(gpuAvgs,Main::cellPtr->cpu_avgs,Main::cellPtr->N_blocks*WID3*sizeof(Real),cudaMemcpyHostToDevice,stream);
	 transfersToGpuEvents.insert(streamID,stream,cellIndex);
	 allTasksCompleted = false;
      }
      if (transfersToGpu.empty() == false) allTasksCompleted = false;
	 
      // Check for completed CPU -> GPU transfers. If transfer(s) have completed, one or more 
      // local cell may have become ready. In such a case insert the local cell to readyGpuCells.
      while (transfersToGpuEvents.getCompleted(streamID,cellID) == true) {
	 if (streamPool.insert(streamID) == false)
	   cerr << "ERROR occurred while inserting idle streams to streampool!" << endl;
      
	 // Get a list of all local cells that need the transferred remote neighbour data.
	 // Increase arrived counter and test if all required data has been transferred. If 
	 // all data has been transferred, make the cell ready:
	 for (multimap<uint,uint>::const_iterator it=remoteToLocalAvgs.lower_bound(cellID); it!=remoteToLocalAvgs.upper_bound(cellID); ++it) {
	    map<uint,StencilSpatialFluxes>::iterator localCell = stencilSpatialFluxes.find(it->second);
	    ++localCell->second.received;
	    if (localCell->second.received == localCell->second.incoming.size()) {
	       readyGpuCells.insert(it->second,mpiGrid.getNumberOfRemoteNeighbours(it->second));
	    }
	 }
      }
      if (transfersToGpuEvents.empty() == false) allTasksCompleted = false;
   } // while (allTasksCompleted == true)

   mpiGrid.singleModeWaitAllSends(); // wait for MPI sends made in acceleration (avgs)
}

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) { 
   mpiGrid.getCells(Main::cells);

   // Post receives for avgs, needed in fluxes at the beginning of following timestep:
   cuint avgsByteSize = mpiGrid[Main::cells[0]]->N_blocks*WID3*sizeof(Real);
   if (secondStep == true) {
      mpiGrid.startSingleMode();
      for (map<uint,StencilSpatialFluxes>::iterator it=stencilSpatialFluxes.begin(); it!=stencilSpatialFluxes.end(); ++it) it->second.received = 0;
      for (set<uint>::const_iterator it=avgsReceives.begin(); it!=avgsReceives.end(); ++it) {
	 mpiGrid.singleReceive(*it,*it,avgsByteSize,reinterpret_cast<char*>(mpiGrid[*it]->cpu_avgs),*it);
      }
   }
   
   // Clear Velocity Moments:
   dim3 gridSize;
   dim3 blockSize(WID,1,1);
   if (secondStep == true) {
      gridSize.x = Main::cells.size();
      cuda_clearVelocityMoments<<<gridSize,blockSize>>>(gpu_cellParams.getArray<Real>());
   }
   
   // Propagate distribution function in spatial dimensions and 
   // calculate velocity moments:
   blockSize.x = WID;
   blockSize.y = WID;
   blockSize.z = WID;
   
   dim3 redGridSize1(16,1,1);
   dim3 redBlockSize1(32,1,1);
   dim3 redGridSize2(1,1,1);
   dim3 redBlockSize2(16,1,1);

   cuint fluxUpdateSize     = mpiGrid[Main::cells[0]]->N_blocks*WID3;
   cuint fluxUpdateByteSize = mpiGrid[Main::cells[0]]->N_blocks*WID3*sizeof(Real);
   uint priority;
   cudaStream_t stream;
   uint64_t cellID;
   uint streamID;
   uint cellIndex;
   uint computedCells = 0;
   cuint cellsToCompute = mpiGrid.getNumberOfLocalCells();
   cuint dynamicMemory = 0;
   cuint momentBufferSize = 1024;

   int counter = 0;
   bool succeeded;
   uint offsetAvgs  = numeric_limits<uint>::max();
   uint offsetFlux  = numeric_limits<uint>::max();
   uint offsetCellParams = numeric_limits<uint>::max();
   uint offsetBlockParams = numeric_limits<uint>::max();
   uint offsetRho   = numeric_limits<uint>::max();
   uint offsetRhovx = numeric_limits<uint>::max();
   uint offsetRhovy = numeric_limits<uint>::max();
   uint offsetRhovz = numeric_limits<uint>::max();
   uint offsetFluxBuffer = numeric_limits<uint>::max();
   Real* gpuAvgs = NULL;
   Real* gpuFlux = NULL;
   Real* gpuCellParams = NULL;
   Real* gpuRho = NULL;
   Real* gpuRhovx = NULL;
   Real* gpuRhovy = NULL;
   Real* gpuRhovz = NULL;
   Real* gpuFluxBuffer = NULL;
   map<uint,StencilSpatialFluxes>::iterator stencilFluxesIt;
   map<uint,StencilSpatialProp>::iterator stencilSpatialIt;
   set<uint>::const_iterator receiveBufferIt;

   list<uint> transfersToGpu;
   
   // Push all inner cells to readyGpuCells:
   for (map<uint,StencilSpatialProp>::const_iterator it=stencilSpatialProp.begin(); it!=stencilSpatialProp.end(); ++it) {
      if (it->second.incoming.size() == 0) {
	 readyGpuCells.insert(it->first,stencilAcceleration[it->first].outgoing.size());
      }
   }
   
   bool allTasksCompleted = false;
   while (allTasksCompleted == false) {
      allTasksCompleted = true;
      
      // Check if remote flux updates have arrived:
      mpiGrid.singleModeWaitSome2();
      while (mpiGrid.getReadyCell2(cellIndex) == true) {
	 stencilSpatialIt = stencilSpatialProp.find(cellIndex);
	 ++stencilSpatialIt->second.received;
	 if (stencilSpatialIt->second.received == stencilSpatialIt->second.incoming.size()) {

	    // Sum dfdt contributions in CPU into a single array.
	    stencilFluxesIt = stencilSpatialFluxes.find(cellIndex);
	    if (stencilFluxesIt->second.receiveBuffers.size() > 1) {
	       receiveBufferIt = stencilFluxesIt->second.receiveBuffers.begin();
	       Real* transferBuffer = fluxReceiveBuffers[*receiveBufferIt];
	       while (receiveBufferIt != stencilFluxesIt->second.receiveBuffers.end()) {
		  ++receiveBufferIt;
		  Real* sumBuffer = fluxReceiveBuffers[*receiveBufferIt];
		  for (size_t i=0; i<fluxUpdateSize; ++i) transferBuffer[i] += sumBuffer[i];
	       }
	    }
	    transfersToGpu.push_back(cellIndex);
	 }
	 allTasksCompleted = false;
      }
      
      while (transfersToGpu.empty() == false) {
	 if (streamPool.available() == 0) break;
	 
	 cellIndex = transfersToGpu.front();
	 transfersToGpu.pop_front();
	 streamPool.get(streamID,stream); // available() > 0, so this will succeed

	 stencilFluxesIt = stencilSpatialFluxes.find(cellIndex);
	 receiveBufferIt = stencilFluxesIt->second.receiveBuffers.begin();
	 Real* cpuBuffer = fluxReceiveBuffers[*receiveBufferIt];
	 
	 offsetFluxBuffer = gpuFluxBuffers.getOffset(cellIndex);
	 gpuFluxBuffer    = gpuFluxBuffers.getArray<Real>(offsetFluxBuffer);
	 
	 if (cudaMemcpyAsync(gpuFluxBuffer,cpuBuffer,fluxUpdateByteSize,cudaMemcpyHostToDevice,stream) != cudaSuccess) 
	   cerr << "ERROR could not perform async memcpy in spatial fluxes!" << endl;
	 transfersToGpuEvents.insert(streamID,stream,cellIndex);
	 allTasksCompleted = false;
      }
      if (transfersToGpu.empty() == false) allTasksCompleted = false;

      // Check for completed CPU to GPU transfers:
      if (transfersToGpuEvents.getCompleted(streamID,cellID) == true) {
	 streamPool.insert(streamID);
	 readyGpuCells.insert(cellID,stencilAcceleration[cellID].outgoing.size());
	 allTasksCompleted = false;
      }
      if (transfersToGpuEvents.empty() == false) allTasksCompleted = false;

      if (transfersToCpuEvents.getCompleted(streamID,cellID) == true) {
	 streamPool.insert(streamID);
	 for (set<int>::const_iterator it=stencilAcceleration[cellID].outgoingHosts.begin(); it!=stencilAcceleration[cellID].outgoingHosts.end(); ++it) {
	    mpiGrid.singleSend(*it,cellID,avgsByteSize,reinterpret_cast<char*>(mpiGrid[cellID]->cpu_avgs),cellID);
	 }
      }
      if (transfersToCpuEvents.empty() == false) allTasksCompleted = false;
      
      // Check for completed computations:
      if (computationEvents.getCompleted(streamID,cellID) == true) {
	 if (streamPool.insert(streamID) == false)
	   cerr << "ERROR occurred while inserting idle streams to streampool!" << endl;
	 gpuRhoBuffers.releaseArray(cellID);
	 gpuRhovxBuffers.releaseArray(cellID);
	 gpuRhovyBuffers.releaseArray(cellID);
	 gpuRhovzBuffers.releaseArray(cellID);
	 
	 // Divide velocity moments by spatial cell volume in order to get 1/m^3 units:
	 Main::cellPtr = mpiGrid[cellID];
	 creal dV = Main::cellPtr->cpu_cellParams[CellParams::DX]*Main::cellPtr->cpu_cellParams[CellParams::DY]*Main::cellPtr->cpu_cellParams[CellParams::DZ];
	 Main::cellPtr->cpu_cellParams[CellParams::RHO  ] /= dV;
	 Main::cellPtr->cpu_cellParams[CellParams::RHOVX] /= dV;
	 Main::cellPtr->cpu_cellParams[CellParams::RHOVY] /= dV;
	 Main::cellPtr->cpu_cellParams[CellParams::RHOVZ] /= dV;
	 
	 // REMOTE COUNTERS
	 if (secondStep == true) if (stencilAcceleration[cellID].outgoing.size() > 0) {
	    Main::cellPtr = mpiGrid[cellID];
	    offsetAvgs = gpu_avgs.getOffset(cellID);
	    gpuAvgs    = gpu_avgs.getArray<Real>(offsetAvgs);
	    
	    streamPool.get(streamID,stream);
	    cudaMemcpyAsync(Main::cellPtr->cpu_avgs,gpuAvgs,Main::cellPtr->N_blocks*WID3*sizeof(Real),cudaMemcpyDeviceToHost,stream);
	    transfersToCpuEvents.insert(streamID,stream,cellID);
	    allTasksCompleted = false;
	 }
	 ++computedCells;
      }
      if (computedCells != cellsToCompute) allTasksCompleted = false;
      
      // Check for cells that are ready for computing:
      counter = 0;
      while (readyGpuCells.empty() == false) {
	 if (streamPool.get(streamID,stream,false) == false) goto exitComputation; // jump to end of if (readyGpuCells.empty)
	 readyGpuCells.pop(cellIndex,priority);
	 
	 // First check that buffers for reducing velocity moments could be allocated. 
	 // If not, then we need to wait for previous computations to complete.
	 offsetRho   = gpuRhoBuffers.reserveArray(cellIndex,momentBufferSize);
	 offsetRhovx = gpuRhovxBuffers.reserveArray(cellIndex,momentBufferSize);
	 offsetRhovy = gpuRhovyBuffers.reserveArray(cellIndex,momentBufferSize);
	 offsetRhovz = gpuRhovzBuffers.reserveArray(cellIndex,momentBufferSize);
	 
	 succeeded = true;
	 if (offsetRho   == numeric_limits<uint>::max()) succeeded = false;
	 if (offsetRhovx == numeric_limits<uint>::max()) succeeded = false;
	 if (offsetRhovy == numeric_limits<uint>::max()) succeeded = false;
	 if (offsetRhovz == numeric_limits<uint>::max()) succeeded = false;
	 if (succeeded == false) {
	    readyGpuCells.insert(cellIndex,priority);
	    streamPool.insert(streamID);
	    goto exitComputation; // jump to end of if (readyGpuCells.empty)
	 }
	 Main::cellPtr = mpiGrid[cellIndex];
	 gridSize.x = Main::cellPtr->N_blocks;
	 gridSize.y = 1;
	    
	 offsetAvgs        = gpu_avgs.getOffset(cellIndex);
	 offsetFlux        = gpu_flux.getOffset(cellIndex);
	 offsetCellParams  = gpu_cellParams.getOffset(cellIndex);
	 offsetBlockParams = gpu_blockParams.getOffset(cellIndex);
	 offsetFluxBuffer  = gpuFluxBuffers.getOffset(cellIndex);
	 
	 gpuAvgs           = gpu_avgs.getArray<Real>(offsetAvgs);
	 gpuFlux           = gpu_flux.getArray<Real>(offsetFlux);
	 gpuCellParams     = gpu_cellParams.getArray<Real>(offsetCellParams);
	 
	 if (offsetFluxBuffer == numeric_limits<uint>::max()) gpuFluxBuffer = NULL;
	 else gpuFluxBuffer = gpuFluxBuffers.getArray<Real>(offsetFluxBuffer);
	 gpuRho   = gpuRhoBuffers.getArray<Real>(offsetRho);
	 gpuRhovx = gpuRhovxBuffers.getArray<Real>(offsetRhovx);
	 gpuRhovy = gpuRhovyBuffers.getArray<Real>(offsetRhovy);
	 gpuRhovz = gpuRhovzBuffers.getArray<Real>(offsetRhovz);
	    
	 // Ghost cells do not need to be propagated. Velocity moments, however, need to be calculated
	 // on the second propagation step.
	 if (cellStates[cellIndex].ghostFlag == true) {
	    cuda_calcVelocityMoments<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuAvgs,gpuCellParams,offsetBlockParams,gpuRho,gpuRhovx,gpuRhovy,gpuRhovz);
	 } else {
	    if (secondStep == false) {
	       cuda_trans_propagate<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuAvgs,gpuFlux,gpuCellParams,offsetBlockParams,offsetFluxBuffer,gpuFluxBuffer);
	    } else {
	       cuda_trans_propagateWithMoments<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuAvgs,gpuFlux,gpuCellParams,offsetBlockParams,gpuRho,gpuRhovx,gpuRhovy,gpuRhovz,offsetFluxBuffer,gpuFluxBuffer);
	    }
	 }

	 // Calculate velocity moments with reduction if we are on the second propagation step:
	 if (secondStep == true) {
	    cuda_calcVelocityMoments<<<gridSize,blockSize,dynamicMemory,stream>>>(gpuAvgs,gpuCellParams,offsetBlockParams,gpuRho,gpuRhovx,gpuRhovy,gpuRhovz);
	    
	    cuda_trans_reduce1<<<redGridSize1,redBlockSize1,dynamicMemory,stream>>>(gpuRho);
	    cuda_trans_reduce1<<<redGridSize1,redBlockSize1,dynamicMemory,stream>>>(gpuRhovx);
	    cuda_trans_reduce1<<<redGridSize1,redBlockSize1,dynamicMemory,stream>>>(gpuRhovy);
	    cuda_trans_reduce1<<<redGridSize1,redBlockSize1,dynamicMemory,stream>>>(gpuRhovz);
	    
	    cuda_trans_reduce2<<<redGridSize2,redBlockSize2,dynamicMemory,stream>>>(gpuRho  ,gpuCellParams + CellParams::RHO  );
	    cuda_trans_reduce2<<<redGridSize1,redBlockSize2,dynamicMemory,stream>>>(gpuRhovx,gpuCellParams + CellParams::RHOVX);
	    cuda_trans_reduce2<<<redGridSize1,redBlockSize2,dynamicMemory,stream>>>(gpuRhovy,gpuCellParams + CellParams::RHOVY);
	    cuda_trans_reduce2<<<redGridSize1,redBlockSize2,dynamicMemory,stream>>>(gpuRhovz,gpuCellParams + CellParams::RHOVZ);
	    
	    // Velocity moments need to be transferred to CPU on every timestep:
	    if (cudaMemcpyAsync(Main::cellPtr->cpu_cellParams,gpuCellParams,SIZE_CELLPARAMS*sizeof(Real),cudaMemcpyDeviceToHost,stream) != cudaSuccess)
	      cerr << "Failed to insert async memcopy of cell #" << cellIndex << " cellParams!" << endl;
	     
	    // Distribution function needs to be transferred to CPU when data is written to disk:
	    if (transferAvgs == true) // TEST
	      if (cudaMemcpyAsync(Main::cellPtr->cpu_avgs,gpuAvgs,Main::cellPtr->N_blocks*WID3*sizeof(Real),cudaMemcpyDeviceToHost,stream) != cudaSuccess)
		cerr << "Failed to async memcpy avgs of cell #" << cellIndex << " to CPU" << endl;
	 }
	 computationEvents.insert(streamID,stream,cellIndex);
      
	 ++counter;
	 if (counter == MAX_COMPUTE_SUBMITS) break;
	 allTasksCompleted = false;
      }
      exitComputation:
      
   } // while (allTasksCompleted == false)

   if (transfersToCpuEvents.empty() == false) cerr << "Proc #" << mpiGrid.rank() << " has GPU->CPU transfers!" << endl;
   if (transfersToGpuEvents.empty() == false) cerr << "Proc #" << mpiGrid.rank() << " has CPU->GPU transfers!" << endl;
   
   mpiGrid.singleModeWaitAllSends2();
}

void calculateSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, Real& dt) { }
void calculateCellParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, ID::type cell) { }
void initialLoadBalance(ParGrid<SpatialCell>& mpiGrid) { }

void calculateVelocityMoments(ParGrid<SpatialCell>& mpiGrid) { 
   mpiGrid.getCells(Main::cells);
   
   // Clear velocity moments:
   dim3 gridSize;
   dim3 blockSize(WID,1,1);
   gridSize.x = Main::cells.size();
   cuda_clearVelocityMoments<<<gridSize,blockSize>>>(gpu_cellParams.getArray<Real>());

   // Calculate velocity moments:
   blockSize.x = WID;
   blockSize.y = WID;
   blockSize.z = WID;
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellIndex = Main::cells[cell];
      Main::cellPtr   = mpiGrid[cellIndex];
      
      cuint offsetAvgs        = gpu_avgs.getOffset(cellIndex);
      cuint offsetCellParams  = gpu_cellParams.getOffset(cellIndex);
      cuint offsetBlockParams = gpu_blockParams.getOffset(cellIndex);
      Real* gpuAvgs           = gpu_avgs.getArray<Real>(offsetAvgs);
      Real* gpuCellParams     = gpu_cellParams.getArray<Real>(offsetCellParams);

      cuint offsetRho   = gpuRhoBuffers.reserveArray(cellIndex,1024);
      cuint offsetRhovx = gpuRhovxBuffers.reserveArray(cellIndex,1024);
      cuint offsetRhovy = gpuRhovyBuffers.reserveArray(cellIndex,1024);
      cuint offsetRhovz = gpuRhovzBuffers.reserveArray(cellIndex,1024);
      
      Real* gpuRho   = gpuRhoBuffers.getArray<Real>(offsetRho);
      Real* gpuRhovx = gpuRhovxBuffers.getArray<Real>(offsetRhovx);
      Real* gpuRhovy = gpuRhovyBuffers.getArray<Real>(offsetRhovy);
      Real* gpuRhovz = gpuRhovzBuffers.getArray<Real>(offsetRhovz);
      
      gridSize.x = Main::cellPtr->N_blocks;
      gridSize.y = 1;
      cuda_calcVelocityMoments<<<gridSize,blockSize>>>(gpuAvgs,gpuCellParams,offsetBlockParams,gpuRho,gpuRhovx,gpuRhovy,gpuRhovz);
      
      dim3 redGridSize1(16,1,1);
      dim3 redBlockSize1(32,1,1);
      dim3 redGridSize2(1,1,1);
      dim3 redBlockSize2(16,1,1);
      
      cuda_trans_reduce1<<<redGridSize1,redBlockSize1>>>(gpuRho);
      cuda_trans_reduce1<<<redGridSize1,redBlockSize1>>>(gpuRhovx);
      cuda_trans_reduce1<<<redGridSize1,redBlockSize1>>>(gpuRhovy);
      cuda_trans_reduce1<<<redGridSize1,redBlockSize1>>>(gpuRhovz);
      
      cuda_trans_reduce2<<<redGridSize2,redBlockSize2>>>(gpuRho  ,gpuCellParams + CellParams::RHO  );
      cuda_trans_reduce2<<<redGridSize1,redBlockSize2>>>(gpuRhovx,gpuCellParams + CellParams::RHOVX);
      cuda_trans_reduce2<<<redGridSize1,redBlockSize2>>>(gpuRhovy,gpuCellParams + CellParams::RHOVY);
      cuda_trans_reduce2<<<redGridSize1,redBlockSize2>>>(gpuRhovz,gpuCellParams + CellParams::RHOVZ);
      
      gpuRhoBuffers.releaseArray(cellIndex);
      gpuRhovxBuffers.releaseArray(cellIndex);
      gpuRhovyBuffers.releaseArray(cellIndex);
      gpuRhovzBuffers.releaseArray(cellIndex);
      
      if (cudaMemcpy(Main::cellPtr->cpu_cellParams,gpuCellParams,SIZE_CELLPARAMS*sizeof(Real),cudaMemcpyDeviceToHost) != cudaSuccess) 
	 cerr << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while copying cellParams to CPU!" << endl;

      creal dV = Main::cellPtr->cpu_cellParams[CellParams::DX]*Main::cellPtr->cpu_cellParams[CellParams::DY]*Main::cellPtr->cpu_cellParams[CellParams::DZ];
      for (int i=CellParams::RHO; i<SIZE_CELLPARAMS; ++i) Main::cellPtr->cpu_cellParams[i] /= dV;
   }
}

