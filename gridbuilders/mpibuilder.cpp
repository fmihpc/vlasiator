#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <mpi.h>
#include <set>
#include <time.h>
#include "../definitions.h"
#include "../common.h"
#include "../parameters.h"
#include "../mpiconversion.h"
#include "mpibuilder.h"

// Perhaps a better solution is in order for obtaining the 
// initial phase space density, but for now the old file is reused
#include "../project.h"

using namespace std;
namespace VC = VirtualCell;

// *********************************
// ***** CLASS BUFFERALLOCATOR *****
// *********************************

template<typename T> class BufferAllocator {
 public:
   BufferAllocator(MPI_Comm comm,const size_t& N_buffers);
   ~BufferAllocator();
   void freeAll();
   T* getBuffer(const size_t& size);
   bool send(const int& dest,const int& tag);
 private:
   MPI_Request* requests;
   MPI_Comm comm;
   size_t N_buffers;
   size_t* bufferSizes;
   set<size_t> freeBuffers;
   T** buffers;
   size_t currentBuffer;
};

template<typename T>
BufferAllocator<T>::BufferAllocator(MPI_Comm comm,const size_t& N_buffers): comm(comm),N_buffers(N_buffers) {
   currentBuffer = 0;
   requests = new MPI_Request[N_buffers];
   buffers  = new T*[N_buffers];
   bufferSizes = new size_t[N_buffers];
   for (size_t i=0; i<N_buffers; ++i) requests[i] = MPI_REQUEST_NULL;
   for (size_t i=0; i<N_buffers; ++i) freeBuffers.insert(i);
   for (size_t i=0; i<N_buffers; ++i) buffers[i] = NULL;
   for (size_t i=0; i<N_buffers; ++i) bufferSizes[i] = 0;
}

template<typename T>
BufferAllocator<T>::~BufferAllocator() {
   delete [] buffers;
   delete bufferSizes;
   delete requests; 
   buffers = NULL;
   bufferSizes = NULL;
   requests = NULL;
}

template<typename T>
T* BufferAllocator<T>::getBuffer(const size_t& size) {
   // If all buffers are full, wait for some to become available:
   if (freeBuffers.empty() == true) {
      int readySends;
      int* readyIndices = new int[N_buffers];
      MPI_Waitsome(N_buffers,requests,&readySends,readyIndices,MPI_STATUSES_IGNORE);
      for (int rs=0; rs<readySends; ++rs) freeBuffers.insert(readyIndices[rs]);
      delete readyIndices;
   }
   // Remove buffer from list of available buffers:
   set<size_t>::const_iterator it = freeBuffers.begin();
   currentBuffer = *it;
   freeBuffers.erase(currentBuffer);
   
   // Allocate buffer:
   delete buffers[currentBuffer];
   buffers[currentBuffer] = new T[size];
   bufferSizes[currentBuffer] = size;
   return buffers[currentBuffer];
}

template<typename T>
void BufferAllocator<T>::freeAll() {
   MPI_Waitall(N_buffers,requests,MPI_STATUSES_IGNORE);
}

template<typename T>
bool BufferAllocator<T>::send(const int& dest,const int& tag) {
   if (MPI_Isend(buffers[currentBuffer],bufferSizes[currentBuffer],MPI_Type<T>(),dest,tag,comm,&(requests[currentBuffer])) != MPI_SUCCESS) return false;
   return true;
}

// *******************************************
// ***** DEFINITION FOR CLASS MPIBUILDER *****
// *******************************************

MPIBuilder::MPIBuilder(): GridBuilder(),initialized(false) {
   rqstDataCounter = 0;
   N_rqstDataRequests = 0;
   N_rqstDataReceives = 0;
   rqstDataReceives = NULL;
   rqstDataRequests = NULL;
}

MPIBuilder::~MPIBuilder() { 
   delete rqstDataReceives;
   delete rqstDataRequests;
}

// Be careful here, we really can run out of memory in this function....
bool MPIBuilder::addCellBlockDataRequests(VirtualCell::ID& totalCells,VirtualCell::ID* cellIDs,uint* blocksPerCell,
			      Real** avgs,Real** blockParams,uint** nbrsVel) {
   if (initialized == false) return false;
   bool success = true;
   
   // Gather the number of cells requested per process to master:
   VC::ID* cellsPerProcess = NULL;
   if (mpiRank == mpiMasterRank) {
      cellsPerProcess = new VC::ID[N_processes];
   }
   if (MPI_Gather(&totalCells,1,MPI_Type<VC::ID>(),cellsPerProcess,1,MPI_Type<VC::ID>(),mpiMasterRank,comm) != MPI_SUCCESS) success = false;

   // Slave processes post receives for data and exit:
   if (mpiRank != mpiMasterRank) {
      if (rqstDataRequests == NULL) {N_rqstDataRequests = 2;            rqstDataRequests = new MPI_Request[N_rqstDataRequests];}
      if (rqstDataReceives == NULL) {N_rqstDataReceives = 3*totalCells; rqstDataReceives = new MPI_Request[N_rqstDataReceives];}
      MPI_Isend(cellIDs      ,totalCells,MPI_Type<VC::ID>(),mpiMasterRank,mpiRank,comm,&(rqstDataRequests[0]));
      MPI_Isend(blocksPerCell,totalCells,MPI_Type<uint>()  ,mpiMasterRank,mpiRank,comm,&(rqstDataRequests[1]));
      
      // Post separate receive for each spatial cell (required):
      for (VC::ID i=0; i<totalCells; ++i) {
	 MPI_Irecv(avgs[i]       ,blocksPerCell[i]*SIZE_VELBLOCK   ,MPI_Type<Real>(),mpiMasterRank,mpiRank,comm,&(rqstDataReceives[3*i+0]));
	 MPI_Irecv(blockParams[i],blocksPerCell[i]*SIZE_BLOCKPARAMS,MPI_Type<Real>(),mpiMasterRank,mpiRank,comm,&(rqstDataReceives[3*i+1]));
	 MPI_Irecv(nbrsVel[i]    ,blocksPerCell[i]*SIZE_NBRS_VEL   ,MPI_Type<uint>(),mpiMasterRank,mpiRank,comm,&(rqstDataReceives[3*i+2]));
      }
      return success;
   }
   
   VC::ID* cellIDsBuffer      = NULL;
   uint*  blocksPerCellBuffer = NULL;
   BufferAllocator<Real> freeAvgsBuffers(comm,blockBufferSize);
   BufferAllocator<Real> freeParamsBuffers(comm,blockBufferSize);
   BufferAllocator<uint> freeNbrsBuffers(comm,blockBufferSize);
   for (int i=0; i<N_processes; ++i) {
      // Master process copies its own data manually after everyone else has been served,
      // may help with memory management
      if (i == mpiRank) continue;
      
      // Receive cell IDs and blocks per cell from process i:
      delete cellIDsBuffer;
      delete blocksPerCellBuffer;
      try {
	 cellIDsBuffer       = new VC::ID[cellsPerProcess[i]];
	 blocksPerCellBuffer = new uint[cellsPerProcess[i]];
      } catch (exception e) {
	 cerr << "FATAL ERROR: Master process #" << mpiRank << " could not allocate memory for ";
	 cerr << "cellIDsBuffer,blocksPerCellBuffer in addCellBlockDataRequests!" << endl;
	 exit(1);
      }
      if (MPI_Recv(cellIDsBuffer      ,cellsPerProcess[i],MPI_Type<VC::ID>(),i,i,comm,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      if (MPI_Recv(blocksPerCellBuffer,cellsPerProcess[i],MPI_Type<uint>()  ,i,i,comm,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;

      for (VC::ID c=0; c<cellsPerProcess[i]; ++c) {
	 // Get available buffers:
	 Real* const avgsBuffer = freeAvgsBuffers.getBuffer(blocksPerCellBuffer[c]*SIZE_VELBLOCK);
	 Real* const blockParamsBuffer = freeParamsBuffers.getBuffer(blocksPerCellBuffer[c]*SIZE_BLOCKPARAMS);
	 uint* const nbrsVelBuffer = freeNbrsBuffers.getBuffer(blocksPerCellBuffer[c]*SIZE_NBRS_VEL);
	 
	 // Get data to buffers and send:
	 if (getCellBlockData(cellIDsBuffer[c],blocksPerCellBuffer[c],avgsBuffer,blockParamsBuffer,nbrsVelBuffer) == false) success = false;
	 const int dest = i;
	 const int tag  = i;
	 if (freeAvgsBuffers.send(dest,tag)   == false) success = false;
	 if (freeParamsBuffers.send(dest,tag) == false) success = false;
	 if (freeNbrsBuffers.send(dest,tag)   == false) success = false;
      }
   }
   // Wait for the rest of sends to complete and deallocate buffers:
   freeAvgsBuffers.freeAll();
   freeParamsBuffers.freeAll();
   freeNbrsBuffers.freeAll();
   delete cellsPerProcess;
   delete blocksPerCellBuffer;
   delete cellIDsBuffer;
   
   // Master copies its own data directly:
   for (VC::ID c=0; c<totalCells; ++c) {
      if (getCellBlockData(cellIDs[c],blocksPerCell[c],avgs[c],blockParams[c],nbrsVel[c]) == false) success = false;
   }
   return success;
}

bool MPIBuilder::addCellBlockNumberRequests(VirtualCell::ID& totalCells,VirtualCell::ID* cellIDs,uint* N_blocks) {
   if (initialized == false) return false;
   bool success = true;
   
   // Gather the number of cells requested per process to master:
   VC::ID* cellsPerProcess = NULL;
   if (mpiRank == mpiMasterRank) cellsPerProcess = new VC::ID[N_processes];
   if (MPI_Gather(&totalCells,1,MPI_Type<VC::ID>(),cellsPerProcess,1,MPI_Type<VC::ID>(),mpiMasterRank,comm) != MPI_SUCCESS) success = false;
   
   // Slave processes post receives for data and exit:
   if (mpiRank != mpiMasterRank) {
      if (rqstDataRequests == NULL) {N_rqstDataRequests = 1; rqstDataRequests = new MPI_Request[N_rqstDataRequests];}
      if (rqstDataReceives == NULL) {N_rqstDataReceives = 1; rqstDataReceives = new MPI_Request[N_rqstDataReceives];}
      MPI_Isend(cellIDs ,totalCells,MPI_Type<VC::ID>(),mpiMasterRank,mpiRank,comm,&(rqstDataRequests[0]));
      MPI_Irecv(N_blocks,totalCells,MPI_Type<uint>()  ,mpiMasterRank,mpiRank,comm,&(rqstDataReceives[0]));
      return success;
   }
   
   VC::ID* cellIDsBuffer  = NULL;
   BufferAllocator<uint> freeBuffers(comm,sendBufferSize);
   for (int i=0; i<N_processes; ++i) {
      // Master process copies its own data manually
      if (i == mpiRank) {
	 if (getCellNumberOfBlocks(totalCells,cellIDs,N_blocks) == false) success = false;
	 continue;
      }
      
      // Receive cell IDs from process i:
      delete cellIDsBuffer;
      cellIDsBuffer = new VC::ID[cellsPerProcess[i]];
      if (MPI_Recv(cellIDsBuffer,cellsPerProcess[i],MPI_Type<VC::ID>(),i,i,comm,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      uint* const blockBuffer = freeBuffers.getBuffer(cellsPerProcess[i]);
      
      // Get data to buffers and send:
      if (getCellNumberOfBlocks(cellsPerProcess[i],cellIDsBuffer,blockBuffer) == false) success = false;
      const int dest = i;
      const int tag  = i;
      if (freeBuffers.send(dest,tag) == false) success = false;
   }
   // Wait for rest of sends to complete:
   freeBuffers.freeAll();
   delete cellIDsBuffer;
   delete cellsPerProcess;
   return success;
}

/** Send a message to master process to post a request for data for the given cell, 
 * then post receives for the incoming data.
 * @param cellID Global ID of the cell whose data is requested, or 
 * std::numeric_limits<VirtualCell::ID>::max() to indicate that all requests have been added.
 * @param coords Allocated buffer for receiving cell coordinate data. The size of this buffer 
 * must be (number of dimensions) times two.
 * @param nbrIDs Allocated buffer for receiving the global IDs of cell's spatial neighbours.
 * The size of this buffer can be determined from data previously sent by master process.
 * @param nbrTypes Allocated buffer for receiving cell's neighbours' type identifiers. The 
 * size of this buffer must equal the size of nbrIDs.
 * @return If true, the request was posted successfully. If any process returns false at any 
 * point, program execution should be aborted.
 */
bool MPIBuilder::addCellNbrRequests(VirtualCell::ID& totalCells,VirtualCell::ID& totalNbrs,VirtualCell::ID* cellIDs,
				    uchar* nbrsPerCell,Real* coords,VirtualCell::ID* nbrIDs,uchar* nbrTypes) {
   if (initialized == false) return false;
   bool success = true;
   
   // Gather the number of cells requested per process to master:
   VC::ID* cellsPerProcess = NULL;
   VC::ID* nbrsPerProcess = NULL;
   if (mpiRank == mpiMasterRank) {
      cellsPerProcess = new VC::ID[N_processes];
      nbrsPerProcess  = new VC::ID[N_processes];
   }
   if (MPI_Gather(&totalCells,1,MPI_Type<VC::ID>(),cellsPerProcess,1,MPI_Type<VC::ID>(),mpiMasterRank,comm) != MPI_SUCCESS) success = false;
   if (MPI_Gather(&totalNbrs ,1,MPI_Type<VC::ID>(),nbrsPerProcess ,1,MPI_Type<VC::ID>(),mpiMasterRank,comm) != MPI_SUCCESS) success = false;
   
   // Slave processes post receives for data and exit:
   if (mpiRank != mpiMasterRank) {
      if (rqstDataRequests == NULL) {N_rqstDataRequests = 1; rqstDataRequests = new MPI_Request[N_rqstDataRequests];}
      if (rqstDataReceives == NULL) {N_rqstDataReceives = 3; rqstDataReceives = new MPI_Request[N_rqstDataReceives];}
      MPI_Isend(cellIDs    ,totalCells  ,MPI_Type<VC::ID>(),mpiMasterRank,mpiRank,comm,&(rqstDataRequests[0]));
      MPI_Irecv(coords     ,6*totalCells,MPI_Type<Real>()  ,mpiMasterRank,mpiRank,comm,&(rqstDataReceives[0]));
      MPI_Irecv(nbrIDs     ,totalNbrs   ,MPI_Type<VC::ID>(),mpiMasterRank,mpiRank,comm,&(rqstDataReceives[1]));
      MPI_Irecv(nbrTypes   ,totalNbrs   ,MPI_Type<uchar>() ,mpiMasterRank,mpiRank,comm,&(rqstDataReceives[2]));
      return success;
   }
   
   VC::ID* cellIDsBuffer  = NULL;
   BufferAllocator<Real>   freeCoordsBuffers(comm,sendBufferSize);
   BufferAllocator<VC::ID> freeNbrIDsBuffers(comm,sendBufferSize);
   BufferAllocator<uchar>  freeNbrTypesBuffers(comm,sendBufferSize);
   for (int i=0; i<N_processes; ++i) {
      // Master process copies its own data manually
      if (i == mpiRank) {
	 if (getCellNbrData(totalCells,cellIDs,coords,nbrIDs,nbrTypes) == false) success = false;
	 continue;
      }
      // Receive cell IDs from process i:
      delete cellIDsBuffer;
      cellIDsBuffer = new VC::ID[cellsPerProcess[i]];
      if (MPI_Recv(cellIDsBuffer,cellsPerProcess[i],MPI_Type<VC::ID>(),i,i,comm,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      // Get available buffers:
      Real* const coordsBuffer   = freeCoordsBuffers.getBuffer(6*cellsPerProcess[i]);
      VC::ID* const nbrIDsBuffer = freeNbrIDsBuffers.getBuffer(nbrsPerProcess[i]);
      uchar* const nbrTypesBuffer = freeNbrTypesBuffers.getBuffer(nbrsPerProcess[i]);
      
      // Get data to buffers and send:
      if (getCellNbrData(cellsPerProcess[i],cellIDsBuffer,coordsBuffer,nbrIDsBuffer,nbrTypesBuffer) == false) success = false;
      const int dest = i;
      const int tag  = i;
      if (freeCoordsBuffers.send(dest,tag) == false) success = false;
      if (freeNbrIDsBuffers.send(dest,tag) == false) success = false;
      if (freeNbrTypesBuffers.send(dest,tag) == false) success = false;
   }
   // Wait for rest of sends to complete:
   freeCoordsBuffers.freeAll();
   freeNbrIDsBuffers.freeAll();
   freeNbrTypesBuffers.freeAll();
   delete cellsPerProcess;
   delete nbrsPerProcess;
   delete cellIDsBuffer;
   return success;
}

bool MPIBuilder::addCellParamsRequests(VirtualCell::ID& totalCells,VirtualCell::ID* cellIDs,Real* cellParams) {
   bool success = true;
   if (initialized == false) return false;

   // Gather the number of cells requested per process to master:
   VC::ID* cellsPerProcess = NULL;
   if (mpiRank == mpiMasterRank) cellsPerProcess = new VC::ID[N_processes];   
   if (MPI_Gather(&totalCells,1,MPI_Type<VC::ID>(),cellsPerProcess,1,MPI_Type<VC::ID>(),mpiMasterRank,comm) != MPI_SUCCESS) success = false;
   
   // Slave processes post receives for data and exit:
   if (mpiRank != mpiMasterRank) {
      if (rqstDataRequests == NULL) {N_rqstDataRequests = 1; rqstDataRequests = new MPI_Request[N_rqstDataRequests];}
      if (rqstDataReceives == NULL) {N_rqstDataReceives = 1; rqstDataReceives = new MPI_Request[N_rqstDataReceives];}
      MPI_Isend(cellIDs    ,totalCells                ,MPI_Type<VC::ID>(),mpiMasterRank,mpiRank,comm,&(rqstDataRequests[0]));
      MPI_Irecv(cellParams ,totalCells*SIZE_CELLPARAMS,MPI_Type<Real>()  ,mpiMasterRank,mpiRank,comm,&(rqstDataReceives[0]));
      return success;
   }
   
   VC::ID* cellIDsBuffer  = NULL;
   BufferAllocator<Real> freeBuffers(comm,sendBufferSize);
   for (int i=0; i<N_processes; ++i) {
      // Master process copies its own data manually
      if (i == mpiRank) {
	 if (getCellParams(totalCells,cellIDs,cellParams) == false) success = false;
	 continue;
      }
      
      // Receive cell IDs from process i:
      delete cellIDsBuffer;
      cellIDsBuffer = new VC::ID[cellsPerProcess[i]];
      if (MPI_Recv(cellIDsBuffer,cellsPerProcess[i],MPI_Type<VC::ID>(),i,i,comm,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;

      // Get a free buffer for sending data to process i:
      Real* const paramsBuffer = freeBuffers.getBuffer(cellsPerProcess[i]*SIZE_CELLPARAMS);
      
      // Get data to buffers and send:
      if (getCellParams(cellsPerProcess[i],cellIDsBuffer,paramsBuffer) == false) success = false;
      const int dest = i;
      const int tag  = i;
      if (freeBuffers.send(dest,tag) == false) success = false;
   }
   // Wait for rest of sends to complete:
   freeBuffers.freeAll();
   delete cellsPerProcess;
   delete cellIDsBuffer;
   return success;
}

bool MPIBuilder::finalize() {return true;}

bool MPIBuilder::initialize(MPI_Comm comm,const int& MASTER_RANK) {
   typedef Parameters P;
   initialized = true;

   MPI_Comm_rank(comm,&mpiRank);
   MPI_Comm_size(comm,&N_processes);
   mpiMasterRank = MASTER_RANK;
   this->comm = comm;
   
   // Only master process reads parameters:
   if (mpiRank == mpiMasterRank) {
      // Define the required parameters and query their values from Parameters:
      P::add("mpibuilder.send_buffer_size","Size of send buffer (measured in number of cells), defaults to value 100.",sendBufferSize,100);
      P::add("mpibuilder.block_buffer_size","Size of send buffer (measured in number of spatial cells) when distributing velocity blocks, defaults to value 10.",blockBufferSize,10);
      P::add("gridbuilder.q","Charge of simulated particle species, in Coulombs.",q,numeric_limits<Real>::max());
      P::add("gridbuilder.m","Mass of simulated particle species, in kilograms.",m,numeric_limits<Real>::max());
      P::add("gridbuilder.dt","Timestep in seconds.",dt,numeric_limits<Real>::max());
      P::add("gridbuilder.t_min","Simulation time at timestep 0, in seconds.",t_min,numeric_limits<Real>::max());
      P::add("gridbuilder.timestep","Timestep when grid is loaded. Defaults to value zero.",timestep,0);
      P::add("gridbuilder.max_timesteps","Max. value for timesteps. Defaults to value zero.",max_timesteps,0);
      P::parse();
      
      // Check that we got sane values:
      if (q == numeric_limits<Real>::max()) initialized = false;
      if (m == numeric_limits<Real>::max()) initialized = false;
      if (dt == numeric_limits<Real>::max()) initialized = false;
      if (t_min == numeric_limits<Real>::max()) initialized = false;
   }
   
   // Master process lets everyone know if everything is ok:
   bool globalSuccess = initialized;
   if (MPI_Bcast(&globalSuccess,1,MPI_BYTE,MASTER_RANK,comm) != MPI_SUCCESS) initialized = false;
   if (globalSuccess == false) initialized = false;
   if (initialized == false) return initialized;
   
   // Broadcast the values of internal variables, if needed:
   if (MPI_Bcast(&sendBufferSize,1,MPI_Type<uint>(),MASTER_RANK,comm) != MPI_SUCCESS) initialized = false;
   return initialized;
}

bool MPIBuilder::processCellBlockDataRequests() {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return true;
   bool success = true;
   return success;
}

bool MPIBuilder::processCellBlockNumberRequests() {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return true;
   bool success = true;
   return success;
}

bool MPIBuilder::processCellNbrRequests() {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return true;
   bool success = true;
   return success;
}

bool MPIBuilder::processCellParamsRequests() {
   if (initialized == false) return false;
   if (mpiRank != mpiMasterRank) return true;
   bool success = true;
   return success;
}

bool MPIBuilder::waitCellBlockDataRequests() {return waitCellNbrRequests();}

bool MPIBuilder::waitCellBlockNumberRequests() {return waitCellNbrRequests();}

bool MPIBuilder::waitCellNbrRequests() {
   bool success = true;
   
   if (mpiRank != mpiMasterRank) {
      // wait for all requests for neighbour data to complete:
      if (MPI_Waitall(N_rqstDataRequests,rqstDataRequests,MPI_STATUSES_IGNORE) != MPI_SUCCESS) success = true;
   
      // Wait for all data receive requests to complete:
      mpiStatuses.resize(N_rqstDataReceives);
      if (MPI_Waitall(N_rqstDataReceives,rqstDataReceives,&(mpiStatuses[0])) != MPI_SUCCESS) success = true;
   }
   // Free memory and return:
   rqstDataCounter = 0;
   N_rqstDataRequests = 0;
   N_rqstDataReceives = 0;
   mpiStatuses.resize(0);

   delete rqstDataRequests;
   rqstDataRequests = NULL;
   delete rqstDataReceives;
   rqstDataReceives = NULL;
   return success;
}

bool MPIBuilder::waitCellParamsRequests() {return waitCellNbrRequests();}

