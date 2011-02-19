#ifndef MPIBUILDER_H
#define MPIBUILDER_H

#include <vector>
#include <map>
#include "../gridbuilder.h"

class MPIBuilder: public GridBuilder {
 public:
   MPIBuilder();
   virtual ~MPIBuilder();

   // **************************************************
   // ***** VIRTUAL FUNCTIONS OF CLASS GRIDBUILDER *****
   // *****        DEFINED BY MPIBUILDER.          *****
   // **************************************************
   
   virtual bool addCellBlockDataRequests(VirtualCell::ID& totalCells,VirtualCell::ID* cellIDs,uint* blocksPerCell,
					 Real** avgsBuffer,Real** blockParamsBuffer,uint** nbrsVelBuffer);
   virtual bool addCellBlockNumberRequests(VirtualCell::ID& totalCells,VirtualCell::ID* cellIDs,uint* N_blocks);
   virtual bool addCellNbrRequests(VirtualCell::ID& totalCells,VirtualCell::ID& totalNbrs,VirtualCell::ID* cellIDs,
				   uchar* nbrsPerCell,Real* coords,VirtualCell::ID* nbrIDs,uchar* nbrTypes);
   virtual bool addCellParamsRequests(VirtualCell::ID& totalCells,VirtualCell::ID* cellIDs,Real* cellParams);
   virtual bool finalize();
   virtual bool initialize(MPI_Comm comm,const int& MASTER_RANK);
   virtual bool processCellBlockDataRequests();
   virtual bool processCellBlockNumberRequests();
   virtual bool processCellNbrRequests();
   virtual bool processCellParamsRequests();
   virtual bool waitCellBlockDataRequests();
   virtual bool waitCellBlockNumberRequests();
   virtual bool waitCellNbrRequests();
   virtual bool waitCellParamsRequests();
   
   // **************************************************
   // ***** DECLARATIONS OF PURE VIRTUAL FUNCTIONS *****
   // *****  IMPLEMENTING CLASS MUST DEFINE THESE  *****
   // **************************************************

   virtual bool getCellBlockData(const VirtualCell::ID& cellID,cuint& N_blocks,Real* avgs,Real* blockParams,uint* nbrsVel) = 0;
   virtual bool getCellNumberOfBlocks(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,uint* N_blocks) = 0;
   virtual bool getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams) = 0;
   
   // **************************************************
   // ***** VIRTUAL FUNCTIONS OF CLASS GRIDBUILDER *****
   // *****   THAT ARE NOT DEFINED BY MPIBUILDER.  *****
   // *****  IMPLEMENTING CLASS MUST DEFINE THESE  *****
   // **************************************************
   
   //virtual bool getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs) = 0;
   //virtual bool getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* nbrTypes) = 0;
   //virtual bool getParameter(const std::string& parameterName,std::string& value) = 0;
   //virtual bool getTotalNumberOfCells(VirtualCell::ID& N_cells) = 0;
   
 protected:
   MPI_Comm comm;                         /**< MPI Communicator used by this process.*/
   bool initialized;                      /**< If true, MPIBuilder initialized successfully.*/
   int mpiRank;                           /**< MPI rank of this process in the given communicator.*/
   int mpiMasterRank;                     /**< MPI rank of master process in MPI communicator comm.*/
   int N_processes;                       /**< Number of processes in MPI communicator comm.*/
   uint sendBufferSize;                   /**< Size of send buffer (in number of spatial cells, not bytes).*/
   uint blockBufferSize;                  /**< Size of send buffer (in number of spatial cells, not bytes) when 
					   * sending velocity block data to slave processes.*/
   std::vector<MPI_Status>  mpiStatuses;  /**< MPI_Status messages tracked when processing requests.*/
   
   int rqstDataCounter;                   /**< Counter tracking the number of placed data requests.*/
   MPI_Request* rqstDataReceives;         /**< MPI_Request messages related to receiving the requested data.*/
   MPI_Request* rqstDataRequests;         /**< MPI_Request messages related to sending requests of data.*/
   uint N_rqstDataReceives;               /**< Size of rqstDataReceives.*/
   uint N_rqstDataRequests;               /**< Size of rqstDataRequests.*/
};

#endif
