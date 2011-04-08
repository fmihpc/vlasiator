#ifndef PARGRID_H
#define PARGRID_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <limits>
#include <list>
#include <vector>
#include <mpi.h>
#include <zoltan_cpp.h>
#include <ctime>
#include <set>
#include <cmath>

#include "definitions.h"
#include "parameters.h"
#include "mpiconversion.h"
#include "mpilogger.h"
extern MPILogger mpilogger;

namespace ID {
   typedef unsigned int type;
}

// Partitioners which can be used.
enum LBM {
   Block,Random,RCB,RIB,HSFC,Graph,Hypergraph,Hierarchical
};

template<class C> struct ParCell {
   uint unrefInd_i; /**< The i index of the cell.*/
   uint unrefInd_j; /**< The j index of the cell.*/
   uint unrefInd_k; /**< The k index of the cell.*/
   Real xcrd;       /**< The x-coordinate of the midpoint of the cell.*/
   Real ycrd;       /**< The y-coordinate of the midpoint of the cell.*/
   Real zcrd;       /**< The z-coordinate of the midpoint of the cell.*/
   Real dx;
   Real dy;
   Real dz;
   uint N_blocks;   /**< Number of velocity blocks in the cell.*/
   
   bool hasRemoteNeighbours; /**< If true, at least one of the cell's neighbours are 
			      * assigned to a different MPI process.*/
   uint N_remNbrs;           /**< Total number of remote neighbours this cell has.*/
   uint N_receivedRemNbrs;   /**< Number of this cell's remote neighbours that have been received from other MPI processes.*/
   std::map<uchar,ID::type> neighbours; /**< Global IDs of cell's neighbours. uchar is a user-defined number
					 * identifying which neighbour, e.g. +x or -y, each one is.*/
   
   unsigned char refLevel;   /**< Refinement level of the cell, 0 = base grid (unrefined).*/
   C* dataptr;               /**< Pointer to user data.*/
   ParCell(): dataptr(NULL),N_blocks(0),hasRemoteNeighbours(false),refLevel(0) { }
   ~ParCell() {delete dataptr; dataptr=NULL;}
};

template<class C> class ParGrid {
 public: 
   
   ParGrid(const LBM& method,int argn,char* args[]);
   ~ParGrid();

   // Some functions which have the same name as in dccrg:
   Real get_cell_x_size(const ID::type& globalID) const;
   Real get_cell_y_size(const ID::type& globalID) const;
   Real get_cell_z_size(const ID::type& globalID) const;
   Real get_cell_x(const ID::type& globalID) const;
   Real get_cell_y(const ID::type& globalID) const;
   Real get_cell_z(const ID::type& globalID) const;
   Real get_cell_x_max(const ID::type& globalID) const;
   Real get_cell_x_min(const ID::type& globalID) const;
   Real get_cell_y_max(const ID::type& globalID) const;
   Real get_cell_y_min(const ID::type& globalID) const;
   Real get_cell_z_max(const ID::type& globalID) const;
   Real get_cell_z_min(const ID::type& globalID) const;
   ID::type get_cell(const Real& x,const Real& y,const Real& z) const;

   template<typename T> bool addCell(const T& cellID,const Real* const coords,const Real* const sizes,
				     const std::vector<T>& nbrIDs,const std::vector<uchar>& nbrTypes);
   void barrier() {MPI_Barrier(MPI_COMM_WORLD);}
   template<class CONT> void getAllCells(CONT& rlist) const;
   template<class CONT> void getBoundaryCells(CONT& rlist) const;
   void getCellDistribution(std::map<ID::type,int>& rlist) const {rlist=hostProcesses;}
   template<class CONT> void getCells(CONT& rlist) const;
   template<class CONT> void getExistingNeighbours(CONT& rlist,const ID::type& globalID) const;
   bool getHost(const ID::type& cellID,int& host) const;
   template<class CONT> void getInnerCells(CONT& rlist) const;
   ID::type getNeighbour(const ID::type globalID,cuchar& nbrTypeID) const;
   uint getNumberOfLocalCells() const {return localCells.size();}
   uint getNumberOfReceives() const {return receiveList.size();}
   uint getNumberOfRemoteCells() const {return remoteCells.size();}
   uint getNumberOfRemoteNeighbours(const ID::type& globalID) const; 
   uint getNumberOfSends() const {return sendList.size();}
   bool getReadyCell(ID::type& globalID);
   bool getReadyCell2(ID::type& globalID);
   template<class CONT> void getReceiveList(CONT& rlist) const;
   unsigned char getRefinementLevel(const ID::type& globalID) const;
   uint getRemainingReceives() const;
   uint getRemainingReceives2() const;
   template<class CONT> void getRemoteCells(CONT& rlist) const;
   ID::type getRemoteNeighbour(const ID::type& globalID,const uchar& nbrTypeID) const;
   template<class CONT> void getRemoteNeighbours(const ID::type& globalID,CONT& rlist) const;
   template<class CONT> void getRemoteToLocalMapping(const ID::type& remoteID,CONT& rlist) const;
   template<class CONT> void getSendList(CONT& rlist) const;
   bool hasRemoteNeighbours(const ID::type& globalID) const;
   bool initialize();
   bool initialize(cuint& xsize,cuint& ysize,cuint& zsize,creal& xmin,creal& ymin,creal& zmin,
		   creal& xmax,creal& ymax,creal& zmax);
   bool initialLoadBalance();
   bool isInitialized() const {return initialized;}
   bool loadBalance();
   C* operator[](const ID::type& id) const;
   void print() const;
   void printTransfers() const;
   int processes() const;
   int rank() const;
   bool returnReadyCell(const ID::type& globalID);
   bool returnReadyCell2(const ID::type& globalID);
   bool setLoadBalancingMethod(const LBM& method);
   uint singleModeTestSome();
   uint singleModeTestSome2();
   bool singleModeWaitAllSends();
   bool singleModeWaitAllSends2();
   uint singleModeWaitSome();
   uint singleModeWaitSome2();
   bool singleReceive(const ID::type& sourceID,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID);
   bool singleReceive2(const int& hostID,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID);
   bool singleSend(const int& destHost,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID);
   bool singleSend2(const int& destHost,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID);
   bool startNeighbourExchange(const uint& identifier);
   void startSingleMode(const int& transfers=0);
   void startSingleMode2(const int& transfers=0);
   bool uncalculatedCells() const;
   bool waitAll();
   bool waitAllReceives();
   bool waitAllSends();
   bool waitAnyReceive();
   void waitForReceives() const;
   bool waitSomeReceives();
   void writeLoadDistribution();

   // Zoltan callback functions
   static void getCellCoordinates(void* data,int N_globalEntries,int N_localEntries,ZOLTAN_ID_PTR globalID,
				  ZOLTAN_ID_PTR localID,double* geometryData,int* rcode);
   static void getLocalCellList(void* data,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
				ZOLTAN_ID_PTR localIDs,int weightDim,float* cellWeights,int* rcode);
   static int getGridDimensions(void* data,int* rcode);
   static void getEdgeList(void* parGridPtr,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			   ZOLTAN_ID_PTR localID,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			   int N_weights,float* weight,int* rcode);
   static void getHierarchicalParameters(void* parGridPtr,int level,Zoltan_Struct* zs,int* rcode);
   static int getHierarchicalPartNumber(void* parGridPtr,int level,int* rcode);
   static void getHyperedges(void* parGridPtr,int N_globalIDs,int N_vtxedges,int N_pins,int format,
			     ZOLTAN_ID_PTR vtxedge_GID,int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* rcode);
   static void getHyperedgeWeights(void* parGridPtr,int N_globalIDs,int N_localIDs,int N_edges,int N_weights,
				   ZOLTAN_ID_PTR edgeGlobalID,ZOLTAN_ID_PTR edgeLocalID,float* edgeWeights,int* rcode);
   static int getNumberOfEdges(void* parGridPtr,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
			       ZOLTAN_ID_PTR localIDs,int* rcode);
   static int getNumberOfHierarchicalLevels(void* parGridPtr,int* rcode);
   static void getNumberOfHyperedges(void* parGridPtr,int* N_lists,int* N_pins,int* format,int* rcode);
   static void getNumberOfHyperedgeWeights(void* parGridPtr,int* N_edges,int* rcode);
   static int getNumberOfLocalCells(void* data,int* rcode);

 private:
   static LBM balanceMethod;          /**< The load balance method currently in use.*/
   Real grid_xmin;             /**< x-coordinate of the lower left corner of the grid.*/
   Real grid_xmax;             /**< x-coordinate of the upper right corner of the grid.*/
   Real grid_ymin;             /**< y-coordinate of the lower left corner of the grid.*/
   Real grid_ymax;             /**< y-coordinate of the upper right corner of the grid.*/
   Real grid_zmin;             /**< z-coordinate of the lower left corner of the grid.*/
   Real grid_zmax;             /**< z-coordinate of the upper right corner of the grid.*/
   static std::string imbalanceTolerance; /**< Imbalance tolerance of the load balancing.*/
   bool initialized;           /**< If true, ParGrid was initialized successfully and is ready for use.*/
   MPI_Datatype MPIdataType;   /**< The MPI datatype which is currently used in communications.*/
   bool MPItypeFreed;          /**< If true, MPIdataType has been deallocated and it is safe to create a new one.*/
   static int myrank;                 /**< The rank if this MPI process.*/
   static uint N_hierarchicalLevels;
   int N_processes;            /**< Total number of MPI processes.*/
   static uint N_processesPerPart;
   uint N_receivesRemaining;   /**< Number of remote cells this process has received.*/
   uint N_receivesRemaining2;
   std::string N_weights_cell; /**< Number of weights assigned to each cell.*/
   std::string N_weights_edge; /**< Number of weights assigned to each edge.*/
   bool periodic_x;            /**< If true, x-direction is periodic.*/
   bool periodic_y;            /**< If true, y-direction is periodic.*/
   bool periodic_z;            /**< If true, z-direction is periodic.*/
   Real unref_dx;              /**< Unrefined cell size in x-direction.*/
   Real unref_dy;              /**< Unrefined cell size in y-direction.*/
   Real unref_dz;              /**< Unrefined cell size in z-direction.*/
   uint unrefSize_x;           /**< Number of x-cells in unrefined grid.*/
   uint unrefSize_y;           /**< Number of y-cells in unrefined grid.*/
   uint unrefSize_z;           /**< Number of z-cells in unrefined grid.*/
   static float weightCell;    /**< Weight of each cell (Zoltan).*/
   static float weightEdge;    /**< Weight of each edge (Zoltan).*/
   Zoltan* zoltan;             /**< Pointer to Zoltan load balancing library.*/
   
   static std::map<ID::type,int> hostProcesses;            /**< For each cell, the host process that contains that cell.
							    * The contents of this array should be the same on each MPI
							    * process.*/
   static std::map<ID::type,ParCell<C> > localCells;       /**< Associative container containing all cells that are currently 
							    * assigned to this process.*/
   std::vector<uint> localReceiveIDs;
   std::vector<uint> localReceiveIDs2;
   std::map<ID::type,uint> nbrReferences;
   std::multimap<ID::type,ID::type> remoteToLocal;
   std::vector<MPI_Request> MPIrecvRequests;               /**< Container for active MPI_Requests due to receives.*/
   std::vector<MPI_Request> MPIrecvRequests2;
   std::vector<MPI_Request> MPIsendRequests;               /**< Container for active MPI_Requests due to sends.*/
   std::vector<MPI_Request> MPIsendRequests2;
   std::vector<MPI_Status> MPIstatuses;                    /**< Container for holding MPI_Status messages from sends and receives.*/
   std::list<ID::type> readyCells;
   std::list<ID::type> readyCells2;
   static std::map<ID::type,int> receiveList;              /**< A list of cells, identified by their global ID, that are 
							    * sent to this process during neighbour data exchange, and the 
							    * rank of the MPI process that sends the cell.*/
   static std::map<ID::type,ParCell<C> > remoteCells;
   static std::map<std::pair<ID::type,int>,char> sendList; /**< A list of cells, identified by their global ID, that this 
							    * MPI process has to send during neighbour data exchange, and 
							    * the rank of the receiving MPI process.*/

   void allocateCells();
   void buildExchangeLists();
   void buildInitialGrid();
   void buildUnrefNeighbourLists();
   static float calculateCellWeight(const ID::type& globalID);
   static float calculateEdgeWeight(const ID::type& globalID);
   static float calculateHyperedgeWeight(const ID::type& globalID);
   uint calculateNbrIndex(const ID::type& globalID,const ID::type& nbrID);
   ID::type calculateUnrefinedIndex(const ID::type& i,const ID::type& j,const ID::type& k) const;
   void calculateUnrefinedIndices(const ID::type& index,ID::type& i,ID::type& j,ID::type& k) const;
   static std::string loadBalanceMethod(const LBM& method); // TEST: const removed
   bool syncCellAssignments();
   void syncCellCoordinates();
};

// ***************************************************************
// ************** DEFINITIONS FOR STATIC MEMBERS *****************
// ***************************************************************
template<class C> float ParGrid<C>::weightCell;
template<class C> float ParGrid<C>::weightEdge;
template<class C> std::map<ID::type,int> ParGrid<C>::hostProcesses;
template<class C> std::map<ID::type,ParCell<C> > ParGrid<C>::localCells;
template<class C> std::map<ID::type,int> ParGrid<C>::receiveList;
template<class C> std::map<ID::type,ParCell<C> > ParGrid<C>::remoteCells;
template<class C> std::map<std::pair<ID::type,int>,char> ParGrid<C>::sendList;

template<class C> uint ParGrid<C>::N_processesPerPart;
template<class C> uint ParGrid<C>::N_hierarchicalLevels;
template<class C> int ParGrid<C>::myrank;
template<class C> std::string ParGrid<C>::imbalanceTolerance;
template<class C> LBM ParGrid<C>::balanceMethod;

// ***************************************************************
// ************** BEGIN MEMBER FUNCTION DEFITINIONS **************
// ***************************************************************
template<class C> ParGrid<C>::ParGrid(const LBM& method,int argn,char* args[]) {
   initialized = true;
   unref_dx = 0.0;
   unref_dy = 0.0;
   unref_dz = 0.0;
   MPItypeFreed = true;
   
   // Set load balancing parameters:
   balanceMethod = method;
   N_weights_cell = "1";
   N_weights_edge = "0";
   imbalanceTolerance = "1.05";
   weightCell = 1.0;
   weightEdge = 10.0;
   N_hierarchicalLevels = 2;
   N_processesPerPart = 12;

   // Get the rank of this process, and the total number of MPI processes:
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
   // Attempt to init Zoltan:
   float zoltanVersion;
   int rvalue = Zoltan_Initialize(argn,args,&zoltanVersion);
   if (rvalue != ZOLTAN_OK) {
      std::cerr << "ParGrid: Zoltan init failed!" << std::endl;
      initialized = false;
   }
   if (initialized == false) return;
   
   // Create a new Zoltan object:
   zoltan = new Zoltan(MPI_COMM_WORLD);
   
   // Set Zoltan parameters:
   zoltan->Set_Param("NUM_GID_ENTRIES","1");
   zoltan->Set_Param("LB_METHOD",loadBalanceMethod(balanceMethod).c_str());
   zoltan->Set_Param("LB_APPROACH","PARTITION");
   zoltan->Set_Param("RETURN_LISTS","ALL");
   zoltan->Set_Param("OBJ_WEIGHT_DIM",N_weights_cell.c_str());
   zoltan->Set_Param("EDGE_WEIGHT_DIM",N_weights_edge.c_str());
   zoltan->Set_Param("DEBUG_LEVEL","0");
   zoltan->Set_Param("IMBALANCE_TOL",imbalanceTolerance.c_str());
   zoltan->Set_Param("HIER_CHECKS","0");
   zoltan->Set_Param("HIER_DEBUG_LEVEL","0");
   zoltan->Set_Param("PHG_CUT_OBJECTIVE","CONNECTIVITY");
   
   // Register generic callback functions:
   zoltan->Set_Num_Obj_Fn(&ParGrid<C>::getNumberOfLocalCells,this);
   zoltan->Set_Obj_List_Fn(&ParGrid<C>::getLocalCellList,this);   
   // Register geometry-based load balancing callback functions:
   zoltan->Set_Num_Geom_Fn(&ParGrid<C>::getGridDimensions,this);
   zoltan->Set_Geom_Fn(&ParGrid<C>::getCellCoordinates,this);
   // Register graph-based load balancing callback functions:
   zoltan->Set_Num_Edges_Fn(&ParGrid<C>::getNumberOfEdges,this);
   zoltan->Set_Edge_List_Fn(&ParGrid<C>::getEdgeList,this);
   // Register hypergraph-based load balancing callback functions:
   zoltan->Set_HG_Size_CS_Fn(getNumberOfHyperedges,this);
   zoltan->Set_HG_CS_Fn(getHyperedges,this);
   zoltan->Set_HG_Size_Edge_Wts_Fn(getNumberOfHyperedgeWeights,this);
   zoltan->Set_HG_Edge_Wts_Fn(getHyperedgeWeights,this);
   // Register hierarchical load balancing callback functions:
   zoltan->Set_Hier_Num_Levels_Fn(getNumberOfHierarchicalLevels,this);
   zoltan->Set_Hier_Part_Fn(getHierarchicalPartNumber,this);
   zoltan->Set_Hier_Method_Fn(getHierarchicalParameters,this);
}

template<class C> ParGrid<C>::~ParGrid() {
   // Call destructor for each user-defined data cell:
   typename std::map<ID::type,ParCell<C> >::iterator it = localCells.begin();
   while (it != localCells.end()) {
      delete it->second.dataptr;
      it->second.dataptr = NULL;
      ++it;
   }   
   it = remoteCells.begin();
   while (it != remoteCells.end()) {
      delete it->second.dataptr;
      it->second.dataptr = NULL;
      ++it;
   }
   
   // Delete Zoltan before calling MPI_Finalize():
   delete zoltan;
   zoltan = NULL;
   
   MPI_Finalize();
}

template<class C> template<typename T>
bool ParGrid<C>::addCell(const T& cellID,const Real* const coords,const Real* sizes,
			 const std::vector<T>& nbrIDs,const std::vector<uchar>& nbrTypes) {
   /*
   std::cerr << "ParGrid: proc #" << myrank << " received cell #" << cellID << ' ';
   for (int i=0; i<3; ++i) std::cerr << coords[i] << ' ';
   for (int i=0; i<3; ++i) std::cerr << sizes[i] << ' ';   
   for (size_t i=0; i<nbrIDs.size(); ++i) std::cerr << nbrIDs[i] << ":" << (int)nbrTypes[i] << ' ';
   std::cerr << std::endl;
   */
   localCells[cellID];
   localCells[cellID].xcrd = coords[0];
   localCells[cellID].ycrd = coords[1];
   localCells[cellID].zcrd = coords[2];
   localCells[cellID].dx = sizes[0];
   localCells[cellID].dy = sizes[1];
   localCells[cellID].dz = sizes[2];
   for (size_t i=0; i<nbrIDs.size(); ++i) {
      localCells[cellID].neighbours[nbrTypes[i]] = nbrIDs[i];
      if (nbrReferences.find(nbrIDs[i]) == nbrReferences.end()) nbrReferences[nbrIDs[i]] = 0;
      ++nbrReferences[nbrIDs[i]];
   }
   return true;
}

/** Call the constructor of each user-defined data cell. Here we also allocate 
 * memory for remote cells.
 */
template<class C> void ParGrid<C>::allocateCells() {
   // Local cells have been pushed to localCells, but remoteCells is empty. Insert remote cells.
   for (std::map<ID::type,int>::const_iterator it=receiveList.begin(); it!=receiveList.end(); ++it) {
      remoteCells[it->first];
   }
   #ifndef NDEBUG
      std::cout << "ParGrid::allocateCells for process " << myrank << " requires " << localCells.size() << " local cells, ";
      std::cout << remoteCells.size() << " remote cells." << std::endl;
   #endif
   
   // Now ask user to allocate memory for each cell (if necessary):
   for (typename std::map<ID::type,ParCell<C> >::iterator it=localCells.begin(); it!=localCells.end(); ++it) {
      it->second.dataptr = new C;
   }
   for (typename std::map<ID::type,ParCell<C> >::iterator it=remoteCells.begin(); it!=remoteCells.end(); ++it) {
      it->second.dataptr = new C;
   }
}

/**
 * Prerequisite: the neighour lists must be up-to-date.
 */
template<class C>
void ParGrid<C>::buildExchangeLists() {
   // Delete old data
   receiveList.clear();
   sendList.clear();
   remoteToLocal.clear();
   
   bool hasRemotes;
   typename std::map<ID::type,ParCell<C> >::iterator it = localCells.begin();
   while (it != localCells.end()) {
      hasRemotes = false;
      it->second.N_remNbrs = 0;
      // Go through the cell's neighbour list
      for (std::map<uchar,ID::type>::const_iterator itt=it->second.neighbours.begin(); itt!=it->second.neighbours.end(); ++itt) {
	 // Check that the neighbour exists and that it is not assigned to this process:
	 const ID::type nbrID = itt->second;
	 if (nbrID == std::numeric_limits<ID::type>::max()) continue;
	 #ifndef NDEBUG
	 if (hostProcesses.find(nbrID) == hostProcesses.end()) {
	    std::cerr << "ParGrid CRITICAL ERROR: ";
	    std::cerr << "Proc #" << myrank << " nbrID " << nbrID << " not in hostProcesses!" << std::endl;
	    exit(1);
	 }
	 #endif
	 const int hostID = hostProcesses[nbrID];
	 if (hostID == myrank) continue;
	 // Cell has a remote neighbour. Add this cell to sendList, and its remote 
	 // neighbour to receiveList.
	 hasRemotes = true;
	 sendList[std::pair<ID::type,int>(it->first,hostID)] = 1;
	 receiveList[nbrID] = hostID;
	 
	 // Add association between a remote cell and this local cell 
	 // (this is for testSome):
	 ++(it->second.N_remNbrs);
	 remoteToLocal.insert(std::pair<ID::type,ID::type>(nbrID,it->first));
      }
      it->second.hasRemoteNeighbours = hasRemotes;
      ++it;
   }
}

template<class C> void ParGrid<C>::buildInitialGrid() {
   cuint N_total = unrefSize_x*unrefSize_y*unrefSize_z;
   
   // Determine how many cells should be built on this process:
   cuint N_cells = N_total / N_processes;
   uint N_extra = 0;
   if (N_total % N_processes > 0) {
      if (myrank < N_total % N_processes) ++N_extra;
   }
   
   // Determine the global ID of each cell to build. Note that the number of 
   // unrefined cells is probably not exactly divisible by the number of MPI
   // processes, and thus some processes should get an extra cell whose index
   // is calculated to i_extra.
   cuint i_start = myrank*(N_total / N_processes);
   cuint i_end   = i_start + N_cells;
   uint i_extra = 0;
   if (N_extra > 0) i_extra = N_processes*(N_total/N_processes) + (myrank);
   // Push all cells to localCells:
   for (uint i=i_start; i<i_end; ++i) {
      ParCell<C> p;
      calculateUnrefinedIndices(i,p.unrefInd_i,p.unrefInd_j,p.unrefInd_k);
      p.xcrd = grid_xmin + p.unrefInd_i*unref_dx;
      p.ycrd = grid_ymin + p.unrefInd_j*unref_dy;
      p.zcrd = grid_zmin + p.unrefInd_k*unref_dz;
      p.dx = unref_dx;
      p.dy = unref_dy;
      p.dz = unref_dz;
      localCells[i] = p;
   }
   // If an extra cell was assigned to this process, push it to localCells:
   if (i_extra > 0) {
      ParCell<C> p;
      calculateUnrefinedIndices(i_extra,p.unrefInd_i,p.unrefInd_j,p.unrefInd_k);
      p.xcrd = grid_xmin + p.unrefInd_i*unref_dx;
      p.ycrd = grid_ymin + p.unrefInd_j*unref_dy;
      p.zcrd = grid_zmin + p.unrefInd_k*unref_dz;
      p.dx = unref_dx;
      p.dy = unref_dy;
      p.dz = unref_dz;
      localCells[i_extra] = p;
   }
}

template<class C> void ParGrid<C>::buildUnrefNeighbourLists() {
   ID::type globalID;
   ID::type i,j,k;
   ID::type x_neg,x_pos,y_neg,y_pos,z_neg,z_pos;
   for (typename std::map<ID::type,ParCell<C> >::iterator it=localCells.begin(); it!=localCells.end(); ++it) {
      globalID = it->first;
      calculateUnrefinedIndices(globalID,i,j,k);
      // Determine x-indices of left and right x-neighbours:
      if (i == 0) {
	 if (periodic_x == true) x_neg = unrefSize_x-1;
	 else x_neg = std::numeric_limits<ID::type>::max();
      } else
	x_neg = i-1;
      if (i == unrefSize_x-1) {
	 if (periodic_x == true) x_pos = 0;
	 else x_pos = std::numeric_limits<ID::type>::max();
      } else
	x_pos = i+1;
      // Determine y-indices of left and right y-neighbours:
      if (j == 0) {
	 if (periodic_y == true) y_neg = unrefSize_y-1;
	 else y_neg = std::numeric_limits<ID::type>::max();
      } else
	y_neg = j-1;
      if (j == unrefSize_y-1) {
	 if (periodic_y == true) y_pos = 0;
	 else y_pos = std::numeric_limits<ID::type>::max();
      } else
	y_pos = j+1;
      // Determine z-indices of left and right z-neighbours:
      if (k == 0) {
	 if (periodic_z == true) z_neg = unrefSize_z-1;
	 else z_neg = std::numeric_limits<ID::type>::max();
      } else
	z_neg = k-1;
      if (k == unrefSize_z-1) {
	 if (periodic_z == true) z_pos = 0;
	 else z_pos = std::numeric_limits<ID::type>::max();
      } else
	z_pos = k+1;
      
      // Store calculated neighbour indices:
      if (x_neg != std::numeric_limits<ID::type>::max()) x_neg = calculateUnrefinedIndex(x_neg,j,k);
      if (x_pos != std::numeric_limits<ID::type>::max()) x_pos = calculateUnrefinedIndex(x_pos,j,k);
      if (y_neg != std::numeric_limits<ID::type>::max()) y_neg = calculateUnrefinedIndex(i,y_neg,k);
      if (y_pos != std::numeric_limits<ID::type>::max()) y_pos = calculateUnrefinedIndex(i,y_pos,k);
      if (z_neg != std::numeric_limits<ID::type>::max()) z_neg = calculateUnrefinedIndex(i,j,z_neg);
      if (z_pos != std::numeric_limits<ID::type>::max()) z_pos = calculateUnrefinedIndex(i,j,z_pos);
      
      if (x_neg != std::numeric_limits<ID::type>::max()) it->second.neighbours[ 0] = x_neg;
      if (x_pos != std::numeric_limits<ID::type>::max()) it->second.neighbours[ 8] = x_pos;
      if (y_neg != std::numeric_limits<ID::type>::max()) it->second.neighbours[16] = y_neg;
      if (y_pos != std::numeric_limits<ID::type>::max()) it->second.neighbours[24] = y_pos;
      if (z_neg != std::numeric_limits<ID::type>::max()) it->second.neighbours[32] = z_neg;
      if (z_pos != std::numeric_limits<ID::type>::max()) it->second.neighbours[40] = z_pos;
   }
}

/** Calculate the statistical weight of the given cell. This weight is passed to Zoltan to be 
 * used in load balancing. Here each cell is given the same weight, as given by parameter 
 * ParGrid::weightCell.
 * @param globalID The global ID of the spatial cell whose weight is requested.
 * @return The weight of the given spatial cell.
 */
template<class C> float ParGrid<C>::calculateCellWeight(const ID::type& globalID) {return weightCell;}

/** Calculate the statistical weight of the given edge. This weight is passed to Zoltan to 
 * be used in graph-based load balancing. Here each edge is given the same weight, as given by 
 * parameter ParGrid::weightEdge. The sum of all edge weights for a given spatial cell is then 
 * the number of existing neighbours multiplied by weightEdge, in accordance with 
 * ParGrid::calculateHyperedgeWeight.
 * @param globalID The global ID of the spatial cell whose edge weight is requested.
 * @return The weight of the requested edge.
 */
template<class C> float ParGrid<C>::calculateEdgeWeight(const ID::type& globalID) {return weightEdge;}

/** Calculate the statistical weight of the given hyperedge. This weight is passed to Zoltan 
 * to be used in hypergraph-based load balancing. Here the weight of a hyperedge is calculated 
 * as the number of existing neighbours multiplied by parameter ParGrid::weightEdge.
 * @param globalID The global ID of the hyperedge whose weight is requested.
 * @return The weight of the requested hyperedge.
 */
template<class C> float ParGrid<C>::calculateHyperedgeWeight(const ID::type& globalID) {
   uint N_neighbours = 0;
   N_neighbours = localCells[globalID].neighbours.size();
   return weightEdge*N_neighbours;
}

template<class C> uint ParGrid<C>::calculateNbrIndex(const ID::type& cellID,const ID::type& nbrID) {
   // Get pointers to cell and its neighbour:
   ParCell<C>* cellPtr = &(localCells[cellID]);
   ParCell<C>* nbrPtr;
   if (localCells.find(nbrID) != localCells.end()) nbrPtr = &(localCells[nbrID]);
   else if (remoteCells.find(nbrID) != remoteCells.end()) nbrPtr = &(remoteCells[nbrID]);
   else {
      std::cerr << "ParGrid ERROR: Could not find neighbour!" << std::endl << std::flush;
      exit(1);
   }
   
   if (cellPtr->refLevel == nbrPtr->refLevel) {
      Real dx = unref_dx / (cellPtr->refLevel + 1);
      Real dy = unref_dy / (cellPtr->refLevel + 1);
      Real dz = unref_dz / (cellPtr->refLevel + 1);
      if (nbrPtr->xcrd < cellPtr->xcrd - 0.1*dx) return  0;
      if (nbrPtr->xcrd > cellPtr->xcrd + 0.1*dx) return  4;
      if (nbrPtr->ycrd < cellPtr->ycrd - 0.1*dy) return  8;
      if (nbrPtr->ycrd > cellPtr->ycrd + 0.1*dy) return 12;
      if (nbrPtr->zcrd < cellPtr->zcrd - 0.1*dz) return 16;
      if (nbrPtr->zcrd > cellPtr->zcrd + 0.1*dz) return 20;
   }
   return 254;
}

template<class C> ID::type ParGrid<C>::calculateUnrefinedIndex(const ID::type& i,const ID::type& j,const ID::type& k) const {
   return k*unrefSize_y*unrefSize_x + j*unrefSize_x + i;
}

template<class C> void ParGrid<C>::calculateUnrefinedIndices(const ID::type& index,ID::type& i,ID::type& j,ID::type& k) const {
   cuint N_total = unrefSize_x*unrefSize_y*unrefSize_z;
   uint ind = index;
   k = ind / (unrefSize_x*unrefSize_y);
   ind -= k*unrefSize_x*unrefSize_y;
   j = ind / unrefSize_x;
   i = ind - j*unrefSize_x;
}

template<class C> bool ParGrid<C>::startNeighbourExchange(cuint& identifier) {
   bool rvalue = true;
   // Get the MPI_Datatype used in communication from user and pass it to MPI:
   // (This does not work with variable-length user data)
   if (MPItypeFreed == false) MPI_Type_free(&MPIdataType);
   C::getMPIdatatype(identifier,MPIdataType);
   int rcode = MPI_Type_commit(&MPIdataType);
   MPItypeFreed = false;
   #ifndef NDEBUG
      if (rcode != MPI_SUCCESS) {std::cerr << "ParGrid::startNeighbourExchange MPI_Type_commit failed!" << std::endl; rvalue=false;}
   #endif

   /*
   for (int i=0; i<N_processes; ++i) {
      if (i == myrank) {
	 std::cerr << "PROC #" << myrank << " RECEIVES:" << std::endl;
	 for (std::map<ID::type,int>::const_iterator it=receiveList.begin(); it!=receiveList.end(); ++it) {
	    std::cerr << "\t recv cell#" << it->first << " from proc#" << it->second << std::endl;
	 }
	 std::cerr << "PROC #" << myrank << " SENDS:" << std::endl;
	 for (std::map<std::pair<ID::type,int>,char>::const_iterator it=sendList.begin(); it!=sendList.end(); ++it) {
	    std::cerr << "\t send cell#" << it->first.first << " to proc #" << it->first.second << std::endl;
	 }
      }
      barrier();
   }
   barrier();
   */
   // Post receives for each remote cell:
   for (std::map<ID::type,int>::const_iterator it=receiveList.begin(); it!=receiveList.end(); ++it) {
      void* const buffer = (remoteCells[it->first].dataptr)->getBaseAddress(identifier); // Starting address of receive buffer
      const int count = 1;                            // Number of elements in receive buffer
      const int source = it->second;                  // Rank of source MPI process
      const int tag = it->first;                      // Message tag (use global ID)
      MPIrecvRequests.push_back(MPI_Request());
      if (MPI_Irecv(buffer,count,MPIdataType,source,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back())) != MPI_SUCCESS) rvalue=false;
   }

   // Post sends for each boundary cell:
   for (std::map<std::pair<ID::type,int>,char>::const_iterator it=sendList.begin(); it!=sendList.end(); ++it) {
      void* const buffer = (localCells[it->first.first].dataptr)->getBaseAddress(identifier); // Starting address of send buffer
      const int count = 1;                            // Number of elements in send buffer
      const int dest = it->first.second;              // Rank of destination MPI process
      const int tag = it->first.first;                // Message tag
      MPIsendRequests.push_back(MPI_Request());
      if (MPI_Issend(buffer,count,MPIdataType,dest,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) rvalue=false;
   }

   // Clear the number of received neighbours:
   N_receivesRemaining = receiveList.size();
   for (typename std::map<ID::type,ParCell<C> >::iterator it=localCells.begin(); it!=localCells.end(); ++it) 
     it->second.N_receivedRemNbrs = 0;
   #ifndef NDEBUG
      if (readyCells.size() != 0) {
	 std::cerr << "ParGrid::startNeighbourExchange readyCells is not empty!" << std::endl << std::flush;
	 exit(1);
      }
   #endif
   readyCells.clear();
   
   return rvalue;
}

/** Return the global ID of the unrefined cell in which the given point belongs.
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @param z Z-coordinate.
 * @return The global ID of the cell containing the given point. If the cell does not exists, 
 * value numeric_limits<ID::type>::max() is returned.
 */
template<class C> ID::type ParGrid<C>::get_cell(const Real& x,const Real& y,const Real& z) const {
   // Check that the cell exists:
   if (x < grid_xmin || x > grid_xmax) return std::numeric_limits<ID::type>::max();
   if (y < grid_ymin || y > grid_ymax) return std::numeric_limits<ID::type>::max();
   if (z < grid_zmin || z > grid_zmax) return std::numeric_limits<ID::type>::max();
   // Calculate unrefined indices:
   const ID::type I = static_cast<ID::type>(floor((x - grid_xmin) / unref_dx));
   const ID::type J = static_cast<ID::type>(floor((y - grid_ymin) / unref_dy));
   const ID::type K = static_cast<ID::type>(floor((z - grid_zmin) / unref_dz));
   return calculateUnrefinedIndex(I,J,K);
}

template<class C> Real ParGrid<C>::get_cell_x(const ID::type& globalID) const {
   // Try to find the cell from maps:
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.xcrd + 0.5*unref_dx;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.xcrd + 0.5*unref_dx;
      
   // Calculate the x-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_xmin + (i+0.5)*unref_dx;
}

template<class C> Real ParGrid<C>::get_cell_x_size(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return std::numeric_limits<Real>::max();
   return it->second.dx;
}

template<class C> Real ParGrid<C>::get_cell_y(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.ycrd + 0.5*unref_dy;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.ycrd + 0.5*unref_dy;
   
   // Calculate the y-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_ymin + (j+0.5)*unref_dy;
}

template<class C> Real ParGrid<C>::get_cell_y_size(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return std::numeric_limits<Real>::max();
   return it->second.dy;
}

template<class C> Real ParGrid<C>::get_cell_z(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.zcrd + 0.5*unref_dz;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.zcrd + 0.5*unref_dz;
   
   // Calculate the z-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_zmin + (k+0.5)*unref_dz;
}

template<class C> Real ParGrid<C>::get_cell_z_size(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return std::numeric_limits<Real>::max();
      return it->second.dz;
}

template<class C> Real ParGrid<C>::get_cell_x_max(const ID::type& globalID) const {
   // Try to find the cell from maps:
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.xcrd + unref_dx;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.xcrd + unref_dx;
   
   // Calculate the x-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_xmin + (i+1)*unref_dx;
}

template<class C> Real ParGrid<C>::get_cell_x_min(const ID::type& globalID) const {
   // Try to find the cell from maps:
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.xcrd;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.xcrd;

   // Calculate the x-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_xmin;
}

template<class C> Real ParGrid<C>::get_cell_y_max(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.ycrd + unref_dy;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.ycrd + unref_dy;
   
   // Calculate the y-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_ymin + (j+1)*unref_dy;
}

template<class C> Real ParGrid<C>::get_cell_y_min(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.ycrd;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.ycrd;
   
   // Calculate the y-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_ymin;
}

template<class C> Real ParGrid<C>::get_cell_z_max(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.zcrd + unref_dz;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.zcrd + unref_dz;
   
   // Calculate the z-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_zmin + (k+1)*unref_dz;
}

template<class C> Real ParGrid<C>::get_cell_z_min(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it != localCells.end()) return it->second.zcrd;
   it = remoteCells.find(globalID);
   if (it != remoteCells.end()) return it->second.zcrd;
   
   // Calculate the z-coordinate of the unrefined cell:
   ID::type i,j,k;
   calculateUnrefinedIndices(globalID,i,j,k);
   return grid_zmin;
}

template<class C> template<class CONT> void ParGrid<C>::getAllCells(CONT& rlist) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.begin();
   while (it != localCells.end()) {
      rlist.push_back(it->first);
      ++it;
   }
   
   it = remoteCells.begin();
   while (it != remoteCells.end()) {
      rlist.push_back(it->first);
      ++it;
   }
}

template<class C> template<class CONT> void ParGrid<C>::getBoundaryCells(CONT& rlist) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.begin();
   while (it != localCells.end()) {
      if (it->second.hasRemoteNeighbours == true) rlist.push_back(it->first);
      ++it;
   }
}

template<class C> template<class CONT> void ParGrid<C>::getCells(CONT& rlist) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.begin();
   while (it != localCells.end()) {
      rlist.push_back(it->first);
      ++it;
   }
}

template<class C> template<class CONT> void ParGrid<C>::getExistingNeighbours(CONT& rlist,const ID::type& globalID) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return;

   for (std::map<uchar,ID::type>::const_iterator itt = it->second.neighbours.begin(); itt!=it->second.neighbours.end(); ++itt) {
      rlist.push_back(itt->second);
   }
}

template<class C> bool ParGrid<C>::getHost(const ID::type& cellID,int& host) const {
   std::map<ID::type,int>::const_iterator it = hostProcesses.find(cellID);
   if (it == hostProcesses.end()) return false;
   host = it->second;
   return true;
}

template<class C> template<class CONT> void ParGrid<C>::getInnerCells(CONT& rlist) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.begin();
   while (it != localCells.end()) {
      if (it->second.hasRemoteNeighbours == false) rlist.push_back(it->first);
      ++it;
   }
}

template<class C>
ID::type ParGrid<C>::getNeighbour(const ID::type globalID,cuchar& nbrTypeID) const {
   // Check if this process has the cell with given global ID:
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return std::numeric_limits<ID::type>::max();

   std::map<uchar,ID::type>::const_iterator itt = it->second.neighbours.find(nbrTypeID);
   if (itt == it->second.neighbours.end()) return std::numeric_limits<ID::type>::max();
   return itt->second;
}

template<class C> uint ParGrid<C>::getNumberOfRemoteNeighbours(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return std::numeric_limits<ID::type>::max();   
   return it->second.N_remNbrs;
}

template<class C> bool ParGrid<C>::getReadyCell(ID::type& globalID) {
   // If there are no ready cells, return immediately.
   if (readyCells.size() == 0) return false;
   
   // Pop and return the first cell in readyCells:
   globalID = readyCells.front();
   readyCells.pop_front();
   return true;
}

template<class C> bool ParGrid<C>::getReadyCell2(ID::type& globalID) {
   if (readyCells2.size() == 0) return false;
   
   globalID = readyCells2.front();
   readyCells2.pop_front();
   return true;
}

template<class C> template<class CONT> void ParGrid<C>::getReceiveList(CONT& rlist) const {
   rlist.clear();
   for (std::map<ID::type,int>::const_iterator it = receiveList.begin(); it != receiveList.end(); ++it) {
      rlist.push_back(std::pair<ID::type,int>(it->first,it->second));
   }
}

template<class C> unsigned char ParGrid<C>::getRefinementLevel(const ID::type& globalID) const {
   if (localCells.find(globalID) == localCells.end()) return std::numeric_limits<unsigned char>::max();
   return localCells[globalID].refLevel;
}

template<class C> uint ParGrid<C>::getRemainingReceives() const {
   #ifndef NDEBUG
      if (N_receivesRemaining > 0.1*std::numeric_limits<uint>::max()) {
	 std::cerr << "N_receivesRemaining has a high value " << N_receivesRemaining << std::endl;
      }
   #endif
   //if (myrank == 0) {std::cerr << "ParGrid::getRemainingReceives N_receivesRemaining = " << N_receivesRemaining << std::endl;} // TEST
   return N_receivesRemaining;
}

template<class C> uint ParGrid<C>::getRemainingReceives2() const {
   return N_receivesRemaining2;
}

template<class C> template<class CONT> void ParGrid<C>::getRemoteCells(CONT& rlist) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = remoteCells.begin();
   while (it != remoteCells.end()) {
      rlist.push_back(it->first);
      ++it;
   }
}

template<class C> ID::type ParGrid<C>::getRemoteNeighbour(const ID::type& globalID,const uchar& nbrTypeID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return std::numeric_limits<ID::type>::max();
   
   std::map<uchar,ID::type>::const_iterator neighbour = it->second.neighbours.find(nbrTypeID);
   if (neighbour == it->second.neighbours.end()) return std::numeric_limits<ID::type>::max();
   
   std::map<ID::type,int>::const_iterator host = hostProcesses.find(neighbour->second);
   if (host->second == myrank) return std::numeric_limits<ID::type>::max();
   return neighbour->second;
}

template<class C> template<class CONT> void ParGrid<C>::getRemoteNeighbours(const ID::type& globalID,CONT& rlist) const {
   rlist.clear();
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return;
   
   for (std::map<uchar,ID::type>::const_iterator nbrs = it->second.neighbours.begin(); nbrs != it->second.neighbours.end(); ++nbrs) {
      if (remoteCells.find(nbrs->second) != remoteCells.end()) rlist.push_back(nbrs->second);
   }
}

template<class C> template<class CONT> void ParGrid<C>::getRemoteToLocalMapping(const ID::type& remoteID,CONT& rlist) const {
   rlist.clear();
   for (std::multimap<ID::type,ID::type>::const_iterator it = remoteToLocal.lower_bound(remoteID); it != remoteToLocal.upper_bound(remoteID); ++it) 
     rlist.push_back(it->second);
}

template<class C> template<class CONT> void ParGrid<C>::getSendList(CONT& rlist) const {
   rlist.clear();
   for (std::map<std::pair<ID::type,int>,char>::const_iterator it = sendList.begin(); it != sendList.end(); ++it) {
      rlist.push_back(std::pair<ID::type,int>(it->first.first,it->first.second));
   }
}

template<class C> bool ParGrid<C>::hasRemoteNeighbours(const ID::type& globalID) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID);
   if (it == localCells.end()) return false;
   return it->second.hasRemoteNeighbours;
}

template<class C> bool ParGrid<C>::initialize() {
   initialized = true;
   MPItypeFreed = true;
   // At this point we do not know who owns the boundary cells, 
   // i.e. with which processes I need to exchange data with.
   if (syncCellAssignments() == false) initialized = false;

   buildExchangeLists();
   initialLoadBalance();
   buildExchangeLists();
   // Allocate memory for user data:
   std::cerr << "ParGrid: Allocating " << localCells.size() << " local cells" << std::endl;
   for (typename std::map<ID::type,ParCell<C> >::iterator it=localCells.begin(); it!= localCells.end(); ++it) {
      it->second.dataptr = new C;
   }
   for (std::map<ID::type,int>::const_iterator it=receiveList.begin(); it!=receiveList.end(); ++it) {
      remoteCells[it->first];
   }
   std::cerr << "ParGrid: Allocating " << remoteCells.size() << " remote cells" << std::endl;
   for (typename std::map<ID::type,ParCell<C> >::iterator it=remoteCells.begin(); it!=remoteCells.end(); ++it) {
      it->second.dataptr = new C;
   }
   syncCellCoordinates();
   return initialized;
}

template<class C>
bool ParGrid<C>::initialize(cuint& xsize,cuint& ysize,cuint& zsize,creal& xmin,creal& ymin,creal& zmin,
			    creal& xmax,creal& ymax,creal& zmax) {
   MPItypeFreed = true;
   unrefSize_x = xsize;
   unrefSize_y = ysize;
   unrefSize_z = zsize;
   unref_dx = (xmax - xmin) / xsize;
   unref_dy = (ymax - ymin) / ysize;
   unref_dz = (zmax - zmin) / zsize;
   grid_xmin = xmin;
   grid_xmax = xmax;
   grid_ymin = ymin;
   grid_ymax = ymax;
   grid_zmin = zmin;
   grid_zmax = zmax;
   
   periodic_x = false;
   periodic_y = false;
   periodic_z = false;
   
   buildInitialGrid();    // Build an initial guess for the grid
   buildUnrefNeighbourLists();
   syncCellAssignments();
  
   // Load balance grid. Invalidates contents of hostProcesses
   // and neighbour lists (since no data is actually transmitted).
   initialLoadBalance();
   syncCellAssignments();
   buildExchangeLists();  // Build send/receive lists.

   // Now cells have been assigned to this process, and we also know which 
   // cells this process will receive from other MPI processes. Allocate
   // memory:
   allocateCells();
   syncCellCoordinates();
   initialized = true;
   return initialized;
}

template<class C> bool ParGrid<C>::initialLoadBalance() {
   bool rvalue = true;

   // Request load balance from Zoltan, and get cells which should be imported and exported:
   int changes,N_globalIDs,N_localIDs,N_import,N_export;
   int* importProcesses;
   int* importParts;
   int* exportProcesses;
   int* exportParts;
   ZOLTAN_ID_PTR importGlobalIDs;
   ZOLTAN_ID_PTR importLocalIDs;
   ZOLTAN_ID_PTR exportGlobalIDs;
   ZOLTAN_ID_PTR exportLocalIDs;
   if (zoltan->LB_Partition(changes,N_globalIDs,N_localIDs,N_import,importGlobalIDs,importLocalIDs,
			    importProcesses,importParts,N_export,exportGlobalIDs,exportLocalIDs,
			    exportProcesses,exportParts) != ZOLTAN_OK) { 
      std::cerr << "ParGrid FATAL ERROR: Zoltan failed on load balancing!" << std::endl << std::flush;
      zoltan->LB_Free_Part(&importGlobalIDs,&importLocalIDs,&importProcesses,&importParts);
      zoltan->LB_Free_Part(&exportGlobalIDs,&exportLocalIDs,&exportProcesses,&exportParts);
      rvalue = false;
      exit(1);
   }
   
   // Do a pass of hostProcesses updates. We need to exchange lists of migrating cells 
   // and their new hosts to every neighbouring process. A process is a neighbouring process if that 
   // process has at least one of my local cells' (remote) neighbours.
   std::set<int>                       neighbouringHosts;      // Ranks of neighbouring processes
   for (std::map<ID::type,int>::const_iterator it=hostProcesses.begin(); it!=hostProcesses.end(); ++it) {
      if (it->second != myrank) neighbouringHosts.insert(it->second);
   }

   // Update hostProcesses based on local information available to this process. 
   // This needs to occur after neighbouringHosts has been updated.
   for (int i=0; i<N_export; ++i) hostProcesses[exportGlobalIDs[i]] = exportProcesses[i];
   for (int i=0; i<N_import; ++i) hostProcesses[importGlobalIDs[i]] = myrank;

   // Exchange information on cell migrations between neighbouring processes.
   int* neighbourChanges                    = new int[neighbouringHosts.size()];              // Number of exports nbr. proc. has
   ZOLTAN_ID_TYPE** neighbourMigratingCells = new ZOLTAN_ID_TYPE* [neighbouringHosts.size()]; // Global IDs of cells nbr. proc. is exporting
   int** neighbourMigratingHosts            = new int* [neighbouringHosts.size()];            // New hosts for cells nbr. proc. exports
   ID::type counter = 0;
   MPIrecvRequests.resize(neighbouringHosts.size());
   for (std::set<int>::const_iterator it=neighbouringHosts.begin(); it!=neighbouringHosts.end(); ++it) {
      // Receive nmbr. of exports from nbr. proc. *it
      if (MPI_Irecv(&(neighbourChanges[counter]),1,MPI_INT,*it,*it   ,MPI_COMM_WORLD,&(MPIrecvRequests[counter])) != MPI_SUCCESS) rvalue = false;
      // Send nmbr. of my exports, list of exported cells & new hosts to nbr. proc *it
      MPI_Request mpiRequest1;
      MPI_Request mpiRequest2;
      MPI_Request mpiRequest3;
      if (MPI_Isend(&N_export      ,1       ,MPI_INT                   ,*it,myrank,MPI_COMM_WORLD,&mpiRequest1) != MPI_SUCCESS) rvalue = false;
      if (MPI_Isend(exportGlobalIDs,N_export,MPI_Type<ZOLTAN_ID_TYPE>(),*it,myrank,MPI_COMM_WORLD,&mpiRequest2) != MPI_SUCCESS) rvalue = false;
      if (MPI_Isend(exportProcesses,N_export,MPI_INT                   ,*it,myrank,MPI_COMM_WORLD,&mpiRequest3) != MPI_SUCCESS) rvalue = false;
      MPI_Request_free(&mpiRequest1);
      MPI_Request_free(&mpiRequest2);
      MPI_Request_free(&mpiRequest3);
      ++counter;
   }
   // Wait for the numbers of exported cells to arrive from each rem. nbr.
   if (MPI_Waitall(neighbouringHosts.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE) != MPI_SUCCESS) rvalue = false;

   // Allocate buffers for receiving lists of exported cells & new hosts based on
   // numbers of exported cells, and receive data:
   counter = 0;
   MPIrecvRequests.resize(neighbouringHosts.size()*2);
   for (std::set<int>::const_iterator it=neighbouringHosts.begin(); it!=neighbouringHosts.end(); ++it) {
      neighbourMigratingCells[counter] = new ZOLTAN_ID_TYPE[neighbourChanges[counter]];
      neighbourMigratingHosts[counter] = new int[neighbourChanges[counter]];
      int tag = *it;
      if (MPI_Irecv(neighbourMigratingCells[counter],neighbourChanges[counter],MPI_Type<ZOLTAN_ID_TYPE>(),*it,tag,MPI_COMM_WORLD,&(MPIrecvRequests[2*counter+0])) != MPI_SUCCESS) rvalue = false;
      if (MPI_Irecv(neighbourMigratingHosts[counter],neighbourChanges[counter],MPI_INT                   ,*it,tag,MPI_COMM_WORLD,&(MPIrecvRequests[2*counter+1])) != MPI_SUCCESS) rvalue = false;
      ++counter;
   }
   if (MPI_Waitall(2*neighbouringHosts.size(),&(MPIrecvRequests[0]),MPI_STATUSES_IGNORE) != MPI_SUCCESS) rvalue = false;

   // Go through the lists of migrations received from neighbouring hosts. If those
   // lists contain cells that I have in hostProcesses, update their statuses:
   counter = 0;
   for (std::set<int>::const_iterator it=neighbouringHosts.begin(); it!=neighbouringHosts.end(); ++it) {
      for (int i=0; i<neighbourChanges[counter]; ++i) {
	 std::map<ID::type,int>::iterator tmp = hostProcesses.find(neighbourMigratingCells[counter][i]);
	 if (tmp != hostProcesses.end()) tmp->second = neighbourMigratingHosts[counter][i];
      }
      ++counter;
   }
   delete neighbourChanges;
   delete [] neighbourMigratingCells;
   delete [] neighbourMigratingHosts;
   neighbourChanges = NULL;
   neighbourMigratingCells = NULL;
   neighbourMigratingHosts = NULL;

   // **********************************************************
   // ***** SEND EXPORTED CELLS AND RECEIVE IMPORTED CELLS *****
   // **********************************************************
   
   // Create a suitable MPI type for sending & receiving cells:
   struct CellData {
      ID::type IDs[49];
      Real coords[6];
      uchar nbrTypes[49];
      int nbrHosts[49];
   } cellData;
   
   MPI_Datatype types[] = {MPI_Type<ID::type>(),MPI_Type<Real>(),MPI_Type<uchar>(),MPI_Type<int>()};
   int lengths[] = {49,6,49,49};
   MPI_Aint disp[4];
   disp[0] = 0;
   disp[1] = reinterpret_cast<char*>(cellData.coords)-reinterpret_cast<char*>(cellData.IDs);
   disp[2] = reinterpret_cast<char*>(cellData.nbrTypes)-reinterpret_cast<char*>(cellData.IDs);
   disp[3] = reinterpret_cast<char*>(cellData.nbrHosts)-reinterpret_cast<char*>(cellData.IDs);
   MPI_Type_create_struct(4,lengths,disp,types,&MPIdataType);
   MPI_Type_commit(&MPIdataType);
   MPItypeFreed = false;
   
   // Insert each imported cell to receiveBuffer and post receives:
   CellData* sendBuffer    = new CellData[N_export];
   CellData* receiveBuffer = new CellData[N_import];
   MPIrecvRequests.resize(N_import);
   for (int i=0; i<N_import; ++i) {
      void* const buffer = &(receiveBuffer[i]);
      const int  count   = 1;
      const int  tag     = 0;
      const int  source  = importProcesses[i];
      if (MPI_Irecv(buffer,count,MPIdataType,source,tag,MPI_COMM_WORLD,&(MPIrecvRequests[i])) != MPI_SUCCESS) rvalue=false;
   }      

   // Insert each exported cell to sendBuffer and post send:
   MPIsendRequests.resize(N_export);
   for (int i=0; i<N_export; ++i) {
      sendBuffer[i].IDs[0] = exportGlobalIDs[i];
      sendBuffer[i].nbrTypes[0] = localCells[exportGlobalIDs[i]].refLevel;
      int index=1;
      for (std::map<uchar,ID::type>::const_iterator it=localCells[exportGlobalIDs[i]].neighbours.begin(); it!=localCells[exportGlobalIDs[i]].neighbours.end(); ++it) {
	 sendBuffer[i].IDs[index] = it->second;
	 sendBuffer[i].nbrTypes[index] = it->first;
	 sendBuffer[i].nbrHosts[index] = hostProcesses[it->second];
	 // Reduce reference count to neighbours:
	 --nbrReferences[it->second];
	 
	 ++index;
      }
      for (int j=index; j<49; ++j) {sendBuffer[i].IDs[j] = std::numeric_limits<ID::type>::max();}
      sendBuffer[i].coords[0] = localCells[exportGlobalIDs[i]].xcrd;
      sendBuffer[i].coords[1] = localCells[exportGlobalIDs[i]].ycrd;
      sendBuffer[i].coords[2] = localCells[exportGlobalIDs[i]].zcrd;
      sendBuffer[i].coords[3] = localCells[exportGlobalIDs[i]].dx;
      sendBuffer[i].coords[4] = localCells[exportGlobalIDs[i]].dy;
      sendBuffer[i].coords[5] = localCells[exportGlobalIDs[i]].dz;
      
      void* const buffer = &(sendBuffer[i]);
      const int   count  = 1;
      const int   tag    = 0;
      const int   dest   = exportProcesses[i];
      if (MPI_Issend(buffer,count,MPIdataType,dest,tag,MPI_COMM_WORLD,&(MPIsendRequests[i])) != MPI_SUCCESS) rvalue=false;
   }
   // Wait for sends & receives to complete. waitAll() unallocates MPIdataType.
   waitAll();

   // Insert all received cells:
   for (int i=0; i<N_import; ++i) {
      ParCell<C> p;
      const ID::type globalID = receiveBuffer[i].IDs[0];
      p.refLevel = receiveBuffer[i].nbrTypes[0];
      p.xcrd     = receiveBuffer[i].coords[0];
      p.ycrd     = receiveBuffer[i].coords[1];
      p.zcrd     = receiveBuffer[i].coords[2];
      p.dx       = receiveBuffer[i].coords[3];
      p.dy       = receiveBuffer[i].coords[4];
      p.dz       = receiveBuffer[i].coords[5];
      int j=1;
      while (j < 49 && receiveBuffer[i].IDs[j] != std::numeric_limits<ID::type>::max()) {
	 p.neighbours[receiveBuffer[i].nbrTypes[j]] = receiveBuffer[i].IDs[j];
	 hostProcesses[receiveBuffer[i].IDs[j]] = receiveBuffer[i].nbrHosts[j];
	 // Increase reference count for all neighbours:
	 if (nbrReferences.find(receiveBuffer[i].IDs[j]) == nbrReferences.end()) nbrReferences[receiveBuffer[i].IDs[j]] = 0;
	 ++nbrReferences[receiveBuffer[i].IDs[j]];
	 
	 ++j;
      }
      localCells[globalID] = p;
   }
   delete sendBuffer;
   delete receiveBuffer;
   sendBuffer = NULL;
   receiveBuffer = NULL;
   
   // Erase exported cells:
   for (int i=0; i<N_export; ++i) localCells.erase(exportGlobalIDs[i]);

   // Erase entries from hostProcesses with 0 references remaining:
   std::map<ID::type,uint>::iterator it = nbrReferences.begin();
   while (it != nbrReferences.end()) {
      if (it->second == 0) {
	 std::map<ID::type,uint>::iterator tmp = it;
	 ++it;
	 nbrReferences.erase(tmp);
	 hostProcesses.erase(tmp->first);
      } else {
	 ++it;
      }
   }
   
   // Deallocate Zoltan arrays:
   zoltan->LB_Free_Part(&importGlobalIDs,&importLocalIDs,&importProcesses,&importParts);
   zoltan->LB_Free_Part(&exportGlobalIDs,&exportLocalIDs,&exportProcesses,&exportParts);
   return rvalue;
}
/*
template<class C> bool ParGrid<C>::loadBalance() {
   bool rvalue = true;
   int changes,N_globalIDs,N_localIDs,N_import,N_export;
   int* importProcesses;
   int* importParts;
   int* exportProcesses;
   int* exportParts;
   ZOLTAN_ID_PTR importGlobalIDs;
   ZOLTAN_ID_PTR importLocalIDs;
   ZOLTAN_ID_PTR exportGlobalIDs;
   ZOLTAN_ID_PTR exportLocalIDs;
   
   // Ask Zoltan the cells which should be imported and exported:
   if (zoltan->LB_Partition(changes,N_globalIDs,N_localIDs,N_import,importGlobalIDs,importLocalIDs,
			    importProcesses,importParts,N_export,exportGlobalIDs,exportLocalIDs,
			    exportProcesses,exportParts) == ZOLTAN_OK) {
      
      for (int i=0; i<N_processes; ++i) {
	 if (i == myrank) {
	    std::cerr << "Proc #" << myrank << " exporting cells:" << std::endl;
	    for (int j=0; j<N_export; ++j) {
	       std::cerr << exportGlobalIDs[j] << " to proc #" << exportProcesses[j] << " to part " << exportParts[j] << std::endl;
	    }
	    std::cerr << std::endl;
	    std::cerr << "Proc #" << myrank << " importing cells:" << std::endl;
	    for (int j=0; j<N_import; ++j) {
	       std::cerr << importGlobalIDs[j] << " to proc #" << importProcesses[j] << " to part " << importParts[j] << std::endl;
	    }
	    std::cerr << std::endl;
	 }
	 barrier();
      }
   } else {
      std::cerr << "(ParGrid) ERROR: Zoltan failed to balance load!" << std::endl;
   }
}
*/
template<class C> std::string ParGrid<C>::loadBalanceMethod(const LBM& method) {
   switch (method) {
    case Block:
      return "BLOCK";
      break;
    case Random:
      return "RANDOM";
      break;
    case RCB:
      return "RCB";
      break;
    case RIB:
      return "RIB";
      break;
    case HSFC:
      return "HSFC";
      break;
    case Graph:
      return "GRAPH";
      break;
    case Hypergraph:
      return "HYPERGRAPH";
      break;
    case Hierarchical:
      return "HIER";
      break;
   }
   return std::string("");
}

/** Return a pointer to the cell with the given global index.
 * @param id The global ID of the requested cell.
 * @return If NULL, the cell was not found on this process. Otherwise a pointer to 
 * the cell data is returned.
 */
template<class C> C* ParGrid<C>::operator[](const ID::type& id) const {
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(id);
   if (it != localCells.end()) return it->second.dataptr;
   it = remoteCells.find(id);
   if (it != remoteCells.end()) return it->second.dataptr;
   return NULL;
}

template<class C> void ParGrid<C>::print() const {
   std::cerr << "ParGrid:" << std::endl;
   std::cerr << "\tMPI reports " << N_processes << " processes." << std::endl;
   std::cerr << "\tMPI says my rank is " << myrank << std::endl;
}

template<class C> void ParGrid<C>::printTransfers() const {
   mpilogger << "(PARGRID): Sends and receives per MPI rank:" << std::endl;
   for (int i=0; i<N_processes; ++i) {
      if (i == myrank) continue;
      uint N_sends = 0;
      uint N_recvs = 0;
      for (std::map<std::pair<ID::type,int>,char>::const_iterator it=sendList.begin(); it!=sendList.end(); ++it) {
	 if (it->first.second == i) ++N_sends;
      }
      for (std::map<ID::type,int>::const_iterator it=receiveList.begin(); it!=receiveList.end(); ++it) {
	 if (it->second == i) ++N_recvs;
      }
      if (N_sends == 0 && N_recvs == 0) continue;
      mpilogger << "\tProc " << "\t#" << i << "\tSends: " << N_sends << "\tReceives: " << N_recvs << std::endl;
   }
   mpilogger << "\t Local Cells = " << localCells.size() << "\t Remote Cells = " << remoteCells.size() << std::endl << write;
}

template<class C> int ParGrid<C>::processes() const {return N_processes;}

template<class C> int ParGrid<C>::rank() const {return myrank;}

template<class C> bool ParGrid<C>::returnReadyCell(const ID::type& globalID) {
   #ifndef NDEBUG
   if (localCells.find(globalID) == localCells.end()) return false;
   #endif
   std::cerr << "ParGrid: " << globalID << " returned to readyCells" << std::endl;
   readyCells.push_back(globalID);
   return true;
}

template<class C> bool ParGrid<C>::returnReadyCell2(const ID::type& globalID) {
   readyCells2.push_back(globalID);
   return true;
}

template<class C> uint ParGrid<C>::singleModeTestSome() {
   if (N_receivesRemaining == 0) return 0;
   
   // Ensure that enough status messages have been allocated:
   if (MPIstatuses.size() < MPIrecvRequests.size()) MPIstatuses.resize(MPIrecvRequests.size());
   
   // Test if some of the active receives has completed, and insert all completed 
   // localIDs to readyCells:
   int receives = 0;
   int readyIndices[MPIrecvRequests.size()];
   if (MPI_Testsome(MPIrecvRequests.size(),&(MPIrecvRequests[0]),&receives,readyIndices,&(MPIstatuses[0])) != MPI_SUCCESS) {
      std::cerr << "ParGrid::singleModeTestSome failed on process #" << myrank << std::endl;
      exit(1);
   }
   if (receives == MPI_UNDEFINED) {
      std::cerr << "ParGrid::singleModeTestSome receives is MPI_UNDEFINED!" << std::endl;
      exit(1);
   }
   for (int i=0; i<receives; ++i) readyCells.push_back(localReceiveIDs[readyIndices[i]]);
   N_receivesRemaining -= receives;
   return receives;
}

template<class C>bool ParGrid<C>::singleModeWaitAllSends() {
   bool rvalue = true;
   // Reserve enough space for MPI status messages:
   if (MPIstatuses.size() < MPIsendRequests.size()) MPIstatuses.resize(MPIsendRequests.size());
       
   // Wait for all sends to complete:
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),&(MPIstatuses[0]));
   #ifndef NDEBUG
      for (uint i=0; i<MPIsendRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS) rvalue=false;
   #endif
       
   // Free memory:
   MPIsendRequests.clear();
   return rvalue;
}

template<class C>bool ParGrid<C>::singleModeWaitAllSends2() {
   bool rvalue = true;
   
   if (MPIstatuses.size() < MPIsendRequests2.size()) MPIstatuses.resize(MPIsendRequests2.size());
   
   MPI_Waitall(MPIsendRequests2.size(),&(MPIsendRequests2[0]),&(MPIstatuses[0]));
   
   MPIsendRequests2.clear();
   return rvalue;
}

template<class C> uint ParGrid<C>::singleModeWaitSome() {
   if (N_receivesRemaining == 0) return 0;
   
   // Ensure that enough status messages have been allocated:
   if (MPIstatuses.size() < MPIrecvRequests.size()) MPIstatuses.resize(MPIrecvRequests.size());
   
   // Wait until at least one receive that has been posted with singleReceive has completed:
   int receives = 0;
   int readyIndices[MPIrecvRequests.size()];
   while (receives == 0) {
      // Wait for at least one receive to complete:
      //std::cerr << "ParGrid: proc #" << myrank << " receivesRemaining = " << N_receivesRemaining << ' ' << MPIrecvRequests.size() << std::endl;
      if (MPI_Waitsome(MPIrecvRequests.size(),&(MPIrecvRequests[0]),&receives,readyIndices,&(MPIstatuses[0])) != MPI_SUCCESS) {
	 std::cerr << "ParGrid::singleModeWaitSome failed on process #" << myrank << std::endl << std::flush;
	 exit(1);
      }
      if (receives == MPI_UNDEFINED) {
	 std::cerr << "ParGrid::singleModeWaitSome receives is MPI_UNDEFINED!" << std::endl << std::flush;
	 exit(1);
      }

      for (int i=0; i<receives; ++i) readyCells.push_back(localReceiveIDs[readyIndices[i]]);
      N_receivesRemaining -= receives;
   }
   if (N_receivesRemaining == 0) MPIrecvRequests.clear();
   return receives;
}

template<class C> uint ParGrid<C>::singleModeWaitSome2() {
   if (N_receivesRemaining2 == 0) return 0;
   
   // Ensure that enough status messages have been allocated:
   if (MPIstatuses.size() < MPIrecvRequests2.size()) MPIstatuses.resize(MPIrecvRequests2.size());
   
   // Wait until at least one receive has completed:
   int receives = 0;
   int readyIndices[MPIrecvRequests2.size()];
   while (receives == 0) {
      if (MPI_Waitsome(MPIrecvRequests2.size(),&(MPIrecvRequests2[0]),&receives,readyIndices,&(MPIstatuses[0])) != MPI_SUCCESS) {
	 std::cerr << "ParGrid::singleModeWaitSome2 failed on process #" << myrank << std::endl << std::flush;
	 exit(1);
      }
      if (receives == MPI_UNDEFINED) {
	 std::cerr << "ParGrid::singleModeWaitSome2 receives is MPI_UNDEFINED!" << std::endl << std::flush;
	 exit(1);
      }
      
      for (int i=0; i<receives; ++i) readyCells2.push_back(localReceiveIDs2[readyIndices[i]]);
      N_receivesRemaining2 -= receives;
   }
   if (N_receivesRemaining2 == 0) MPIrecvRequests2.clear();
   return receives;
}

template<class C> bool ParGrid<C>::singleReceive(const ID::type& sourceID,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID) {
   std::map<ID::type,int>::const_iterator host = hostProcesses.find(sourceID);
   if (host == hostProcesses.end()) return false;
   if (host->second == myrank) return false;

   //std::cerr << "Proc #" << myrank << " recv src #" << host->second << " tag #" << tag << std::endl;
   
   ++N_receivesRemaining;
   localReceiveIDs.push_back(localID);
   MPIrecvRequests.push_back(MPI_Request());
   if (MPI_Irecv(buffer,byteSize,MPI_BYTE,host->second,tag,MPI_COMM_WORLD,&(MPIrecvRequests.back())) != MPI_SUCCESS) return false;
   return true;
}

template<class C> bool ParGrid<C>::singleReceive2(const int& hostID,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID) {
   //std::map<ID::type,int>::const_iterator host = hostProcesses.find(sourceID);
   //if (host == hostProcesses.end()) return false;
   //if (host->second == myrank) return false;

   ++N_receivesRemaining2;
   localReceiveIDs2.push_back(localID);
   MPIrecvRequests2.push_back(MPI_Request());
   if (MPI_Irecv(buffer,byteSize,MPI_BYTE,hostID,tag,MPI_COMM_WORLD,&(MPIrecvRequests2.back())) != MPI_SUCCESS) return false;
   return true;
}
   
template<class C> bool ParGrid<C>::singleSend(const int& destHost,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID) {
   /*
   std::map<ID::type,int>::const_iterator host = hostProcesses.find(destID);
   if (host == hostProcesses.end()) return false;
   if (host->second == myrank) return false;
   */
   
   //std::cerr << "Proc #" << myrank << " send dest #" << destHost << " tag #" << tag << std::endl;
   
   MPIsendRequests.push_back(MPI_Request());
   //std::cerr << "ParGrid proc #" << myrank << " MPIsendRequests.size() = " << MPIsendRequests.size() << std::endl;
   if (MPI_Issend(buffer,byteSize,MPI_BYTE,destHost,tag,MPI_COMM_WORLD,&(MPIsendRequests.back())) != MPI_SUCCESS) {
      std::cerr << "ParGrid singleSend failed to send data!" << std::endl;
      return false;   
   }
   return true;
}

template<class C> bool ParGrid<C>::singleSend2(const int& destHost,const int& tag,const size_t& byteSize,char* buffer,const ID::type& localID) {
   MPIsendRequests2.push_back(MPI_Request());
   if (MPI_Issend(buffer,byteSize,MPI_BYTE,destHost,tag,MPI_COMM_WORLD,&(MPIsendRequests2.back())) != MPI_SUCCESS) {
      std::cerr << "ParGrid singleSend failed to send data!" << std::endl;
      return false;
   }
   return true;
}

template<class C> void ParGrid<C>::startSingleMode(const int& transfers) {
   if (MPItypeFreed == false) MPI_Type_free(&MPIdataType);
   MPItypeFreed = true;
   
   if (MPIsendRequests.size() > 0) std::cerr << "ParGrid MPIsendRequests size > 0" << std::endl;
   if (MPIrecvRequests.size() > 0) std::cerr << "ParGrid MPIrecvRequests size > 0" << std::endl;
   
   MPIsendRequests.clear();
   MPIrecvRequests.clear();
   MPIstatuses.clear();
   localReceiveIDs.clear();
   N_receivesRemaining = 0;
   readyCells.clear();
   
   // If the number of transfers is known, request enough capacity from 
   // vectors so that send/request posting is faster:
   MPIrecvRequests.reserve(transfers);
   MPIsendRequests.reserve(transfers);
   MPIstatuses.reserve(transfers);
}

template<class C> void ParGrid<C>::startSingleMode2(const int& transfers) {
   if (MPItypeFreed == false) MPI_Type_free(&MPIdataType);
   MPItypeFreed = true;
   
   if (MPIsendRequests2.size() > 0) std::cerr << "ParGrid MPIsendRequests2 size > 0" << std::endl;
   if (MPIrecvRequests2.size() > 0) std::cerr << "ParGrid MPIrecvRequests2 size > 0" << std::endl;
   
   MPIsendRequests2.clear();
   MPIrecvRequests2.clear();
   MPIstatuses.clear();
   localReceiveIDs2.clear();
   N_receivesRemaining2 = 0;
   readyCells2.clear();
   
   MPIrecvRequests2.reserve(transfers);
   MPIsendRequests2.reserve(transfers);
   MPIstatuses.reserve(transfers);
}

/** Synchronizes the information of cells-to-processes assignments among all 
 * MPI processes. Upon successful completion of this function, the member variable 
 * hostProcesses contains for each global ID the rank of the process which has that 
 * cell.
 */
template<class C> bool ParGrid<C>::syncCellAssignments() {
   bool success = true;
   //hostProcesses.clear();
   // Go through every neighbour of every cell. If I do not own the
   // neighbour, add it to remNbrs:
   std::set<ID::type> remNbrs;
   for (typename std::map<ID::type,ParCell<C> >::const_iterator it=localCells.begin(); it!=localCells.end(); ++it) {
      for (std::map<uchar,ID::type>::const_iterator i=it->second.neighbours.begin(); i!=it->second.neighbours.end(); ++i) {
	 const ID::type nbrID = i->second;
	 if (localCells.find(nbrID) == localCells.end()) remNbrs.insert(nbrID);
	 /*
	 // TEMP
	 if (hostProcesses.find(nbrID) == hostProcesses.end()) {
	    std::cerr << "Proc #" << myrank << " nbr cell #" << nbrID << " not in hostProcs" << std::endl;
	 }*/
      }
      /*
      // TEMP
      if (hostProcesses.find(it->first) == hostProcesses.end()) {
	 std::cerr << "Proc #" << myrank << " local cell #" << it->first << " not in hostProcs" << std::endl;
      } else {
	 if (hostProcesses[it->first] != myrank) {
	    std::cerr << "Proc #" << myrank << " local cell #" << it->first << " = " << hostProcesses[it->first];
	    std::cerr << " but I have it!" << std::endl;
	 }
      }
      // END TEMP
      */
      // Add my own cells to hostProcesses:
      hostProcesses[it->first] = myrank;
   }
   
   // Each process sends everyone the IDs of the cells it owns.
   // We can then update hostProcesses based on this information:
   for (int i=0; i<N_processes; ++i) {
      ID::type N_cells;
      
      if (i == myrank) {
	 // Tell everyone how many cells I have:
   	 N_cells = localCells.size();
	 if (MPI_Bcast(&N_cells,1,MPI_Type<ID::type>(),i,MPI_COMM_WORLD) != MPI_SUCCESS) success = false;
	 // Send my cell IDs to everyone:
	 ID::type counter = 0;
	 ID::type* cellIDs = new ID::type[N_cells];
	 for (typename std::map<ID::type,ParCell<C> >::const_iterator it=localCells.begin(); it!=localCells.end(); ++it) {
	    cellIDs[counter] = it->first;
	    ++counter;
	 }
	 if (MPI_Bcast(cellIDs,N_cells,MPI_Type<ID::type>(),i,MPI_COMM_WORLD) != MPI_SUCCESS) success = false;
	 delete cellIDs;
      } else {
	 // Receive the number of cells proc i has:
	 if (MPI_Bcast(&N_cells,1,MPI_Type<ID::type>(),i,MPI_COMM_WORLD) != MPI_SUCCESS) success = false;
	 // Create a buffer and receive cell IDs from proc i:
	 ID::type* cellIDs = new ID::type[N_cells];
	 if (MPI_Bcast(cellIDs,N_cells,MPI_Type<ID::type>(),i,MPI_COMM_WORLD) != MPI_SUCCESS) success = false;
	 
	 // Go through the list of cells, and if one of my boundary cells is on that list,
	 // add an entry to hostProcesses:
	 for (lluint k=0; k<N_cells; ++k) {
	    if (remNbrs.find(cellIDs[k]) != remNbrs.end()) {
	       /*
	       // TEMP
	       if (hostProcesses.find(cellIDs[k]) == hostProcesses.end()) 
		 std::cerr << "Proc #" << myrank << " nbr cell#" << cellIDs[k] << " not in hostProcs" << std::endl;
	       else
		 if (hostProcesses[cellIDs[k]] != i) {
		    std::cerr << "Proc #" << myrank << " hostProcs nbr #" << cellIDs[k] << " = " << hostProcesses[cellIDs[k]];
		    std::cerr << " should equal to " << i << std::endl;
		 }
	       // END TEMP
	       */
	       hostProcesses[cellIDs[k]] = i;
	    }
	 }
	 delete cellIDs;
      }
      barrier();
   }
   remNbrs.clear();
   /*
   // TEMP
   for (int p=0; p<N_processes; ++p) {
      if (myrank == p) {
	 std::cerr << "Proc #" << myrank << " cell information:" << std::endl;
	 for (typename std::map<ID::type,ParCell<C> >::const_iterator it=localCells.begin(); it!=localCells.end(); ++it) {
	    std::cerr << it->first << ": ";
	    for (std::map<uchar,ID::type>::const_iterator i=it->second.neighbours.begin(); i!=it->second.neighbours.end(); ++i) {
	       std::cerr << i->second << "-";
	       if (hostProcesses.find(i->second) == hostProcesses.end()) std::cerr << "? ";
	       else std::cerr << hostProcesses[i->second] << " ";
	    }
	    std::cerr << std::endl;
	 }
      }
      barrier();
   }
   barrier();
   // END TEMP
    */
   return success;
}

/** Exchances unrefined cell indices and lower left corner coordinates of boundary 
 * cells between between neighbouring MPI processes.
 */
template<class C> void ParGrid<C>::syncCellCoordinates() {
   #ifndef NDEBUG
   if (remoteCells.size() != receiveList.size()) {
      std::cerr << "ParGrid::syncCellCoordinates remoteCells.size() = " << remoteCells.size();
      std::cerr << " receiveList.size() = " << receiveList.size() << std::endl;
      exit(1);
   }
   #endif
   
   // Here we send other MPI processes the coordinates of each cell. What is communicated are 
   // unrefInd_i,unrefInd_j,unrefInd_k,xcrd,ycrd,zcrd.
   // First we need to construct an appropriate MPI datatype:
   MPI_Datatype type[6] = {MPI_UNSIGNED,MPI_UNSIGNED,MPI_UNSIGNED,MPI_Type<Real>(),MPI_Type<Real>(),MPI_Type<Real>()};
   int blockLen[6] = {1,1,1,1,1,1};
   MPI_Aint disp[6]; // Byte displacements of each variable
   disp[0] = 0;
   disp[1] = 4;
   disp[2] = 8;
   disp[3] = 12;
   disp[4] = 16;
   disp[5] = 20;

   if (MPItypeFreed == false) MPI_Type_free(&MPIdataType);
   int result = MPI_Type_create_struct(6,blockLen,disp,type,&MPIdataType);
   #ifndef NDEBUG
   if (result != MPI_SUCCESS) {std::cerr << "ParGrid::syncCellCoordinates MPI_Type_create_struct failed!" << std::endl;}
   #endif
   
   result = MPI_Type_commit(&MPIdataType);
   MPItypeFreed = false;
   #ifndef NDEBUG
   if (result != MPI_SUCCESS) {std::cerr << "ParGrid::syncCellCoordinates MPI_Type_commit failed!" << std::endl;}
   #endif
   
   // Create an array for the MPI requests. Make it large enough for both sends and receives:
   if (MPIrecvRequests.size() < receiveList.size()) MPIrecvRequests.resize(receiveList.size());
   if (MPIsendRequests.size() < sendList.size()) MPIsendRequests.resize(sendList.size());

   // Post receives for each remote cell:
   uint counter = 0;
   for (std::map<ID::type,int>::const_iterator it=receiveList.begin(); it!=receiveList.end(); ++it) {      
      void* const buffer = &(remoteCells[it->first]); // Starting address of receive buffer
      const int count = 1;                            // Number of elements in receive buffer
      const int source = it->second;                  // Rank of source MPI process
      const int tag = it->first;                      // Message tag
      MPI_Irecv(buffer,count,MPIdataType,source,tag,MPI_COMM_WORLD,&(MPIrecvRequests[counter]));
      ++counter;
      //std::cerr << "Proc #" << myrank << " recv cell #" << it->first << std::endl;
   }
   
   // Post sends for each boundary cell:
   counter = 0;
   for (std::map<std::pair<ID::type,int>,char>::const_iterator it=sendList.begin(); it!=sendList.end(); ++it) {
      void* const buffer = &(localCells[it->first.first]); // Starting address of send buffer
      const int count = 1;                                 // Number of elements in send buffer
      const int dest = it->first.second;                   // Rank of destination MPI process
      const int tag = it->first.first;                     // Message tag
      MPI_Isend(buffer,count,MPIdataType,dest,tag,MPI_COMM_WORLD,&(MPIsendRequests[counter]));
      ++counter;
      //std::cerr << "Proc #" << myrank << " send cell #" << it->first.first << std::endl;
   }
   // Wait until all data has been received and deallocate memory:
   waitAll();
}

template<class C> bool ParGrid<C>::uncalculatedCells() const {
   if (N_receivesRemaining > 0 || readyCells.size() > 0) return true;
   return false;
}

/** Wait until all sends and receives have been completed.
 * @return True if all communication was successful.
 */
template<class C> bool ParGrid<C>::waitAll() {
   bool rvalue = true;
   // Reserve enough space for MPI status messages:
   MPIstatuses.resize(std::max(MPIrecvRequests.size(),MPIsendRequests.size()));

   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),&(MPIstatuses[0]));
   #ifndef NDEBUG
      for (uint i=0; i<MPIrecvRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS) rvalue=false;
   #endif

   // Wait for all sends to complete:
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),&(MPIstatuses[0]));
   #ifndef NDEBUG
      for (uint i=0; i<MPIsendRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS) rvalue=false;
   #endif
   
   // Free memory:
   MPI_Type_free(&MPIdataType);
   MPItypeFreed = true;
   MPIrecvRequests.clear();
   MPIsendRequests.clear();
   MPIstatuses.clear();
   return rvalue;
}

/** Wait until all remote cell data has arrived from other MPI processes.
 * Vector MPIrecvRequests is cleared here.
 * @return True if all receives were successful.
 */
template<class C> bool ParGrid<C>::waitAllReceives() {
   bool rvalue = true;
   // Reserve enough space for MPI status messages:
   MPIstatuses.resize(MPIrecvRequests.size());
   
   // Wait for all receives to complete:
   MPI_Waitall(MPIrecvRequests.size(),&(MPIrecvRequests[0]),&(MPIstatuses[0]));
   #ifndef NDEBUG
      for (uint i=0; i<MPIrecvRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS) rvalue=false;
   #endif
   
   // Free memory:
   MPIrecvRequests.clear();
   MPIstatuses.clear();
   return rvalue;
}

/** Wait until all cells have been sent to other MPI processes.
 * Vector MPIsendRequests is cleared here.
 * @return True if all sends were successful.
 */
template<class C> bool ParGrid<C>::waitAllSends() {
   bool rvalue = true;
   // Reserve enough space for MPI status messages:
   MPIstatuses.resize(MPIsendRequests.size());
   
   // Wait for all sends to complete:
   MPI_Waitall(MPIsendRequests.size(),&(MPIsendRequests[0]),&(MPIstatuses[0]));
   #ifndef NDEBUG
      for (uint i=0; i<MPIsendRequests.size(); ++i) if (MPIstatuses[i].MPI_ERROR != MPI_SUCCESS) rvalue=false;
   #endif
   
   // Free memory:
   MPIsendRequests.clear();
   MPIstatuses.clear();
   return rvalue;
}

template<class C> bool ParGrid<C>::waitAnyReceive() {
   if (N_receivesRemaining == 0) return false;
   
   // Assure that a status message is allocated:
   if (MPIstatuses.size() == 0) MPIstatuses.resize(1);
   
   // Wait for arriving data until at least one local cell requiring 
   // remote data is ready to be calculated:
   int readyIndex;
   uint N_readyCells = 0;
   while (N_readyCells == 0) {
      // Wait until a receive has completed:
      if (MPI_Waitany(MPIrecvRequests.size(),&(MPIrecvRequests[0]),&readyIndex,&(MPIstatuses[0])) != MPI_SUCCESS) {
	 std::cerr << "ParGrid::testSomeReceives MPI_Waitany failed on process #" << myrank << std::endl;
	 exit(1);
      }
      // Get global ID of the received cell:
      const ID::type globalID = MPIstatuses[0].MPI_TAG;
      
      // Go through all local cells which need the just arrived cell data.
      // Mark that a required remote cell has arrived, and if the local cell is 
      // ready to be calculated add it to readyCells. Note that the received 
      // remote cell may mark several local cells ready.
      for (std::multimap<ID::type,ID::type>::iterator it=remoteToLocal.lower_bound(globalID); it != remoteToLocal.upper_bound(globalID); ++it) {
	 const ID::type localCellID = it->second;
	 ++(localCells[localCellID].N_receivedRemNbrs);
	 // If all neighbours have been received, the cell is ready:
	 if (localCells[localCellID].N_receivedRemNbrs == localCells[localCellID].N_remNbrs) {
	    readyCells.push_back(localCellID);
	    ++N_readyCells;
	 }
      }
      --N_receivesRemaining;
   }
   // If all cells have been received, clear MPIrecvRequests so that it can be reused:
   if (N_receivesRemaining == 0) MPIrecvRequests.clear();
   return true;
}

template<class C> void ParGrid<C>::waitForReceives() const {
   if (myrank == 0) {std::cerr << "ParGrid::waitForReceives" << std::endl;} // TEST
   std::clock_t value = std::clock() + 1.0*CLOCKS_PER_SEC;
   while (std::clock() < value) { }
}

template<class C> bool ParGrid<C>::waitSomeReceives() {
   if (N_receivesRemaining == 0) return false;
   
   // Assure that enough status messages have been allocated:
   if (MPIstatuses.size() < MPIrecvRequests.size()) MPIstatuses.resize(MPIrecvRequests.size());
   
   // Wait until at least one local cell requiring
   // remote data is ready to be calculated:
   uint N_readyCells = 0;
   int receives;
   int readyIndices[MPIrecvRequests.size()];
   ID::type globalID;
   while (N_readyCells == 0) {
      // Wait for at least one receive to complete:
      if (MPI_Waitsome(MPIrecvRequests.size(),&(MPIrecvRequests[0]),&receives,readyIndices,&(MPIstatuses[0])) != MPI_SUCCESS) {
	 std::cerr << "ParGrid::waitSomeReceives failed on process #" << myrank << std::endl;
	 exit(1);
      }
      
      if (receives == MPI_UNDEFINED) {
	 std::cerr << "ParGrid::waitSomeReceives receives is MPI_UNDEFINED!" << std::endl;
	 exit(1);
      }

      for (int i=0; i<receives; ++i) {
	 // Go through all local cells which need the arrived cell data.
	 // Mark that a required remote cell has arrived, and if the local cell is
	 // ready to be calculated add it to readyCells. Note that a single received
	 // remote cell may mark several local cells ready.
	 globalID = MPIstatuses[i].MPI_TAG;
	 for (std::multimap<ID::type,ID::type>::iterator it=remoteToLocal.lower_bound(globalID); it != remoteToLocal.upper_bound(globalID); ++it) {
	    const ID::type localCellID = it->second;
	    ++(localCells[localCellID].N_receivedRemNbrs);
	    // If all neighbours have been received, the cell is ready:
	    if (localCells[localCellID].N_receivedRemNbrs == localCells[localCellID].N_remNbrs) {
	       readyCells.push_back(localCellID);
	       ++N_readyCells;
	    }
	 }
      }
      N_receivesRemaining -= receives;
   }   
   // If all cells have been received, clear MPIrecvRequests so that it can be reused:
   if (N_receivesRemaining == 0) MPIrecvRequests.clear();
   return true;
}

template<class C>
void ParGrid<C>::writeLoadDistribution() {
   uint N_localCells[N_processes];
   uint N_remoteCells[N_processes];
   uint N_sends[N_processes];
   uint N_receives[N_processes];
   
   uint N_myLocals = localCells.size();
   uint N_myRemotes = remoteCells.size();
   uint N_mySends = sendList.size();
   uint N_myReceives = receiveList.size();
   
   MPI_Gather(&N_myLocals,  1,MPI_UNSIGNED,N_localCells, 1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
   MPI_Gather(&N_myRemotes, 1,MPI_UNSIGNED,N_remoteCells,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
   MPI_Gather(&N_mySends,   1,MPI_UNSIGNED,N_sends,      1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
   MPI_Gather(&N_myReceives,1,MPI_UNSIGNED,N_receives,   1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
   
   if (myrank == 0) {
      std::cerr << "# local cells\tremote cells\tsends\treceives\tsends+receives" << std::endl;
      for (int i=0; i<N_processes; ++i) {
	 std::cerr << N_localCells[i] << '\t' << N_remoteCells[i] << '\t' << N_sends[i] << '\t' << N_receives[i] << '\t';
	 std::cerr << N_sends[i]+N_receives[i] << std::endl;
      }
   }
}

// *************************************************************
// ***************** ZOLTAN CALLBACK FUNCTIONS *****************
// *************************************************************

// GENERIC CALLBACK FUNCTIONS
// The functions below are required for everything.

/** Definition for Zoltan callback function ZOLTAN_NUM_OBJ_FN. This function 
 * is required to use Zoltan. The purpose is to tell Zoltan how many cells 
 * are currently assigned to this process.
 * @param parGridPtr A pointer to ParGrid.
 * @param rcode The return code. Upon success this should be ZOLTAN_OK.
 * @return The number of cells assigned to this process.
 */
template<class C>
int ParGrid<C>::getNumberOfLocalCells(void* parGridPtr,int* rcode) {     
   ParGrid<C>* ptr = reinterpret_cast<ParGrid<C>*>(parGridPtr);
   *rcode = ZOLTAN_OK;
   return ptr->localCells.size();
}

/** Definition for Zoltan callback function ZOLTAN_OBJ_LIST_FN. This function 
 * is required to use Zoltan. The purpose is to tell Zoltan the global and local 
 * IDs of the cells assigned to this process, as well as their weights.
 * @param parGridPtr A pointer to ParGrid.
 * @param N_globalIDs The number of array entries used to describe one global ID.
 * @param N_localIDs The number of array entries used to describe one local ID.
 * @param globalIDs An array which is to be filled with the global IDs of the cells 
 * currently assigned to this process. This array has been allocated by Zoltan.
 * @param localIDs An array which is to be filled with the local IDs of the cells 
 * currently assigned to this process. This array has been allocated by Zoltan.
 * @param N_weights
 * @param cellWeights An array which is to be filled with the cell weights.
 * @param rcode The return code. Upon success this should be ZOLTAN_OK.
 */
template<class C>
void ParGrid<C>::getLocalCellList(void* parGridPtr,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
				  ZOLTAN_ID_PTR localIDs,int N_weights,float* cellWeights,int* rcode) {
   ParGrid<C>* parGrid = reinterpret_cast<ParGrid<C>*>(parGridPtr);

   int i=0;
   if (N_weights == 0) { // No cell weights
      for (typename std::map<ID::type,ParCell<C> >::const_iterator it = parGrid->localCells.begin(); it!=localCells.end(); ++it) {
	 globalIDs[i] = it->first;
	 ++i;
      }
   } else { // Cell weights are used
      for (typename std::map<ID::type,ParCell<C> >::const_iterator it = parGrid->localCells.begin(); it!=localCells.end(); ++it) {
	 globalIDs[i] = it->first;
	 cellWeights[i] = calculateCellWeight(it->first);
	 ++i;
      }
   }
   *rcode = ZOLTAN_OK;
}


// GEOMETRY-BASED LOAD BALANCING:                                                          
// The functions below are for geometry-based load balancing
// functions (BLOCK,RCB,RIB,HSFC,Reftree).                         
// -------------------------------------------------------------- 

/** Definition for Zoltan callback function ZOLTAN_NUM_GEOM_FN. This function is 
 * required for geometry-based load balancing. The purpose is 
 * to tell Zoltan the dimensionality of the grid. Here we use three-dimensional 
 * grid.
 * @param parGridPtr A pointer to ParGrid.
 * @param rcode The return code. Upon success this should be ZOLTAN_OK.
 * @return The number of physical dimensions.
 */
template<class C>
int ParGrid<C>::getGridDimensions(void* parGridPtr,int* rcode) {   
   *rcode = ZOLTAN_OK;
   return 3;
}

/** Definition for Zoltan callback function ZOLTAN_GEOM_FN. This function is required for
 * geometry-based load balancing. The purpose of this function is to tell 
 * Zoltan the physical x/y/z coordinates of a given cell.
 * @param parGridPtr Pointer to ParGrid.
 * @param N_globalEntries The size of array globalID.
 * @param N_localEntries The size of array localID.
 * @param globalID The global ID of the cell whose coordinates are queried.
 * @param localID The local ID of the cell whose coordinates are queried.
 * @param geometryData Array where the coordinates should be written. Zoltan has reserved this array, its size 
 * is determined by the getGridDimensions function.
 * @param rcode The return code. Upon success this should be ZOLTAN_OK.
 */
template<class C> 
void ParGrid<C>::getCellCoordinates(void* parGridPtr,int N_globalEntries,int N_localEntries,ZOLTAN_ID_PTR globalID,
				    ZOLTAN_ID_PTR localID,double* geometryData,int* rcode) {
   ParGrid<C>* parGrid = reinterpret_cast<ParGrid<C>*>(parGridPtr);

   // Try to find the given cell:
   typename std::map<ID::type,ParCell<C> >::const_iterator it = parGrid->localCells.find(globalID[0]);
   #ifndef NDEBUG
   if (it == parGrid->localCells.end()) {
      *rcode = ZOLTAN_FATAL;
      std::cerr << "ParGrid getCellCoordinates queried non-existing local cell." << std::endl;
      return;
   }
   #endif
   geometryData[0] = (it->second).xcrd;
   geometryData[1] = (it->second).ycrd;
   geometryData[2] = (it->second).zcrd;
   rcode = ZOLTAN_OK;
}

// GRAPH-BASED LOAD BALANCING

/** Definition for Zoltan callback function ZOLTAN_NUM_EDGES_FN. This function is required 
 * for graph-based load balancing. The purpose is to tell how many edges a given cell has, i.e. 
 * how many neighbours it has to share data with.
 * @param parGridPtr Pointer to ParGrid.
 * @param N_globalIDs The size of array globalID.
 * @param N_localIDs The size of array localID.
 * @param globalID The global ID of the cell whose edges are queried.
 * @param localID The local ID of the cell whose edges are queried.
 * @param rcode The return code. Upon success this should be ZOLTAN_OK.
 * @return The number of edges the cell has. For three-dimensional box grid this is between 3 and 6, 
 * depending on if the cell is on the edge of the simulation volume.
 */
template<class C>
int ParGrid<C>::getNumberOfEdges(void* parGridPtr,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
		     ZOLTAN_ID_PTR localID,int* rcode) {

   ParGrid<C>* parGrid = reinterpret_cast<ParGrid<C>*>(parGridPtr);
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID[0]);
   #ifndef NDEBUG
   if (it == localCells.end()) {
      *rcode = ZOLTAN_FATAL;
      std::cerr << "ParGrid::getNumberOfEdges Zoltan tried to query nonexisting cell!" << std::endl;
      return 0;
   }
   #endif
   int edges = it->second.neighbours.size();
   
   *rcode = ZOLTAN_OK;
   return edges;
}

/** Definition for Zoltan callback function ZOLTAN_EDGE_LIST_FN. This function is required 
 * for graph-based load balancing. The purpose is to give the global IDs of each neighbour of 
 * a given cell, as well as the ranks of the MPI processes which have the neighbouring cells.
 * MAKE SURE THAT hostProcesses IS UP TO DATE BEFORE CALLING THIS!
 * @param parGridPtr A pointer to ParGrid.
 * @param N_globalIDs The size of array globalID.
 * @param N_localIDs The size of array localID.
 * @param globalID The global ID of the cell whose neighbours are queried.
 * @param localID The local ID of the cell whose neighbours are queried.
 * @param nbrGlobalIDs An array where the global IDs of the neighbours are written. Note that
 * Zoltan has already allocated this array based on a call to getNumberOfEdges function.
 * @param nbrHosts For each neighbour, the rank of the MPI process which has the cell.
 * @param N_weights The size of array weight.
 * @param weight The weight of each edge.
 * @param rcode The return code. Upon success should be ZOLTAN_OK.
 */
template<class C>
void ParGrid<C>::getEdgeList(void* parGridPtr,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			     ZOLTAN_ID_PTR localID,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			    int N_weights,float* weight,int* rcode) {
   ParGrid<C>* parGrid = reinterpret_cast<ParGrid<C>*>(parGridPtr);
   typename std::map<ID::type,ParCell<C> >::const_iterator it = localCells.find(globalID[0]);
   #ifndef NDEBUG
   if (it == localCells.end()) {
      *rcode = ZOLTAN_FATAL;
      std::cerr << "ParGrid::getEdgeList Zoltan tried to query nonexisting cell!" << std::endl;
      return;
   }
   #endif
   int index = 0;
   if (N_weights == 0) { // No edge weights
      for (std::map<uchar,ID::type>::const_iterator i=it->second.neighbours.begin(); i!=it->second.neighbours.end(); ++i) {
	 nbrGlobalIDs[index] = i->second;
	 nbrHosts[index] = hostProcesses[i->second];
	 ++index;
      }
   } else { // Edge weights are used
      if (N_weights > 1) std::cerr << "ParGrid::getEdgeList returns incorrect number of weights!" << std::endl;
      for (std::map<uchar,ID::type>::const_iterator i=it->second.neighbours.begin(); i!=it->second.neighbours.end(); ++i) {
	 nbrGlobalIDs[index] = i->second;
	 nbrHosts[index] = hostProcesses[i->second];
	 weight[index] = calculateEdgeWeight(it->first);
	 ++index;
      }
   }
   *rcode = ZOLTAN_OK;
}

// HYPERGRAPH-BASED LOAD BALANCING:
// --------------------------------------------------------------

/** Definition for Zoltan callback function ZOLTAN_HG_SIZE_CS_FN. This function is required 
 * for hypergraph-based load balancing. The purpose is to tell Zoltan which hypergraph format 
 * is used (ZOLTAN_COMPRESSED_EDGE or ZOLTAN_COMPRESSED_VERTEX), how many hyperedges and 
 * vertices there will be, and how many pins.
 * @param parGridPtr A pointer to ParGrid.
 * @param N_lists The total number of vertices or hyperedges (depending on the format) 
 * is written to this variable.
 * @param N_pins The total number of pins (connections between vertices and hyperedges) is 
 * written to this variable.
 * @param format The chosen hyperedge storage format is written to this variable.
 * @param rcode The return code. Upon success should be ZOLTAN_OK.
 */
template<class C>
void ParGrid<C>::getNumberOfHyperedges(void* parGridPtr,int* N_lists,int* N_pins,int* format,int* rcode) {
   *N_lists = localCells.size();
   *format = ZOLTAN_COMPRESSED_VERTEX;
   
   // Calculate the total number of pins:
   unsigned int totalNumberOfPins = 0;
   for (typename std::map<ID::type,ParCell<C> >::const_iterator it=localCells.begin(); it!=localCells.end(); ++it) {
      // Every cell has its own hyperedge:
      ++totalNumberOfPins;
      // Every existing neighbour belongs to cell's hyperedge:
      totalNumberOfPins += it->second.neighbours.size();
   }
   *N_pins = totalNumberOfPins;
   *rcode = ZOLTAN_OK;
}

/** Definition for Zoltan callback function ZOLTAN_HG_CS_FN. This function is required for 
 * hypergraph-based load balancing. The purpose is to give Zoltan the hypergraph in a compressed format.
 * @param parGridPtr A pointer to ParGrid.
 * @param N_globalIDs The size of globalID.
 * @param N_vtxedges The number of entries that need to be written to vtxedge_GID.
 * @param N_pins The number of pins that need to be written to pin_GID.
 * @param format The format that is used to represent the hypergraph, either ZOLTAN_COMPRESSED_EDGE or ZOLTAN_COMPRESSED_VERTEX.
 * @param vtxedge_GID An array where the hypergraph global IDs are written into.
 * @param vtxedge_ptr An array where, for each hyperedge, an index into pin_GID is given from where the pins for that 
 * hyperedge are given.
 * @param pin_GID An array where the pins are written to.
 * @param rcode The return code. Upon success should be ZOLTAN_OK.
 */
template<class C>
void ParGrid<C>::getHyperedges(void* parGridPtr,int N_globalIDs,int N_vtxedges,int N_pins,int format,ZOLTAN_ID_PTR vtxedge_GID,
			       int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* rcode) {
   uint edgeCounter = 0;
   uint pinCounter = 0;
   if (format == ZOLTAN_COMPRESSED_VERTEX) {
      // Go through every local cell:
      for (typename std::map<ID::type,ParCell<C> >::const_iterator it=localCells.begin(); it!=localCells.end(); ++it) {
	 vtxedge_GID[edgeCounter] = it->first;  // The hyperedge has the same global ID as the cell
	 vtxedge_ptr[edgeCounter] = pinCounter; // An index into pin_GID where the pins for this cell are written
	 pin_GID[pinCounter] = it->first;       // Every cell belong in its own hyperedge
	 ++pinCounter;
	 // Add pin to every existing neighbour:
	 for (std::map<uchar,ID::type>::const_iterator i=it->second.neighbours.begin(); i!=it->second.neighbours.end(); ++i) {
	    pin_GID[pinCounter] = i->second;
	    ++pinCounter;
	 }
	 ++edgeCounter;
      }
   } else {
      std::cerr << "ParGrid: ZOLTAN_COMPRESSED_EDGE not implemented, exiting." << std::endl;
      exit(1);
   }
   *rcode = ZOLTAN_OK;
}

/** Definition for Zoltan callback function ZOLTAN_HG_SIZE_EDGE_WTS_FN. This is an optional function 
 * for hypergraph-based load balancing. The purpose is to tell Zoltan how many hyperedges will have 
 * a weight factor. Here we give a weight to each hyperedge.
 * @param parGridPtr A pointer to ParGrid.
 * @param N_edges A parameter where the number of weight-supplying hyperedges is written into.
 * @param rcode The return code. Upon success should be ZOLTAN_OK.
 */
template<class C>
void ParGrid<C>::getNumberOfHyperedgeWeights(void* parGridPtr,int* N_edges,int* rcode) {
   *N_edges = localCells.size();
   *rcode = ZOLTAN_OK;
}

/** Definition for Zoltan callback function ZOLTAN_HG_EDGE_WTS_FN. This is an optional function 
 * for hypergraph-based load balancing. The purpose is to tell Zoltan the weight of each hyperedge.
 * @param parGridPtr A pointer to ParGrid.
 * @param N_globalIDs The size of edgeGlobalID entry.
 * @param N_localIDs The size of edgeLocalID entry.
 * @param N_edges The number of hyperedge weights that need to be written to edgeWeights.
 * @param N_weights Number of weights per hyperedge that need to be written.
 * @param edgeGlobalID An array where the global IDs of each weight-supplying hyperedge are to be written.
 * @param edgeLocalID An array where the local IDs of each weight-supplying hyperedge are to be written.
 * This array can be left empty.
 * @param edgeWeights An array where the hyperedge weights are written into.
 * @param rcode The return code. Upon success should be ZOLTAN_OK.
 */
template<class C>
void ParGrid<C>::getHyperedgeWeights(void* parGridPtr,int N_globalIDs,int N_localIDs,int N_edges,int N_weights,
				     ZOLTAN_ID_PTR edgeGlobalID,ZOLTAN_ID_PTR edgeLocalID,float* edgeWeights,int* rcode) {
   uint counter = 0;
   for (typename std::map<ID::type,ParCell<C> >::const_iterator it=localCells.begin(); it!=localCells.end(); ++it) {
      edgeGlobalID[counter] = it->first;
      edgeWeights[counter] = calculateHyperedgeWeight(it->first);
      ++counter;
   }
   *rcode = ZOLTAN_OK;
}

template<class C>
int ParGrid<C>::getNumberOfHierarchicalLevels(void* parGridPtr,int* rcode) {   
   *rcode = ZOLTAN_OK;
   return N_hierarchicalLevels;
}

template<class C>
int ParGrid<C>::getHierarchicalPartNumber(void* parGridPtr,int level,int* rcode) {
   int rvalue;
   *rcode = ZOLTAN_OK;
   switch (level) {
    case 0:
      // Rank of superpartition containing this process:
      rvalue = myrank / N_processesPerPart;
      break;
    case 1:
      // Rank of this process within the superpartition:
      rvalue = myrank % N_processesPerPart;
      break;
    default:
      // Incorrect number of partitioning levels:
      *rcode = ZOLTAN_FATAL;
      break;
   }
   return rvalue;
}

template<class C>
void ParGrid<C>::getHierarchicalParameters(void* parGridPtr,int level,Zoltan_Struct* zs,int* rcode) {
   *rcode = ZOLTAN_OK;
   switch (level) {
    case 0:
      // Refinement parameters for superpartitions:
      Zoltan_Set_Param(zs,"LB_METHOD","HYPERGRAPH");
      Zoltan_Set_Param(zs,"PHG_CUT_OBJECTIVE","CONNECTIVITY");
      Zoltan_Set_Param(zs,"OBJ_WEIGHT_DIM","1");
      Zoltan_Set_Param(zs,"EDGE_WEIGHT_DIM","0");
      Zoltan_Set_Param(zs,"IMBALANCE_TOL",imbalanceTolerance.c_str());
      break;
    case 1:
      // Refinement parameters for partitioning a superpartition:
      Zoltan_Set_Param(zs,"LB_METHOD","HYPERGRAPH");
      Zoltan_Set_Param(zs,"PHG_CUT_OBJECTIVE","CONNECTIVITY");
      Zoltan_Set_Param(zs,"OBJ_WEIGHT_DIM","1");
      Zoltan_Set_Param(zs,"EDGE_WEIGHT_DIM","0");
      Zoltan_Set_Param(zs,"IMBALANCE_TOL",imbalanceTolerance.c_str());
      break;
    default:
      *rcode = ZOLTAN_FATAL;
      break;
   }
}

#endif

