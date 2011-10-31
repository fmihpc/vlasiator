/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include <stdint.h>
#include <vector>
#include <mpi.h>
#include "definitions.h"
#include "cell_spatial.h"

bool buildSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,
		      creal& zmin,creal& dx,creal& dy,creal& dz,
		      const bool& isRemote);

bool buildGrid(MPI_Comm,const int& MASTER_RANK);

// NOTE FOR VLASOV SIMULATION NEIGHBOUR TYPE IDENTIFIERS
// 
// The type IDs are defined as follows:
// 0-7   are for -x neighbours, where  0 = is for unrefined -x nbr
// 8-15          +x                    8                    +x
// 16-23         -y                   16                    -y
// 24-31         +y                   24                    +y
// 32-39         -z                   32                    -z
// 40-47         +z                   40                    +z   

namespace VirtualCell {
   typedef uint64_t ID;        /**< Datatype which is used for obtaining cell global IDs.*/
}

/** Base class for GridBuilder. GridBuilder is used to 
 * create spatial simulation cells transparently.
 */
class GridBuilder {
 public:
   GridBuilder();
   virtual ~GridBuilder();

   /** Add a request to obtain the average values of distribution function, parameters, and 
    * neighbour lists of all velocity blocks of all spatial cell allocated to this process. 
    * When process has placed all its requests, 
    * GridBuilder::processCellBlockDataRequests and GridBuilder::waitCellBlockDataRequests 
    * need to be called (in that order).
    * @param totalCells Number of spatial cells this process has.
    * @param blockOffset Offset into arrays where this process reads its velocity block data.
    * This parameter is mainly required for restarts.
    * @param cellIDs Global IDs of spatial cells allocated to this process.
    * @param blocksPerCell Number of velocity blocks in each spatial cell.
    * @param avgsBuffer Pointer to arrays in which the velocity block data is to be written.
    * @param blockParams Pointers to arrays in which velocity block parameters are to be written.
    * @param nbrsVels Pointers to arrays in which velocity block neighbour lists are to be written.
    * @return If true, the request added successfully.
    * @see GridBuilder::processCellBlockDataRequests
    * @see GridBuilder::waitCellBlockDataRequests
    */
   virtual bool addCellBlockDataRequests(VirtualCell::ID& totalCells,VirtualCell::ID& blockOffset,VirtualCell::ID* cellIDs,uint* blocksPerCell,
					 Real** avgsBuffer,Real** blockParamsBuffer,uint** nbrsVelBuffer) = 0;
   
   /** Add a request to obtain the number of velocity blocks in each spatial cell allocated 
    * to this process. When a process has placed all its requests,
    * GridBuilder::processCellBlockNumberRequests and GridBuilder::waitCellBlockNumberRequests
    * need to be called (in that order).
    * @param totalCells Number of spatial cells this process has.
    * @param cellOffset Offset into arrays where this process reads its per-spatial-cell data.
    * This parameter is mainly required for restarts.
    * @param cellID Global IDs of spatial cell allocated to this process.
    * @param N_blocks Array in which the number of velocity blocks in each cell is to be written.
    * @return If true, the request was added successfully.
    * @see GridBuilder::processCellBlockNumberRequests
    * @see GridBuilder::waitCellBlockNumberRequests
    */
   virtual bool addCellBlockNumberRequests(VirtualCell::ID& totalCells,VirtualCell::ID& cellOffset,VirtualCell::ID* cellIDs,uint* N_blocks) = 0;
   
   /** Add a request to obtain spatial cell coordinates and neighbour data.
    * When process has placed all its requests,
    * GridBuilder::processCellNbrRequests and GridBuilder::waitCellNbrRequests
    * need to be called (in that order).
    * 
    * NOTE: Cell coordinates and sizes are requested here only because some partitioners 
    * are based on coordinates. If partitioners such as Zoltans RCB or HSFC are dumped, the 
    * cell coordinates and sizes can be requested along with cell parameters. The information 
    * obtained from this step in the grid building process can be used to do an initial 
    * partitioning of cells (without cells actually containing velocity block data) from a 
    * horribly bad partition into a better one.
    * @param totalCells Number of cells allocated to this process.
    * @param totalNbrs Neighbour count over all cells, i.e. sum(i) sum(j) cell[i].neighbours[j]. 
    * @param cellOffset Offset into arrays where this process reads its per-spatial-cell data. 
    * This parameter is mainly required for restarts.
    * @param nbrOffset Offset into arrays where this process reads its spatial cell neighbour data.
    * This parameter is mainly required for restarts.
    * @param cellID Global IDs of all spatial cells allocated to this process.
    * @param nbrsPerCell Array containing the number of spatial neighbours each cell has.
    * @param coords A buffer in which the cells' coordinates and sizes will be written, the 
    * size of this array must be 6*totalCells for 3D data.
    * @param nbrIDs A buffer in which the global IDs of each cell's neighbours are to be written.
    * @param nbrTypes A buffer in which the neighbour type identifiers of each cell are to be written.
    * @return If true, the request was added successfully.
    * @see GridBuilder::processCellNbrRequests
    * @see GridBuilder::waitCellNbrRequests
    */
   virtual bool addCellNbrRequests(VirtualCell::ID& totalCells,VirtualCell::ID& totalNbrs,VirtualCell::ID& cellOffset,
				   VirtualCell::ID& nbrOffset,VirtualCell::ID* cellIDs,
				   uchar* nbrsPerCell,Real* coords,VirtualCell::ID* nbrIDs,uchar* nbrTypes) = 0;

   /** Add a request to obtain parameters of a given spatial cell.
    * After process has placed all its requests,
    * GridBuilder::processCellParamsRequests and GridBuilder::waitCellParamsRequests
    * need to be called (in that order). 
    * @param totalCells Number of spatial cells this process has.
    * @param cellOffset Offset into arrays where this process reads its per-spatial-cell data.
    * This parameter is mainly required for restarts.
    * @param cellID Global IDs of spatial cell allocated to this process.
    * @param cellParams A buffer in which the parameters of each cell are to be written. 
    * The size of this buffer is SIZE_CELLPARAMS*totalCells.
    * @return If true, the request was added successfully.
    * @see GridBuilder::processCellParamsRequests
    * @see GridBuilder::waitCellParamsRequests
    */
   virtual bool addCellParamsRequests(VirtualCell::ID& totalCells,VirtualCell::ID& cellOffset,VirtualCell::ID* cellIDs,Real* cellParams) = 0;

   /** Query if parallel grid is allowed to do a load balance before 
    * all data has been created by a GridBuilder. The default behaviour is 
    * to allow a load balance.
    * @return If true, grid is allowed to do a load balance. 
    */
   virtual bool doInitialLoadBalance();
   
   /** Deallocate all memory reserved by GridBuilder. This function is called 
    * just before GridBuilder is deleted.
    * @return If true, GridBuilder was finalized successfully.
    */
   virtual bool finalize() =0;

   /** Get the global IDs of all cells GridBuilder will create, as well as 
    * the numbers of spatial neighbours each cell has. This function is meant to be 
    * called by master process, which will then distribute the cells to other 
    * processes based on this information. For any other process this 
    * function will immediately return and the given buffers contain invalid data.
    * After master process has distributed cells, every process should allocate memory 
    * based on the data received from master, and then request the 
    * cell neighbour data by calling GridBuilder::addCellNbrRequest for every cellID 
    * it has received.
    * 
    * NOTE: GridBuilder may or may not read the data from an input file, thus this 
    * function should only work successfully on master process.
    * @param cellIDs Vector in which the global IDs of all existing cells are to be written.
    * @param N_nbrs Vector in which the number of spatial neighbours each cell has are to be writtten.
    * @return If true, data was successfully written to buffers. Return value false may 
    * also indicate that this function was called by a process other than the master.
    */
   virtual bool getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs) = 0;

   /** Get the coordinates and sizes, global IDs of neighbours, and type identifiers of 
    * each neighbour of the given spatial cell. 
    * 
    * NOTE: GridBuilder may or may not read the data from an input file, thus this
    * function should only work successfully on master process.
    * @param N_cells Number of spatial cells allocated to this process.
    * @param cellIDs Global IDs of cells.
    * @param spatNbrIDs Array in which global IDs of each cell's neighbours are to be written.
    * @param spatNbrTypes Array in which type identifiers of each cell's neighbours are to be written.
    * @return If true, neighbour data was successfully written to given buffers. Return value 
    * false may also indicate that this function was called by a process other than the master.
    */
   //virtual bool getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* nbrTypes) = 0;
   
   /** Request the value of an input parameter. This function is provided in 
    * cases where the user might want to know the values of grid-related input 
    * parameters. Master process should distribute the values of such parameters 
    * to all processes.
    * @param parameterName The name of the parameter.
    * @param value A string in which the value of the parameter is to be written.
    * @return If true, GridBuilder returned a valid parameter value. If false, 
    * the parameter was not found.
    */
   virtual bool getParameter(const std::string& parameterName,std::string& value) =0;
   
   /** Query the total number of cells the GridBuilder creates. This 
    * function is meant to be called by the master process, which uses the 
    * information to distribute the cells to other processes.
    * 
    * NOTE: GridBuilder may or may not read the data from an input file, thus this
    * function should only work successfully on master process.
    * @param N_cells A variable in which the total number of existing cells is to be written.
    * @return If true, GridBuilder wrote a valid value to N_cells.
    */
   virtual bool getTotalNumberOfCells(VirtualCell::ID& N_cells) =0;
   
   /** Initialize the GridBuilder. If the GridBuilder, for example, reads the 
    * cells from a file, this is the correct place to open that file.
    * This function gets called before any other GridBuilder member functions 
    * are called, except for the constructor.
    * %Parameters and options read from the simulation input files, via calls to 
    * class Parameters member functions, should also be defined and read here. 
    * GridBuilder input file parameters must be defined in a section called "gridbuilder".
    * @return If true, GridBuilder initialized successfully. If any process 
    * returns false, program execution should be aborted.
    */
   virtual bool initialize(MPI_Comm comm,const int& MASTER_RANK) =0;

   /** Process the requests for cell velocity block data. If the data is, for example,
    * read from a file, then only the master process should do file I/O and
    * send the data to processes which requested the data. On other processes
    * this function must immediately return, so that the processes can begin to
    * wait for data to arrive. GridBuilder::waitCellBlockDataRequests needs to be
    * called immediately after this function returns.
    * @return If true, this process processed the requests successfully
    * Program execution should be aborted if any process returned false.
    * @see GridBuilder::waitCellBlockDataRequests
    */
   virtual bool processCellBlockDataRequests() = 0;
   
   /** Process the requests for numbers of velocity blocks in spatial cells. 
    * If the data is, for example,
    * read from a file, then only the master process should do file I/O and
    * send the data to processes which requested the data. On other processes
    * this function must immediately return, so that the processes can begin to
    * wait for data to arrive. GridBuilder::waitCellBlockNumberRequests needs to be
    * called immediately after this function returns.
    * @return If true, this process processed the requests successfully
    * Program execution should be aborted if any process returned false.
    * @see GridBuilder::waitCellBlockNumberRequests
    */
   virtual bool processCellBlockNumberRequests() = 0;
   
   /** Process the requests for cell neighbour data. If the data is, for example,
    * read from a file, then only the master process should do file I/O and 
    * send the data to processes which requested the data. On other processes 
    * this function must immediately return, so that the processes can begin to 
    * wait for data to arrive. GridBuilder::waitCellNbrRequests needs to be 
    * called immediately after this function returns.
    * @return If true, this process processed the requests successfully.
    * Program execution should be aborted if any process returned false.
    * @see GridBuilder::waitCellNbrRequests
    */
   virtual bool processCellNbrRequests() = 0;

   /** Process the requests for spatial cell parameters. If the data is, for example,
    * read from a file, then only the master process should do file I/O and 
    * send the data to processes which requested the data. On other processes
    * this function must immediately return, so that the processes can begin to
    * wait for data to arrive. GridBuilder::waitCellParamsRequests needs to be
    * called immediately after this function returns.
    * @return If true, this process processed the requests successfully.
    * Program execution should be aborted if any process returned false.
    * @see GridBuilder::processCellParamsRequests
    */
   virtual bool processCellParamsRequests() = 0;

   /** Wait until all requests for spatial cell velocity block data have been fulfilled. After
    * this function returns, the buffers given in GridBuilder::addCellBlockDataRequest
    * contain valid data.
    * @return If true, all requested data was received successfully.
    * Program execution should be aborted if any process returned false.
    */
   virtual bool waitCellBlockDataRequests() = 0;
   
   /** Wait until all requests for numbers of velocity blocks in spatial cells have been fulfilled. After
    * this function returns, the buffers given in GridBuilder::addCellBlockNumberRequest
    * contain valid data.
    * @return If true, all requested data was received successfully.
    * Program execution should be aborted if any process returned false.
    */
   virtual bool waitCellBlockNumberRequests() = 0;
   
   /** Wait until all requests for spatial cell neighbour data have been fulfilled. After 
    * this function returns, the buffers given in GridBuilder::addCellNbrRequest
    * contain valid data.
    * @return If true, all requested data was received successfully. 
    * Program execution should be aborted if any process returned false.
    */
   virtual bool waitCellNbrRequests() = 0;
   
   /** Wait until all requests for spatial cell parameter data have been fulfilled. After
    * this function returns, the buffers given in GridBuilder::addCellParamsRequest
    * contain valid data.
    * @return If true, all requested data was received successfully.
    * Program execution should be aborted if any process returned false.
    */
   virtual bool waitCellParamsRequests() = 0;
   
   // ***************************************************************************
   // ***** MEMBER FUNCTIONS FOR GETTING THE VALUES OF MANDATORY PARAMETERS *****
   // ***************************************************************************
   
   Real getSpeciesCharge() const {return q;}
   Real getSpeciesMass() const {return m;}
   Real getDt() const {return dt;}
   Real getTmin() const {return t_min;}
   luint getTimestep() const {return timestep;}
   luint getMaxTimesteps() const {return max_timesteps;}
 protected:
   
   // *********************************************************
   // *****  ALL GRIDBUILDERS ARE REQUIRED TO BE ABLE TO  *****
   // ***** PROVIDE VALUES FOR THE PARAMETERS GIVEN BELOW *****
   // ***** NOTE:   THESE VALUES ARE ONLY REQUESTED ON    *****
   // *****                MASTER PROCESS                 *****
   // *********************************************************
   
   Real   q;            /**< Charge of simulated particle species, in Coulombs.*/
   Real   m;            /**< Mass of simulated particle species, in kilograms.*/
   Real   dt;           /**< Value of timestep, in seconds.*/
   Real   t_min;        /**< Value of time, in seconds, at timestep 0.*/
   luint timestep;      /**< Timestep when the simulation is started. Non-zero when simulation is re-started.*/
   luint max_timesteps; /**< Maximum number of timesteps taken by the simulation.*/
};

// Typedef of a function pointer
typedef GridBuilder* (*GridBuilderCreator)();

/** A small object factory for creating a GridBuilder. The user must 
 * register a valid GridBuilder into GridBuilderFactory when Vlasov 
 * simulation is initialized. The given GridBuilder is then used to 
 * create spatial cells. Typically the registration process is handled by 
 * the object file containing the GridBuilder (with a static initializer), 
 * and the user just links that object file with the rest of the simulation.
 */ 
class GridBuilderFactory {
 public:
   GridBuilderFactory();
   ~GridBuilderFactory();
   
   static GridBuilder* createBuilder();
   static bool deleteBuilder(GridBuilder*& gb);
   static bool registerBuilder(GridBuilderCreator gbc);
   
 private:
   static GridBuilderCreator gbc; /**< Pointer to function which returns a pointer to a new GridBuilder.*/
};

#endif

