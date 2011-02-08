#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include <vector>
#include "definitions.h"
#include "cell_spatial.h"

bool buildSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,
		      creal& zmin,creal& dx,creal& dy,creal& dz,
		      const bool& isRemote);

#ifdef PARGRID
   #include "pargrid.h"
   bool buildGrid(ParGrid<SpatialCell>& mpiGrid);
#else

#endif

/** Base class for GridBuilder. GridBuilder is used to 
 * create spatial simulation cells transparently.
 */
class GridBuilder {
 public:
   GridBuilder();
   virtual ~GridBuilder();

   /** Deallocate all memory reserved by GridBuilder. This function is called 
    * just before GridBuilder is deleted.
    * @return If true, GridBuilder was finalized successfully.
    */
   virtual bool finalize() =0;
   
   /** Request a new spatial cell. Each cell must have a unique ID calculated by 
    * this function. The actual value of the ID is not important, as long as it is unique.
    * This function gets called until it returns false, in 
    * which case all cells have been created.
    * @param cellID A variable where the unique global ID of the cell is to be written.
    * @param coords An array where the cell coordinates are to be written.
    * @param sizes An array where the cell size per coordinate direction are to be written.
    * @param nbrs A vector where the cell's neighbours' global IDs and type identifiers are to be written.
    * Type identifiers are used by the solver to, for example, deduce in which coordinate 
    * direction (relative to the cell) the neighbour is located.
    * @return If true, GridBuilder returned a valid cell. If false, all cells have been given.
    */
   virtual bool getNextCell(lluint& cellID,Real* coords,Real* sizes,std::vector<std::pair<lluint,uchar> >& nbrs) = 0;
   
   /** Request the value of a parameter from GridBuilder. This function is provided in 
    * cases where the user might want to know the values of grid-related parameters.
    * @param parameterName The name of the parameter.
    * @param value A string where the value of the parameter is to be written.
    * @return If true, GridBuilder returned a valid parameter value. If false, 
    * GridBuilder did not recognize the parameter.
    */
   virtual bool getParameter(const std::string& parameterName,std::string& value) =0;
   
   /** Query the total number of cells the GridBuilder creates. This 
    * function gets called before getNextCell. 
    * @param N_cells A variable where the total number of cells are to be written.
    * @return If true, GridBuilder wrote a valid value to parameter N_cells.
    */
   virtual bool getTotalNumberOfCells(lluint& N_cells) =0;
   
   /** Initialize the GridBuilder. If the GridBuilder, for example, reads the 
    * cells from a file, this is the correct place to open that file.
    * This function gets called before any other GridBuilder member functions 
    * are called, except for the constructor.
    * %Parameters and options read from the simulation input files, via calls to 
    * class Parameters member functions, should also be defined and read here. 
    * GridBuilder input file parameters must be defined in a section called "gridbuilder".
    * @return If true, GridBuilder initialized successfully.
    */
   virtual bool initialize() =0;
   
 protected:

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