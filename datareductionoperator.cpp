#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <limits>

#include "datareductionoperator.h"
#include "vlscommon.h"

using namespace std;

namespace DRO {

   // ************************************************************
   // ***** DEFINITIONS FOR DATAREDUCTIONOPERATOR BASE CLASS *****
   // ************************************************************
   
   /** DataReductionOperator base class constructor. The constructor is empty.*/
   DataReductionOperator::DataReductionOperator() { }

   /** DataReductionOperator base class virtual destructor. The destructor is empty.*/
   DataReductionOperator::~DataReductionOperator() { }

   /** Get info on the data the DataReductionOperator writes on disk. A DataReductionOperator writes 
    * an array on disk. Each element of the array is a vector with n elements. Finally, each 
    * vector element has a byte size, as given by the sizeof function.
    * @param dataType Basic datatype, must be int, uint, or float.
    * @param dataSize Byte size of written datatype, for example double-precision floating points
    * have byte size of sizeof(double).
    * @param vectorSize How many elements are in the vector returned by the DataReductionOperator.
    * @return If true, DataReductionOperator returned sensible values.
    */
   bool DataReductionOperator::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      cerr << "ERROR: DataReductionOperator::getDataVectorInfo called insted of derived class function!" << endl;
      return false;
   }
   
   /** Get the name of the reduced data variable. The name is written to the disk as-is 
    * and is used in visualization.
    * @return The name of the data. The base class function returns an empty string.
    */
   std::string DataReductionOperator::getName() const {
      cerr << "ERROR: DataReductionOperator::getName called instead of derived class function!" << endl;
      return string("");
   }
   
   /** Reduce the data and write the data vector to the given buffer.
    * @param N_blocks Number of velocity blocks in array avgs.
    * @param avgs Array containing distribution function values for each velocity block.
    * @param blockParams Array containing the parameters of each velocity block.
    * @param buffer Buffer in which the reduced data is written.
    * @return If true, DataReductionOperator reduced data successfully.
    */
   bool DataReductionOperator::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function!" << endl;
      return false;
   }
   
   /** Set the SpatialCell whose data is going to be reduced by subsequent calls to 
    * DRO::DataReductionOperator::reduceData. This function is provided so that 
    * variables stored per SpatialCell can be accessed.
    * 
    * Spatial cell variables are stored in array SpatialCell::cpu_cellParams. 
    * The contents of array elements are stored in namespace CellParams. For example, 
    * cell.cpu_cellParams[%CellParams::EX] contains the electric field.
    * @param cell The SpatialCell whose data is to be reduced next.
    * @return If true, the SpatialCell was set correctly.
    */
   bool DataReductionOperator::setSpatialCell(const SpatialCell& cell) {
      cerr << "ERROR: DataReductionOperator::setSpatialCell called instead of derived class function!" << endl;
      return false;
   }
   
   VariableE::VariableE(): DataReductionOperator() { }
   VariableE::~VariableE() { }
   
   bool VariableE::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize = 4;
      vectorSize = 3;
      return true;
   }
   
   std::string VariableE::getName() const {return "E";}
   
   bool VariableE::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(E);
      for (int i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableE::setSpatialCell(const SpatialCell& cell) {
      E = &(cell.cpu_cellParams[CellParams::EX]);
      Ex = cell.cpu_cellParams[CellParams::EX];
      Ey = cell.cpu_cellParams[CellParams::EY];
      Ez = cell.cpu_cellParams[CellParams::EZ];
      return true;
   }
   
   VariableB::VariableB(): DataReductionOperator() { }
   VariableB::~VariableB() { }

   bool VariableB::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize = 4;
      vectorSize = 3;
      return true;
   }
   
   std::string VariableB::getName() const {return "B";}
   
   bool VariableB::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(B);
      for (int i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableB::setSpatialCell(const SpatialCell& cell) {
      B  = &(cell.cpu_cellParams[CellParams::BX]);
      Bx = cell.cpu_cellParams[CellParams::BX];
      By = cell.cpu_cellParams[CellParams::BY];
      Bz = cell.cpu_cellParams[CellParams::BZ];
      return true;
   }
   
   VariableRho::VariableRho(): DataReductionOperator() { }
   VariableRho::~VariableRho() { }
   
   std::string VariableRho::getName() const {return "rho";}
   
   bool VariableRho::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize = 4;
      vectorSize = 1;
      return true;
   }
   
   bool VariableRho::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&rho);
      for (int i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRho::setSpatialCell(const SpatialCell& cell) {
      rho = cell.cpu_cellParams[CellParams::RHO];
      return true;
   }
   
   
   MPIrank::MPIrank(): DataReductionOperator() { }
   MPIrank::~MPIrank() { }
   
   bool MPIrank::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = 4;
      vectorSize = 1;
      return true;
   }
   
   std::string MPIrank::getName() const {return "MPI_rank";}
   
   bool MPIrank::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&mpiRank);
      for (int i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MPIrank::setSpatialCell(const SpatialCell& cell) {
      int intRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&intRank);
      rank = 1.0*intRank;
      mpiRank = intRank;
      return true;
   }

   VariableRhoV::VariableRhoV(): DataReductionOperator() { }
   VariableRhoV::~VariableRhoV() { }
   
   bool VariableRhoV::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize = 4;
      vectorSize = 3;
      return true;
   }
   
   std::string VariableRhoV::getName() const {return "rho_v";}
   
   bool VariableRhoV::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(rhov);
      for (int i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoV::setSpatialCell(const SpatialCell& cell) {
      rhov  = &(cell.cpu_cellParams[CellParams::RHOVX]);
      rhovx = cell.cpu_cellParams[CellParams::RHOVX];
      rhovy = cell.cpu_cellParams[CellParams::RHOVY];
      rhovz = cell.cpu_cellParams[CellParams::RHOVZ];
      return true;
   }
   
} // namespace DRO
