#include <cstdlib>
#include <iostream>
#include <limits>
#include <mpi.h>

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

   /** Write the reduced data into the given byte array. The number of bytes written
    * must equal the value returned by DRO::DataReductionOperator::getOutputByteSize.
    * @param byteArray Pointer to the byte array.
    * @return If true, the bytes were appended successfully. The base class function returns false.
    * @see DRO::DataReductionOperator::getOutputByteSize
    */
   bool DataReductionOperator::appendReducedData(unsigned char* const byteArray) {
      cerr << "ERROR: DataReductionOperator::appendReducedData called instead of derived class function!" << endl;
      return false;
   }

   /** Get the byte size of one element in the array of reduced data calculated by 
    * this DRO::DataReductionOperator. For example, if this DRO::DataReductionOperator 
    * produces a vector of floats, then the size of an entry is sizeof(float). If this 
    * DRO::DataReductionOperator calculates a scalar value, it is treated as an array 
    * of size one.
    * @return The byte size of an entry in reduced data array. The base class function 
    * returns zero.
    */
   unsigned char DataReductionOperator::getElementByteSize() const {
      cerr << "ERROR: DataReductionOperator::getElementByteSize called instead of derived class function!" << endl;
      return 0;
   }

   /** Get the name of the reduced data variable. The name is written to the disk as-is 
    * and is used in visualization.
    * @return The name of the data. The base class function returns an empty string.
    */
   std::string DataReductionOperator::getName() const {
      cerr << "ERROR: DataReductionOperator::getName called instead of derived class function!" << endl;
      return string("");
   }

   /** Get the size of the reduced data. The value returned by this function equals 
    * DRO::DataReductionOperator::getElementByteSize() times the number of elements
    * in the reduced data array, which is deduced from the value returned by 
    * DRO::DataReductionOperator::getVariableType.
    * @return The size of reduced data in bytes.
    */
   unsigned int DataReductionOperator::getOutputByteSize() const {
      unsigned int rvalue = this->getElementByteSize();
      // Multiply element byte size by the number of elements in output array:
      switch (this->getVariableType()) {
       case VlsVariable::NULLVARIABLE:
	 rvalue *= 0;
	 break;
       case VlsVariable::SCALAR:
	 rvalue *= 1;
	 break;
       case VlsVariable::VECTOR2:
	 rvalue *= 2;
	 break;
       case VlsVariable::VECTOR3:
	 rvalue *= 3;
	 break;
       case VlsVariable::TENSOR22:
	 rvalue *= 4;
	 break;
       case VlsVariable::TENSOR23:
	 rvalue *= 6;
	 break;
       case VlsVariable::TENSOR32:
	 rvalue *= 6;
	 break;
       case VlsVariable::TENSOR33:
	 rvalue *= 9;
	 break;
       default:
	 rvalue *= 0;
	 break;
      }
      return rvalue;
}
   
   /** Reduce the given data and store the value into internal variables.
    * @param block A velocity block in the current SpatialCell. The size of the block is 
    * WID times WID times WID.
    * @param blockParams %Parameters of the given velocity block, which can be accessed by 
    * using the values defined in namespace BlockParams. For example, the vx-coordinate of the 
    * lower left corner of the block is blockParams[%BlockParams::VXCRD].
    * @return If true, the data was reduced successfully. The base class function returns false.
    */
   bool DataReductionOperator::reduceData(const Real* const block,const Real* const blockParams) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function!" << endl;
      return false;
   }
   
   /** Get the type of variable returned by this DRO::DataReductionOperator. This 
    * function should return one of the values defined in namespace VlsHeader. 
    * For example, if the intent is to write out a scalar field, then VlsHeader::SCALAR 
    * should be returned.
    * @return The type ID of the reduced data variable written out by this operator.
    */
   unsigned char DataReductionOperator::getVariableType() const {
      cerr << "ERROR: DataReductionOperator::getVariableType called instead of derived class function!" << endl;
      return VlsVariable::NULLVARIABLE;
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
   
   bool VariableE::appendReducedData(unsigned char* const byteArray) {
      unsigned char* ptr;      
      ptr = reinterpret_cast<unsigned char*>(&Ex);
      byteArray[ 0] = ptr[0];
      byteArray[ 1] = ptr[1];
      byteArray[ 2] = ptr[2];
      byteArray[ 3] = ptr[3];
      ptr = reinterpret_cast<unsigned char*>(&Ey);
      byteArray[ 4] = ptr[0];
      byteArray[ 5] = ptr[1];
      byteArray[ 6] = ptr[2];
      byteArray[ 7] = ptr[3];
      ptr = reinterpret_cast<unsigned char*>(&Ez);
      byteArray[ 8] = ptr[0];
      byteArray[ 9] = ptr[1];
      byteArray[10] = ptr[2];
      byteArray[11] = ptr[3];
      return true;
   }
   
   unsigned char VariableE::getElementByteSize() const {return sizeof(Ex);}
   
   std::string VariableE::getName() const {return "E";}
   
   bool VariableE::reduceData(const Real* const avgs,const Real* const blockParams) {return true;}
   
   unsigned char VariableE::getVariableType() const {return VlsVariable::VECTOR3;}
   
   bool VariableE::setSpatialCell(const SpatialCell& cell) {
      Ex = cell.cpu_cellParams[CellParams::EX];
      Ey = cell.cpu_cellParams[CellParams::EY];
      Ez = cell.cpu_cellParams[CellParams::EZ];
      return true;
   }
   
   VariableB::VariableB(): DataReductionOperator() { }
   VariableB::~VariableB() { }
   
   bool VariableB::appendReducedData(unsigned char* const byteArray) {
      unsigned char* ptr;
      ptr = reinterpret_cast<unsigned char*>(&Bx);
      byteArray[ 0] = ptr[0];
      byteArray[ 1] = ptr[1];
      byteArray[ 2] = ptr[2];
      byteArray[ 3] = ptr[3];
      ptr = reinterpret_cast<unsigned char*>(&By);
      byteArray[ 4] = ptr[0];
      byteArray[ 5] = ptr[1];
      byteArray[ 6] = ptr[2];
      byteArray[ 7] = ptr[3];
      ptr = reinterpret_cast<unsigned char*>(&Bz);
      byteArray[ 8] = ptr[0];
      byteArray[ 9] = ptr[1];
      byteArray[10] = ptr[2];
      byteArray[11] = ptr[3];
      return true;
   }
   
   unsigned char VariableB::getElementByteSize() const {return sizeof(Bx);}
   
   std::string VariableB::getName() const {return "B";}
   
   bool VariableB::reduceData(const Real* const avgs,const Real* const blockParams) {return true;}
   
   unsigned char VariableB::getVariableType() const {return VlsVariable::VECTOR3;}
   
   bool VariableB::setSpatialCell(const SpatialCell& cell) {
      Bx = cell.cpu_cellParams[CellParams::BX];
      By = cell.cpu_cellParams[CellParams::BY];
      Bz = cell.cpu_cellParams[CellParams::BZ];
      return true;
   }
   
   VariableRho::VariableRho(): DataReductionOperator() { }
   VariableRho::~VariableRho() { }
   
   bool VariableRho::appendReducedData(unsigned char* const byteArray) {
      unsigned char* ptr = reinterpret_cast<unsigned char*>(&rho);
      for (unsigned int i=0; i<sizeof(rho); ++i) byteArray[i] = ptr[i];
      return true;
   }
   
   unsigned char VariableRho::getElementByteSize() const {return sizeof(rho);}
   
   std::string VariableRho::getName() const {return "rho";}
   
   bool VariableRho::reduceData(const Real* const avgs,const Real* const blockParams) {return true;}
   
   unsigned char VariableRho::getVariableType() const {return VlsVariable::SCALAR;}
   
   bool VariableRho::setSpatialCell(const SpatialCell& cell) {
      rho = cell.cpu_cellParams[CellParams::RHO];
      return true;
   }
   
   
   MPIrank::MPIrank(): DataReductionOperator() { }
   MPIrank::~MPIrank() { }
   
   bool MPIrank::appendReducedData(unsigned char* const byteArray) {
      unsigned char* ptr = reinterpret_cast<unsigned char*>(&rank);
      for (unsigned int i=0; i<sizeof(rank); ++i) byteArray[i] = ptr[i];
      return true;
   }
   
   unsigned char MPIrank::getElementByteSize() const {return sizeof(rank);}
   
   std::string MPIrank::getName() const {return "MPI_rank";}
   
   bool MPIrank::reduceData(const Real* const avgs,const Real* const blockParams) {return true;}
   
   unsigned char MPIrank::getVariableType() const {return VlsVariable::SCALAR;}
   
   bool MPIrank::setSpatialCell(const SpatialCell& cell) {
      int intRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&intRank);
      rank = 1.0*intRank;
      return true;
   }

   VariableRhoV::VariableRhoV(): DataReductionOperator() { }
   VariableRhoV::~VariableRhoV() { }
   
   bool VariableRhoV::appendReducedData(unsigned char* const byteArray) {
      unsigned char* ptr;
      ptr = reinterpret_cast<unsigned char*>(&rhovx);
      byteArray[ 0] = ptr[0];
      byteArray[ 1] = ptr[1];
      byteArray[ 2] = ptr[2];
      byteArray[ 3] = ptr[3];
      ptr = reinterpret_cast<unsigned char*>(&rhovy);
      byteArray[ 4] = ptr[0];
      byteArray[ 5] = ptr[1];
      byteArray[ 6] = ptr[2];
      byteArray[ 7] = ptr[3];
      ptr = reinterpret_cast<unsigned char*>(&rhovz);
      byteArray[ 8] = ptr[0];
      byteArray[ 9] = ptr[1];
      byteArray[10] = ptr[2];
      byteArray[11] = ptr[3];
      return true;
   }
   
   unsigned char VariableRhoV::getElementByteSize() const {return sizeof(rhovx);}
   
   std::string VariableRhoV::getName() const {return "rho_v";}
   
   bool VariableRhoV::reduceData(const Real* const avgs,const Real* const blockParams) {return true;}
   
   unsigned char VariableRhoV::getVariableType() const {return VlsVariable::VECTOR3;}
   
   bool VariableRhoV::setSpatialCell(const SpatialCell& cell) {
      rhovx = cell.cpu_cellParams[CellParams::RHOVX];
      rhovy = cell.cpu_cellParams[CellParams::RHOVY];
      rhovz = cell.cpu_cellParams[CellParams::RHOVZ];
      return true;
   }
   
} // namespace DRO
