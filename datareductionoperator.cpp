#include <cstdlib>
#include <iostream>
#include <limits>
#include <mpi.h>

#include "datareductionoperator.h"
#include "vlswriter.h"

using namespace std;

namespace DRO {

   // ************************************************************
   // ***** DEFINITIONS FOR DATAREDUCTIONOPERATOR BASE CLASS *****
   // ************************************************************
   
   /** DataReductionOperator base class constructor. The constructor is empty, 
    * i.e. it does not initialize anything.
    */
   DataReductionOperator::DataReductionOperator() { }

   /** DataReductionOperator base class virtual destructor. The destructor is empty, 
    * i.e. it does not deallocate anything.
    */
   DataReductionOperator::~DataReductionOperator() { }

   /** Append the reduced data into the given byte array. The number of bytes appended 
    * must equal the value returned by getOutputSizeBytes.
    * @param byteArray Pointer to the byte array.
    * @return If true, the bytes were appended successfully. The base class function returns false.
    * @see getOutputSizeBytes
    */
   bool DataReductionOperator::appendReducedData(unsigned char* const byteArray) {
      cerr << "ERROR: DataReductionOperator::appendReducedData called instead of derived class function!" << endl;
      return false;
   }

   /** Get the byte size of one element in the array containing reduced data. 
    * Typically one wants to return sizeof(float) or sizeof(double).
    * @param Byte size of array element.
    */
   unsigned char DataReductionOperator::getElementByteSize() const {
      cerr << "ERROR: DataReductionOperator::getElementByteSize called instead of derived class function!" << endl;
      return 0;
   }

   /** Get the name of the reduced data. The name is written to the disk as-is and used in 
    * visualization.
    * @return The name of the data. The base class function returns an empty string.
    */
   std::string DataReductionOperator::getName() const {
      cerr << "ERROR: DataReductionOperator::getName called instead of derived class function!" << endl;
      return string("");
   }

   /** Get the size of the reduced data. The value returned by this function equals 
    * getElementByteSize() times the number of returned variables.
    * @return The size of reduced data in bytes.
    */
   unsigned int DataReductionOperator::getOutputByteSize() const {
      unsigned int rvalue = this->getElementByteSize();
      // Multiply element byte size by the number of elements in output array:
      switch (this->getVariableType()) {
       case VlsHeader::NULLVARIABLE:
	 rvalue *= 0;
	 break;
       case VlsHeader::SCALAR:
	 rvalue *= 1;
	 break;
       case VlsHeader::VECTOR2:
	 rvalue *= 2;
	 break;
       case VlsHeader::VECTOR3:
	 rvalue *= 3;
	 break;
       case VlsHeader::TENSOR22:
	 rvalue *= 4;
	 break;
       case VlsHeader::TENSOR23:
	 rvalue *= 6;
	 break;
       case VlsHeader::TENSOR32:
	 rvalue *= 6;
	 break;
       case VlsHeader::TENSOR33:
	 rvalue *= 9;
	 break;
       default:
	 rvalue *= 0;
	 break;
      }
      return rvalue;
}
   
   /** Reduce the given data and store the value into internal variables.
    * @return If true, the data was reduced successfully. The base class function returns false.
    */
   bool DataReductionOperator::reduceData(const Real* const avgs,const Real* const blockParams) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function!" << endl;
      return false;
   }
   
   unsigned char DataReductionOperator::getVariableType() const {
      cerr << "ERROR: DataReductionOperator::getVariableType called instead of derived class function!" << endl;
      return VlsHeader::NULLVARIABLE;
   }
   
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
   
   unsigned char VariableE::getVariableType() const {return VlsHeader::VECTOR3;}
   
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
   
   unsigned char VariableB::getVariableType() const {return VlsHeader::VECTOR3;}
   
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
   
   unsigned char VariableRho::getVariableType() const {return VlsHeader::SCALAR;}
   
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
   
   unsigned char MPIrank::getVariableType() const {return VlsHeader::SCALAR;}
   
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
   
   unsigned char VariableRhoV::getVariableType() const {return VlsHeader::VECTOR3;}
   
   bool VariableRhoV::setSpatialCell(const SpatialCell& cell) {
      rhovx = cell.cpu_cellParams[CellParams::RHOVX];
      rhovy = cell.cpu_cellParams[CellParams::RHOVY];
      rhovz = cell.cpu_cellParams[CellParams::RHOVZ];
      return true;
   }
   
} // namespace DRO
