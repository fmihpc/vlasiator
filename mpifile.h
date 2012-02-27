/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

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

#ifndef MPIWRITER_H
#define MPIWRITER_H

#include <mpi.h>
#include <string>


/** A small class for writing data as byte arrays into a file in parallel.
 * MPI must have been initialized prior to using this class.
 */
class MPIFile {
 public:
   MPIFile();
   ~MPIFile();
   
   bool open(MPI_Comm comm,const std::string& fname,MPI_Info file,const int& accessMode,const bool& deleteFile=false);
   bool close();
   /** Get the MPI status struct related to the previous MPI command.
    * See MPI specification for more details.
    * @return MPI status.
    */
   MPI_Status getStatus() const {return MPIstatus;}
   bool resetPosition();
   
   template<typename T> MPIFile& operator<<(const T& value);
   template<typename T> MPIFile& operator>>(T& value);
   template<typename T> T getCount();
   template<typename T,typename INT> bool read(const INT& count,T* const array);
   template<typename T,typename INT> bool write(const INT& count,const T* const array);
   
 private:
   MPI_File fileptr;     /**< Pointer to file opened with MPI. See MPI specification for more details.*/
   MPI_Info MPIinfo;     /**< Variable for holding MPI file info. See MPI specification for more details.*/
   MPI_Status MPIstatus; /**< Variable for holding MPI status messages. See MPI specification for more details.*/
   void* ptr;            /**< A pointer used when writing datatypes to files as byte arrays.*/
};

// ********************************************************
// ****** DEFINITIONS FOR MPIFILE TEMPLATE FUNCTIONS ******
// ********************************************************

/** Stream input operator for writing a primitive datatype into an MPI file. 
 * This function should work for all datatypes for which sizeof gives an accurate 
 * and portable size, i.e. do not use this function for writing structs or classes.
 * In such cases use MPIFile::write. It is mandatory to call MPIFile::resetPosition 
 * before attempting to access the file.
 * @param value The data to write.
 * @return A reference to MPIFile.
 * @see MPIFile::read.
 * @see MPIFile::resetPosition
 */
template<typename T> inline MPIFile& MPIFile::operator<<(const T& value) {
   ptr = static_cast<void*>(const_cast<T*>(&value));
   MPI_File_write_shared(fileptr,ptr,sizeof(T),MPI_BYTE,&MPIstatus);
   return *this;
}

template<> inline MPIFile& MPIFile::operator<<(const std::string& value) {
   ptr = static_cast<void*>(const_cast<char*>(&(value[0])));
   MPI_File_write_shared(fileptr,ptr,value.size(),MPI_BYTE,&MPIstatus);
   return *this;
}

/** Stream output operator for reading a primitive datatype from an MPI file.
 * This function should work for all datatypes for which sizeof gives an 
 * accurate and portable size, i.e. do not use this function for reading
 * structs or classes. Use MPIFIle::read instead. It is mandatory to call 
 * MPIFile::resetPosition before attampting to access the file.
 * @param value A variable in which the read value is written.
 * @return A reference to MPIFIle.
 * @see MPIFile::read
 * @see MPIFile::resetPosition
 */
template<typename T> inline MPIFile& MPIFile::operator>>(T& value) {
   ptr = static_cast<void*>(&value);
   MPI_File_read_shared(fileptr,ptr,sizeof(T),MPI_BYTE,&MPIstatus);
   return *this;
}

/** Get the number of bytes that were written to file with the last 
 * write command.
 * @return The number of bytes successfully written to file.
 */
template<typename T> inline T MPIFile::getCount() {
   int count;
   MPI_Get_count(&MPIstatus,MPI_BYTE,&count);
   return static_cast<T>(count);
}

/** A function for reading an array of primitive datatypes from an MPI file. This function 
 * should be used to read non-trivial datatypes, such as structs and class data, 
 * by first reading a byte array containing the data and then typecasting it into correct 
 * variables. It is mandatory to call MPIFile::resetPosition before attempting to access 
 * the file.
 * @param count Number of elements to read.
 * @param array An array whose contents will be overwritten with the read data. This array 
 * must be allocated prior to calling this function.
 * @return If true, the data was read successfully.
 * @see MPIFile::resetPosition
 */
template<typename T,typename INT> inline bool MPIFile::read(const INT& count,T* const array) {
   if (MPI_File_read_shared(fileptr,array,sizeof(T)*count,MPI_BYTE,&MPIstatus) == MPI_SUCCESS) return true;
   return false;
}

/** A function for writing an array of primitive datatypes into an MPI file. This function
 * should be used to write non-trivial datatypes, such as structs and class data,
 * by first converting the data into a byte array. It is mandatory to call MPIFile::resetPosition 
 * before attempting to access the file.
 * @param count Number of elements in array.
 * @param array An array of primitive datatypes to write.
 * @return If true, the data was written successfully.
 * @see MPIFile::resetPosition
 */
template<typename T,typename INT> inline bool MPIFile::write(const INT& count,const T* const array) {
   ptr = static_cast<void*>(const_cast<T*>(array));
   bool rvalue = true;
   //if (MPI_File_write_shared(fileptr,ptr,sizeof(T)*count,MPI_BYTE,&MPIstatus) != MPI_SUCCESS) rvalue = false;   
   MPI_Request mpiRequest;
   if (MPI_File_iwrite_shared(fileptr,ptr,sizeof(T)*count,MPI_BYTE,&mpiRequest) != MPI_SUCCESS) rvalue = false;
   MPI_Wait(&mpiRequest,&MPIstatus);
   
   return rvalue;
}

#endif
