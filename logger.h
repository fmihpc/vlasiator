/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <sstream>
#include <mpi.h>

/** A class for writing log messages in parallel. Logger functions in 
 * the same way as standard C++ input streams, i.e. you insert values by 
 * using the stream input operator <<.
 * 
 * There is one major difference between Logger and standard C++ input 
 * streams: the contents of Logger are not flushed/written into the output file 
 * when endl or flush are inserted. Instead, one has to insert manipulator write,
 * defined in logger.h, into the stream. 
 *
 * This is because Logger actually 
 * uses a stringstream buffer internally, and C++ standard manipulators are just 
 * passed on to the stringstream as-is. It is not possible to detect when endl 
 * or flush are inserted, because they have the same function signature as other 
 * manipulators (such as showpos). In principle it is possible to define new manipulators 
 * endl and flush here in the same manner as write, but then wrapper functions would 
 * need to be written for all standard C++ manipulators because compilers think that
 * std::endl and user-defined endl as ambiguous, i.e. compilers cannot decide which 
 * function should be linked.
 */
class Logger {
public:
   Logger();
   ~Logger();
   
   void clear() {outStream.clear();}
   bool close();
   bool flush(bool verbose);
   std::stringstream& getStream() {return outStream;}
   bool open(MPI_Comm comm,const int& MASTER_RANK,const std::string& fname,const bool& deleteFile=true);
   bool print(const std::string& s);
   std::string str() {return outStream.str();}
   
   // ****************************************
   // ****** STREAM INSERTION OPERATORS ******
   // ****************************************
   
   template<typename T> inline Logger& operator<<(const T& value);
   Logger& operator<<(Logger& (*pf)(Logger&));
   Logger& operator<<(std::ostream& (*pf)(std::ostream& ));
   
private:
   bool fileOpen;                       /**< If true, the class has an open MPIFile.*/
   int mpiRank;                         /**< The rank of the process using Logger within a user-defined communicator.*/
   int masterRank;                      /**< MPI rank of the master process.*/
   std::stringstream outStream;         /**< Output buffer.*/
   std::fstream* masterStream;          /**< Output stream for master process only.*/
};

/** Stream insertion operator for inserting values to the input stream.
 * The given value is passed on to internal stream buffer as-is.
 * @param value The value to be inserted into the stream.
 * @return A reference to Logger.
 */
template<typename T>
Logger& Logger::operator<<(const T& value) {
   outStream << value;
   return *this;
}

Logger& write(Logger& logger);
Logger& writeNVerb(Logger& logger);

#endif
