#ifndef MPILOGGER_H
#define MPILOGGER_H

#include "mpifile.h"
#include <sstream>

/** A class for writing log messages in parallel. MPILogger functions in 
 * the same way as standard C++ input streams, i.e. you insert values by 
 * using the stream input operator <<.
 * 
 * There is one major difference between MPILogger and standard C++ input 
 * streams: the contents of MPILogger are not flushed/written into the output file 
 * when endl or flush are inserted. Instead, one has to insert manipulator write,
 * defined in mpilogger.h, into the stream. 
 *
 * This is because MPILogger actually 
 * uses a stringstream buffer internally, and C++ standard manipulators are just 
 * passed on to the stringstream as-is. It is not possible to detect when endl 
 * or flush are inserted, because they have the same function signature as other 
 * manipulators (such as showpos). In principle it is possible to define new manipulators 
 * endl and flush here in the same manner as write, but then wrapper functions would 
 * need to be written for all standard C++ manipulators because compilers think that
 * std::endl and user-defined endl as ambiguous, i.e. compilers cannot decide which 
 * function should be linked.
 */
class MPILogger {
 public:
   MPILogger();
   ~MPILogger();
   
   void clear() {outStream.clear();}
   bool close();
   bool flush();
   std::stringstream& getStream() {return outStream;}
   bool open(MPI_Comm comm,const int& MASTER_RANK,const std::string& fname,const bool& deleteFile=true);
   bool print(const std::string& s);
   std::string str() {return outStream.str();}

   // ****************************************
   // ****** STREAM INSERTION OPERATORS ******
   // ****************************************
   
   template<typename T> inline MPILogger& operator<<(const T& value);
   MPILogger& operator<<(MPILogger& (*pf)(MPILogger&));
   MPILogger& operator<<(std::ostream& (*pf)(std::ostream& ));
   
 private:
   bool fileOpen;                       /**< If true, the class has an open MPIFile.*/
   MPIFile mpiFile;                     /**< MPIFile which is used for parallel I/O.*/
   int mpiRank;                         /**< The rank of the process using MPILogger within a user-defined communicator.*/
   int masterRank;                      /**< MPI rank of the master process.*/
   std::stringstream outStream;         /**< Output buffer.*/
};

/** Stream insertion operator for inserting values to the input stream.
 * The given value is passed on to internal stream buffer as-is.
 * @param value The value to be inserted into the stream.
 * @return A reference to MPILogger.
 */
template<typename T>
MPILogger& MPILogger::operator<<(const T& value) {
   outStream << value;
   return *this;
}

MPILogger& write(MPILogger& logger);

#endif
