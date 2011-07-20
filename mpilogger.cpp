#include <cstdlib>
#include <time.h>
#include "mpilogger.h"
#include <iostream>

using namespace std;

/** Default constructor for class MPILogger.
 */
MPILogger::MPILogger() {
   fileOpen = false;
   masterStream = NULL;
}

/** Destructor for class MPILogger. Destructor only calls MPILogger::close.
 */
MPILogger::~MPILogger() {
   close();
   if (masterStream != NULL) {
      masterStream->close();
      delete masterStream;
   }
}

/** Close a logfile which has been previously opened with 
 * MPILogger::open.
 * @return If true, the logfile was closed successfully.
 */
bool MPILogger::close() {
   if (fileOpen == false) return false;
   #ifndef NDEBUG
      bool success = mpiFile.close();
   #else 
      bool success = true;
      masterStream->close();
      delete masterStream;
      masterStream = NULL;
   #endif
   fileOpen = false;
   return success;
}

/** Write the MPILogger stream buffer into the output file. The stream buffer 
 * is emptied if the write succeeded. The rank of the writing MPI process and 
 * current time and date, as given by ctime, are written before the user-defined 
 * message.
 * @return If true, the buffer was written successfully and buffer was emptied.
 */
bool MPILogger::flush() {
   if (fileOpen == false) return false;
   bool success = true;
   #ifdef NDEBUG
      if (mpiRank != masterRank) return true;
   #endif
   
   // Get the current time. Remove the line break from 
   // the string output given by ctime.
   const time_t rawTime = time(NULL);
   string strTime = ctime(&rawTime);
   
   const size_t lineBreak = strTime.find('\n');
   if (lineBreak != string::npos) {
      strTime.erase(lineBreak,1);
   }

   // Form a new output buffer containing the process number, date 
   // and the user-given message. Attempt to write it to file.
   stringstream tmp;
   tmp << "PROC #" << mpiRank << ' ';
   tmp << strTime << ' ' << endl;
   strTime = outStream.str();
   tmp << strTime << endl;
   #ifndef NDEBUG
      mpiFile << tmp.str();

       // If the write was successful, clear stream buffer:
      if (mpiFile.getCount<size_t>() != tmp.str().size()) return false;
      outStream.str(string(""));
      outStream.clear();
   #else
      (*masterStream) << tmp.str() << std::flush;
      outStream.str(string(""));
      outStream.clear();
   #endif
   return true;
}

/** Open a logfile for parallel writing within the given MPI communicator. 
 * All MPI processes within the given communicator must call this function simultaneously.
 * @param comm Communicator which is used for this logfile.
 * @param fname The name of the logfile.
 * @return If true, the file was opened successfully and MPILogger is ready for use.
 */
bool MPILogger::open(MPI_Comm comm,const int& MASTER_RANK,const std::string& fname,const bool& deleteFile) {
   // Store the MPI rank of this process
   MPI_Comm_rank(comm,&mpiRank);
   masterRank = MASTER_RANK;
   bool rvalue = true;
   
   // If NDEBUG has been defined, only master process writes log messages.
   #ifndef NDEBUG
      if (deleteFile == true) MPI_File_delete(const_cast<char*>(fname.c_str()),MPI_INFO_NULL);
      const int accessMode = (MPI_MODE_WRONLY | MPI_MODE_SEQUENTIAL | MPI_MODE_CREATE);
      if (mpiFile.open(comm,fname,MPI_INFO_NULL,accessMode,true) == false) rvalue = false;
      if (mpiFile.resetPosition() == false) rvalue = false;
      if (rvalue == true) fileOpen = true;
   #else
      if (mpiRank != MASTER_RANK) return rvalue;
      masterStream = new fstream;
      masterStream->open(fname.c_str(), fstream::out);
      if (masterStream->good() == false) rvalue = false;
      fileOpen = true;
   #endif
   return rvalue;
}

bool MPILogger::print(const std::string& s) {
   #ifndef NDEBUG
      mpiFile << s;
      if (mpiFile.getCount<size_t>() != s.size()) return false;
   #else
      if (mpiRank != masterRank) return true;
      (*masterStream) << s;
   #endif
   return true;
}

// *********************************
// ****** STREAM MANIPULATORS ******
// *********************************

/** Function which allows one to insert C++-style stream manipulators, 
 * defined inside MPILogger, to be used in the input stream. This function gets 
 * mainly called when manipulator write is inserted to the stream.
 * @param pf A pointer to a stream manipulator.
 * @return A reference to MPILogger.
 */
MPILogger& MPILogger::operator<<(MPILogger& (*pf)(MPILogger&)) {
   #ifndef NDEBUG
      return (*pf)(*this);
   #else
      if (mpiRank != masterRank) return *this;
      return (*pf)(*this);
   #endif
}

/** Function which allows one to use C++ stream manipulators such as endl or 
 * showpos, with MPILogger. This is just a wrapper function which passes the 
 * manipulators to internal stream buffer as-is.
 * @param pf Function pointer to C++ stream manipulator.
 * @return Reference to MPILogger.
 */
MPILogger& MPILogger::operator<<(std::ostream& (*pf)(std::ostream& )) {
   #ifndef NDEBUG
      (*pf)(outStream);
      return *this;
   #else
      if (mpiRank != masterRank) return *this;
      (*pf)(outStream);
      return *this;
   #endif
}

/** C++-style stream manipulator which tells MPILogger to write the contents 
 * of the stream buffer into the file. The stream buffer is cleared after a 
 * successful write. You need to insert "write" into the 
 * input stream whenever you want to write to the logfile, endl or flush will 
 * not do that.
 * @param logger Reference to MPILogger.
 * @return Reference to MPILogger.
 */
MPILogger& write(MPILogger& logger) {
   if (logger.flush() == false) {
      #ifndef NDEBUG
         cerr << "Logger failed to write!" << endl;
      #endif
   }
   return logger;
}
