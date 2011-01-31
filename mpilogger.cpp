#include <cstdlib>
#include <iostream>
#include <time.h>
#include "mpilogger.h"

using namespace std;

/** Default constructor for class MPILogger.
 */
MPILogger::MPILogger() {
   fileOpen = false;
}

/** Destructor for class MPILogger. Destructor only calls MPILogger::close.
 */
MPILogger::~MPILogger() {
   close();
}

/** Close a logfile which has been previously opened with 
 * MPILogger::open.
 * @return If true, the logfile was closed successfully.
 */
bool MPILogger::close() {
   if (fileOpen == false) return false;
   bool success = mpiFile.close();
   if (success == true) fileOpen = false;
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
   tmp << strTime << ' ' << std::endl;
   strTime = outStream.str();
   tmp << strTime;
   
   mpiFile << tmp.str();

   // If the write was successful, clear stream buffer:
   if (mpiFile.getCount<size_t>() != tmp.str().size()) return false;
   outStream.str(string(""));
   outStream.clear();
   return true;
}

/** Open a logfile for parallel writing within the given MPI communicator. 
 * All MPI processes within the given communicator must call this function simultaneously.
 * @param comm Communicator which is used for this logfile.
 * @param fname The name of the logfile.
 * @return If true, the file was opened successfully and MPILogger is ready for use.
 */
bool MPILogger::open(MPI_Comm comm,const std::string& fname) {
   // Store the MPI rank of this process
   MPI_Comm_rank(comm,&mpiRank);
   
   const int accessMode = (MPI_MODE_WRONLY | MPI_MODE_SEQUENTIAL | MPI_MODE_CREATE);
   bool rvalue = true;
   if (mpiFile.open(comm,fname,MPI_INFO_NULL,accessMode,true) == false) rvalue = false;
   if (mpiFile.resetPosition() == false) rvalue = false;
   if (rvalue == true) fileOpen = true;
   return rvalue;
}

bool MPILogger::print(const std::string& s) {
   mpiFile << s;
   if (mpiFile.getCount<size_t>() != s.size()) return false;
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
   return (*pf)(*this);
}

/** Function which allows one to use C++ stream manipulators such as endl or 
 * showpos, with MPILogger. This is just a wrapper function which passes the 
 * manipulators to internal stream buffer as-is.
 * @param pf Function pointer to C++ stream manipulator.
 * @return Reference to MPILogger.
 */
MPILogger& MPILogger::operator<<(std::ostream& (*pf)(std::ostream& )) {
   (*pf)(outStream);
   return *this;
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
      cerr << "Logger failed to write!" << endl;
   }
   return logger;
}




