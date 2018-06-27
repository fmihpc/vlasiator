/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <time.h>
#include "logger.h"
#include <iostream>

using namespace std;

/** Default constructor for class Logger.
 */
Logger::Logger() {
   fileOpen = false;
   masterStream = NULL;
}

/** Destructor for class Logger. Destructor only calls Logger::close.
 */
Logger::~Logger() {
   close();
   if (masterStream != NULL) {
      masterStream->close();
      delete masterStream;
   }
}

/** Close a logfile which has been previously opened with 
 * Logger::open.
 * @return If true, the logfile was closed successfully.
 */
bool Logger::close() {
   if (fileOpen == false) return false;
   bool success = true;
   masterStream->close();
   delete masterStream;
   masterStream = NULL;
   fileOpen = false;
   return success;
}

/** Write the Logger stream buffer into the output file. The stream buffer 
 * is emptied if the write succeeded. The rank of the writing MPI process and 
 * current time and date, as given by ctime, are written before the user-defined 
 * message.
 * @return If true, the buffer was written successfully and buffer was emptied.
 */
bool Logger::flush(bool verbose) {
   if (fileOpen == false) return false;
   bool success = true;
   if (mpiRank != masterRank) return true;
   stringstream tmp;
   string strTime;
   if (verbose == true)
   {
      // Get the current time. Remove the line break from 
      // the string output given by ctime.
      const time_t rawTime = time(NULL);
      strTime = ctime(&rawTime);
      
      const size_t lineBreak = strTime.find('\n');
      if (lineBreak != string::npos) {
         strTime.erase(lineBreak,1);
      }
      
      // Form a new output buffer containing the process number, date 
      // and the user-given message. Attempt to write it to file.
      
      tmp << strTime << ' ' << endl;
      strTime = outStream.str();
      tmp << strTime << endl;
      
   } else {
      strTime = outStream.str();
      tmp << strTime;
   }
   (*masterStream) << tmp.str() << std::flush;
   outStream.str(string(""));
   outStream.clear();
   return success;
}

/*!
   \brief Open a logfile

   This class can be used to write out logfiles. Only master rank
   writes out anythin, all other processes will not do anything when
   they attempt to use the logger, but it is still perfectly legal for
   them to call it.
   
      
   \param[in] comm  MPI communicator for processes that are able to call this logger
   \param[in] MASTERRANK The rank that does all communication.
   \param[in] fname The name of the logfile.
   \param[in] append Is the file appended (true), or overwritten (false)
   \return If true, the file was opened successfully and Logger is ready for use.
 */
bool Logger::open(MPI_Comm comm,const int& MASTERRANK,const std::string& fname,const bool& append) {
   // Store the MPI rank of this process
   MPI_Comm_rank(comm,&mpiRank);
   masterRank = MASTERRANK;
   bool rvalue = true;
   
   if (mpiRank != MASTERRANK) return rvalue;
   masterStream = new fstream;
   if(append)
      masterStream->open(fname.c_str(), fstream::out|fstream::app);
   else
      masterStream->open(fname.c_str(), fstream::out);

   if (masterStream->good() == false) rvalue = false;
   fileOpen = true;
   return rvalue;
}

bool Logger::print(const std::string& s) {
   if (mpiRank != masterRank) return true;
   (*masterStream) << s;
   return true;
}

// *********************************
// ****** STREAM MANIPULATORS ******
// *********************************

/** Function which allows one to insert C++-style stream manipulators, 
 * defined inside Logger, to be used in the input stream. This function gets 
 * mainly called when manipulator write is inserted to the stream.
 * @param pf A pointer to a stream manipulator.
 * @return A reference to Logger.
 */
Logger& Logger::operator<<(Logger& (*pf)(Logger&)) {
   if (mpiRank != masterRank) return *this;
   return (*pf)(*this);
}

/** Function which allows one to use C++ stream manipulators such as endl or 
 * showpos, with Logger. This is just a wrapper function which passes the 
 * manipulators to internal stream buffer as-is.
 * @param pf Function pointer to C++ stream manipulator.
 * @return Reference to Logger.
 */
Logger& Logger::operator<<(std::ostream& (*pf)(std::ostream& )) {
   if (mpiRank != masterRank) return *this;
   (*pf)(outStream);
   return *this;
}

/** C++-style stream manipulator which tells Logger to write the contents 
 * of the stream buffer into the file. The stream buffer is cleared after a 
 * successful write. You need to insert "write" into the 
 * input stream whenever you want to write to the logfile, endl or flush will 
 * not do that.
 * Argument is true because it writes verbose output.
 * @param logger Reference to Logger.
 * @return Reference to Logger.
 */
Logger& writeVerbose(Logger& logger) {
   if (logger.flush(true) == false) {
      cerr << "Logger failed to write!" << endl;
   }
   return logger;
}

/** C++-style stream manipulator which tells Logger to write the contents 
 * of the stream buffer into the file. The stream buffer is cleared after a 
 * successful write. You need to insert "write" into the 
 * input stream whenever you want to write to the logfile, endl or flush will 
 * not do that.
 * Argument is false because it writes non-verbose output.
 * @param logger Reference to Logger.
 * @return Reference to Logger.
 */
Logger& write(Logger& logger) {
   if (logger.flush(false) == false) {
      cerr << "Logger failed to write!" << endl;
   }
   return logger;
}
