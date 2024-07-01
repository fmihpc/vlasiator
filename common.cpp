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
#include <mpi.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "common.h"

/*! \brief A function to stop the simulation if the boolean condition is true.
 * Raises a flag which gets MPI_Reduced and initiates bailout.
 * [bailout] parameters set details.
 * \param condition if true, bailout is initiated
 * \param message information message printed to cerr
 * \param file code file where bailout was called (use __FILE__ in the call)
 * \param line code line where bailout was called (use __LINE__ in the call)
 */
void bailout(
   const bool condition,
   const std::string& message,
   const char * const file,
   const int line
) {
   #pragma omp critical
   {
      if (condition && (globalflags::bailingOut == 0)) {
         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
         std::cerr << "Process " << myRank << " bailing out";
         if((strcmp(file, "") != 0) && (line != 0)) {
            std::cerr << " at " << file << ":" << line;
         }
         std::cerr << ".";
         if(strcmp(message.c_str(), "") != 0) {
            std::cerr << " " << message;
         }
         std::cerr << std::endl;
         globalflags::bailingOut = 1;
      }
   }
}

/*! \brief A function to stop the simulation if the boolean condition is true.
 * Raises a flag which gets MPI_Reduced and initiates bailout.
 * [bailout] parameters set details.
 * \param condition if true, bailout is initiated
 * \param file code file where bailout was called (use __FILE__ in the call)
 * \param line code line where bailout was called (use __LINE__ in the call)
 */
void bailout(
   const bool condition,
   const char * const file,
   const int line
) {
   bailout(condition, "", file, line);
}

/*! \brief A function to stop the simulation if the boolean condition is true.
 * Raises a flag which gets MPI_Reduced and initiates bailout.
 * [bailout] parameters set details.
 * \param condition if true, bailout is initiated
 * \param message information message printed to cerr
 */
void bailout(
   const bool condition,
   const std::string& message
) {
   bailout(condition, message, "", 0);
}

/*! Helper function for error handling. err_type default to 0.*/
[[ noreturn ]] void abort_mpi(const std::string str, const int err_type) {
   // Single string so output isn't mangled by multiple processes
   std::cerr << (err_type ? std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + str : str) + "\n";
   MPI_Abort(MPI_COMM_WORLD, 1);

   // Dummy abort to convince compiler function doesn't return
   // TODO replace with std::unreachable once we switch to C++23
   abort();
}
