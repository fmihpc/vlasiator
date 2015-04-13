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
   const std::string message,
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
   const std::string message
) {
   bailout(condition, message, "", 0);
}