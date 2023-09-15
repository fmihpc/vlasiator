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

#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <stdint.h>
#include <string>
#include <typeinfo>
#include <vector>

#include "common.h"
#include "version.h"

class Readparameters {
public:
   Readparameters(int cmdargc, char* cmdargv[]);
   ~Readparameters();

   /** Add a new input parameter.
    * Note that parse must be called in order for the input file(s) to be re-read.
    * Only called by the root process.
    * @param name The name of the parameter, as given in the input file(s).
    * @param desc Description for the parameter.
    * @param defValue Default value for variable.
    */
   static void add(const std::string& name, const std::string& desc, const std::string& defValue) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == MASTER_RANK) {
         options[name] = "";
         isOptionParsed[name] = false;
         descriptions->add_options()(
             name.c_str(), boost::program_options::value<std::string>(&(options[name]))->default_value(defValue),
             desc.c_str());
      }
   }

   template <typename T> static void add(const std::string& name, const std::string& desc, const T& defValue) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == MASTER_RANK) {
         std::stringstream ss;

         if constexpr (std::is_floating_point_v<T>) {
            ss << std::setprecision(std::numeric_limits<double>::digits10 + 1) << defValue;
         } else {
            ss << defValue;
         }
         options[name] = "";
         isOptionParsed[name] = false;
         descriptions->add_options()(
             name.c_str(), boost::program_options::value<std::string>(&(options[name]))->default_value(ss.str()),
             desc.c_str());
      }
   }

   /** Get the value of the given parameter added with add().
    * This may be called after having called Parse, and it may be called by any process, in any order.
    * Aborts if given parameter was not found (a parameter passed to get() wasn't add()ed, defaults are ok).
    * @param name The name of the parameter.
    * @param value A variable where the value of the parameter is written.
    */
   static void get(const std::string& name, std::string& value) {
      if (options.find(name) != options.end()) { // check if it exists
         value = options[name];
      } else {
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         if (rank == MASTER_RANK) {
            std::cerr << __FILE__ << ":" << __LINE__ << name + " not declared using the add() function!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
   }

   static void get(const std::string& name, std::vector<std::string>& value) {
      if (vectorOptions.find(name) != vectorOptions.end()) { // check if it exists
         value = vectorOptions[name];
      } else {
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         if (rank == MASTER_RANK) {
            std::cerr << __FILE__ << ":" << __LINE__ << name + " not declared using the add() function!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
   }

   template <typename T> static void get(const std::string& name, T& value) {
      std::string sval;
      get(name, sval);

      try {
         value = boost::lexical_cast<T>(sval);
      } catch (...) {
         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
         if (myRank == MASTER_RANK) {
            std::cerr << __FILE__ << ":" << __LINE__
                      << std::string(" Problems casting ") + name + " " + sval + std::string(" to ") + typeid(T).name()
                      << std::endl;

            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
   }

   /** Get the value of the given parameter added with addComposing().
    * This may be called after having called Parse, and it may be called by any process, in any order.
    * Aborts on failed cast.
    * @param name The name of the parameter.
    * @param value A variable where the value of the parameter is written.
    */
   template <typename T> static void get(const std::string& name, std::vector<T>& value) {
      std::vector<std::string> stringValue;
      get(name, stringValue);

      for (std::vector<std::string>::iterator i = stringValue.begin(); i != stringValue.end(); ++i) {
         try {
            value.push_back(boost::lexical_cast<T>(*i));
         } catch (...) {
            int myRank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            if (myRank == MASTER_RANK) {
               std::cerr << __FILE__ << ":" << __LINE__
                         << std::string(" Problems casting ") + name + *i + std::string(" to ") + typeid(T).name()
                         << std::endl;

               MPI_Abort(MPI_COMM_WORLD, 1);
            }
         }
      }
   }

   // Determine whether a given variable has been set.
   static bool isSet(const std::string& name) {
      return(options.find(name) != options.end());
   }

   static void addComposing(const std::string& name, const std::string& desc);

   static void helpMessage();

   static bool versionMessage();

   static std::string versionInfo();
   
   static std::string configInfo();

   static bool parse(const bool needsRunConfig = true, const bool allowUnknown = true);

   static bool helpRequested;

private:
   static int argc;    /**< How many entries argv contains.*/
   static char** argv; /**< Pointer to char* array containing command line parameters.*/

   static boost::program_options::options_description* descriptions;
   static boost::program_options::variables_map* variables;

   static std::map<std::string, std::string> options;
   static std::map<std::string, bool> isOptionParsed;
   static std::map<std::string, std::vector<std::string>> vectorOptions;
   static std::map<std::string, bool> isVectorOptionParsed;

   static std::string global_config_file_name;
   static std::string user_config_file_name;
   static std::string run_config_file_name;

   static void addDefaultParameters();
};

#endif
