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
#include <limits>
#include <mpi.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "common.h"

using std::string;
using std::vector;

struct Readparameters {
   Readparameters(int argc, char *argv[], MPI_Comm comm);
   static void add(const string &name, const string &desc, const string &defValue);
   static void add(const string &name, const string &desc, const bool &defValue);
   static void add(const string &name, const string &desc, const int &defValue);
   static void add(const string &name, const string &desc, const unsigned int &defValue);
   static void add(const string &name, const string &desc, const float &defValue);
   static void add(const string &name, const string &desc, const double &defValue);

   static void get(const string &name, string &value);
   static void get(const string &name, bool &value);
   static void get(const string &name, int &value);
   static void get(const string &name, unsigned int &value);
   static void get(const string &name, unsigned long &value);
   static void get(const string &name, float &value);
   static void get(const string &name, double &value);

   // Functions for composing options (can be defined multiple times and are all returned as a vector)
   static void addComposing(const string &name, const string &desc);
   static void get(const string &name, vector<string> &value);
   static void get(const string &name, vector<int> &value);
   static void get(const string &name, vector<unsigned int> &value);
   static void get(const string &name, vector<float> &value);
   static void get(const string &name, vector<double> &value);

   static void finalize();
   static void helpMessage();
   static bool versionMessage();
   static bool isInitialized();
   static bool parse(const bool needsRunConfig = true);

   static bool helpRequested;

private:
   static int argc;    /**< How many entries argv contains.*/
   static char **argv; /**< Pointer to char* array containing command line parameters.*/
   static int rank;
   static MPI_Comm comm;

   /** Private default constructor to prevent incorrect initialization.*/
   Readparameters();
   static void addDefaultParameters();
};

#endif
