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

#ifndef READPARAMETERS_H
#define READPARAMETERS_H
#include <limits>
#include <mpi.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "common.h"

struct Readparameters {
    Readparameters(int argc, char* argv[],MPI_Comm comm);
    static bool add(const std::string& name,const std::string& desc,const std::string& defValue);
    static bool add(const std::string& name,const std::string& desc,const bool& defValue);
    static bool add(const std::string& name,const std::string& desc,const int& defValue);
    static bool add(const std::string& name,const std::string& desc,const unsigned int& defValue);
    static bool add(const std::string& name,const std::string& desc,const float& defValue);
    static bool add(const std::string& name,const std::string& desc,const double& defValue);

    static bool get(const std::string& name,std::string& value);
    static bool get(const std::string& name,bool& value);
    static bool get(const std::string& name,int& value);
    static bool get(const std::string& name,unsigned int& value);
    static bool get(const std::string& name,unsigned long& value);
    static bool get(const std::string& name,float& value);
    static bool get(const std::string& name,double& value);

//Functions for composing options (can be defined multiple times and are all returned as a vector)
    static bool addComposing(const std::string& name,const std::string& desc);
    static bool get(const std::string& name,std::vector<std::string>& value);
    static bool get(const std::string& name,std::vector<int>& value);
    static bool get(const std::string& name,std::vector<float>& value);
    static bool get(const std::string& name,std::vector<double>& value);

    
    static bool finalize();
    static bool helpMessage();
    static bool isInitialized();
    static bool parse();
   
private:
    static int argc;                  /**< How many entries argv contains.*/
    static char** argv;              /**< Pointer to char* array containing command line parameters.*/
    static int rank;
    static MPI_Comm comm;

  
    /** Private default constructor to prevent incorrect initialization.*/
    Readparameters();
    static bool addDefaultParameters();
};

#endif
