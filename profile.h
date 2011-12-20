/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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

#ifndef PROFILE_H
#define PROFILE_H

#include "string"
#include "vector"
#include "mpi.h"

namespace profile 
{
   //initialize a timer, with a particular label belonging to some groups    
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const std::string &label);
   int initializeTimer(const std::string &label,const std::string &group1);
   int initializeTimer(const std::string &label,const std::string &group1,const std::string &group2);
   int initializeTimer(const std::string &label,const std::string &group1,const std::string &group2,const std::string &group3);
   int initializeTimer(const std::string &label,const std::vector<std::string> &groups);


   //get id number of a timer, return -1 if it does not exist
   int getId(const std::string &label);
   
   //Start a profiling block with a certain label
   //If timer does not exist, then it is created
   bool start(const std::string &label);
   //start a profiling block with a certain id
   bool start(int id);

   //stop a profiling block with a certain label
   bool stop (const std::string &label,
              double workUnits=-1.0,
              const std::string &workUnitLabel="");
   //stop a profiling block with a certain id
   bool stop (int id,
              double workUnits=-1.0,
              const std::string &workUnitLabel="");

   bool print(MPI_Comm comm,double minParentFraction=0.0);
}

#endif

