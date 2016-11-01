/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
 * 
 * File:   vlsv_util.cpp
 * Author: sandroos
 *
 * Created on July 6, 2015, 12:56 PM
 */

#include <cstdlib>
#include <iostream>
#include <dirent.h>

#include "vlsv_util.h"

using namespace std;

std::vector<std::string> toolutil::getFiles(const std::string& mask) {
   vector<string> fileList;
   const string directory = ".";
   const string suffix = ".vlsv";
   DIR* dir = opendir(directory.c_str());
   if (dir == NULL) {
      cerr << "ERROR in reading directory contents in " << __FILE__ << ":" << __LINE__ << endl;
      closedir(dir);
      return fileList;
   }

   struct dirent* entry = readdir(dir);
   while (entry != NULL) {
      const string entryName = entry->d_name;
      if (entryName.find(mask) == string::npos || entryName.find(suffix) == string::npos) {
         entry = readdir(dir);
         continue;
      }
      fileList.push_back(entryName);
      entry = readdir(dir);
   }
   closedir(dir);
   return fileList;
}
