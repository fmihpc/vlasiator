/* This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute 
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
