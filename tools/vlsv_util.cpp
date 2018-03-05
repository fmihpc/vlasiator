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
 * 
 * File:   vlsv_util.cpp
 * Author: sandroos
 *
 * Created on July 6, 2015, 12:56 PM
 */

#include <cstdlib>
#include <iostream>
#include <glob.h>

#include "vlsv_util.h"

using namespace std;

std::vector<std::string> toolutil::getFiles(const std::string& mask) {
   vector<string> fileList;

   glob_t glob_matches;
   int retval = glob(mask.c_str(), GLOB_TILDE_CHECK, NULL, &glob_matches);
   if(retval != 0) {
      switch(retval) {
         case GLOB_NOSPACE:
            cerr << "ERROR enumerating filenames: out of memory in " << __FILE__ << ":" << __LINE__ << endl;
            break;
         case GLOB_ABORTED:
            cerr << "ERROR in reading directory contents in " << __FILE__ << ":" << __LINE__ << endl;
            break;
         case GLOB_NOMATCH:
            cerr << "ERROR: no matching file found for pattern " << mask << endl;
            break;
      }
      globfree(&glob_matches);
      return fileList;
   }
   for(unsigned int i=0; i<glob_matches.gl_pathc; i++) {
      fileList.push_back(glob_matches.gl_pathv[i]);
   }

   globfree(&glob_matches);
   return fileList;
}
