/* This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute 
 * 
 * File:   vlsv_util.h
 * Author: sandroos
 *
 * Created on July 6, 2015, 12:56 PM
 */

/** @file vlsv_util.h
 * This file contains miscellaneous utility routines that 
 * are used by vlsv tools.*/

#include <string>
#include <vector>

#ifndef VLSV_UTIL_H
#define	VLSV_UTIL_H

namespace toolutil {
   
   std::vector<std::string> getFiles(const std::string& mask);

} // namespace toolutil

#endif	// VLSV_UTIL_H

