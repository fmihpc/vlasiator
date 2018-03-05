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

