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
#include <cstdlib>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unordered_map> // for hasher
#include <limits>
#include "logger.h"
#include "memoryallocation.h"
#include "common.h"
#include "parameters.h"
#ifdef PAPI_MEM
#include "papi.h" 
#endif 

extern Logger logFile, diagnostic;
using namespace std;

#ifdef USE_JEMALLOC
#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

// Declare global new etc. only if using old legacy versions of jemalloc, prior to 5.0.0
#if JEMALLOC_VERSION_MAJOR < 5

// Global new using jemalloc
void *operator new(size_t size)
{
   void *p;
   p =  je_malloc(size);
   if(!p) {
      bad_alloc ba;
      throw ba;
   }
   return p;
}

// Global new[] using jemalloc
void *operator new[](size_t size)
{
   void *p;
   p =  je_malloc(size);
   if(!p) {
      bad_alloc ba;
      throw ba;
   }
   return p;
}

// Global delete using jemalloc
void operator delete(void *p)
{
   je_free(p);
}

// Global delete[] using jemalloc
void operator delete[](void *p)
{
   je_free(p);
}

#if __cpp_sized_deallocation >= 201309
void operator delete(void *ptr, std::size_t size) noexcept {
   je_sdallocx(ptr, size, /*flags=*/0);
}
void operator delete[](void *ptr, std::size_t size) noexcept {
   je_sdallocx(ptr, size, /*flags=*/0);
}
#endif  // __cpp_sized_deallocation
#endif // JEMALLOC_VERSION_MAJOR < 5
#endif // use jemalloc

/*! Purge allocations from all arenas to actually release memory back to system */
void memory_purge() {
#ifdef USE_JEMALLOC
   je_mallctl("arena." STRINGIFY(MALLCTL_ARENAS_ALL) ".purge", NULL, NULL, NULL, 0);
#endif
}

/*! Initialize memory allocator configuration.*/
void memory_configurator() {
#ifdef USE_JEMALLOC
   bool logResult = false;
   bool foo {false};
   size_t bar {1};
   if (logResult) {
      // Read initial value
      je_mallctl("background_thread", &foo, &bar, NULL, 0);
      logFile << "(MEM) mallctl: background_thread value was ";
      if (foo) {
         logFile << "true";
      } else {
         logFile << "false";
      }
   }
   // Set background threads to true
   foo = true;
   bar = 1;
   je_mallctl("background_thread", NULL, NULL, &foo, bar);
   if (logResult) {
      // Read updated value
      je_mallctl("background_thread", &foo, &bar, NULL, 0);
      logFile << ", now set to ";
      if (foo) {
         logFile << "true";
      } else {
         logFile << "false";
      }
      logFile << "." << endl;
   }
#endif
}
