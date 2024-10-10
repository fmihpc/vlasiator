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

/*! Return the amount of free memory on the node in bytes*/
uint64_t get_node_free_memory(){
   uint64_t mem_proc_free = 0;
   FILE * in_file = fopen("/proc/meminfo", "r");
   char attribute_name[200];
   int memory;
   char memory_unit[10];
   const char * memfree_attribute_name = "MemFree:";
   if( in_file ) {
      // Read free memory:
      while( fscanf( in_file, "%s %d %s", attribute_name, &memory, memory_unit ) != EOF ) {
         // Check if the attribute name equals memory free
         if( strcmp(attribute_name, memfree_attribute_name ) == 0 ) {
            //free memory in KB, transform to B
            mem_proc_free = (uint64_t)memory * 1024;
         }
      }
   }
   fclose( in_file );
   
   return mem_proc_free;
}

/*! Measures memory consumption and writes it into logfile. 
 *  Collective operation on MPI_COMM_WORLD
 *  extra_bytes is used for additional buffer for the high water mark, 
 *  for example when estimating refinement memory usage
 */
void report_process_memory_consumption(double extra_bytes){
   /*Report memory consumption into logfile*/

   char nodename[MPI_MAX_PROCESSOR_NAME]; 
   int namelength, nodehash;
   int rank, nProcs, nodeRank, interRank;
   int nNodes;
   const double GiB = pow(2,30);
   const double TiB = pow(2,40);

   hash<string> hasher; 
   MPI_Comm nodeComm;
   MPI_Comm interComm;
   

   MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
   //get name of this node
   MPI_Get_processor_name(nodename,&namelength);   
   nodehash=(int)(hasher(string(nodename)) % std::numeric_limits<int>::max());
   
   //intra-node communicator
   MPI_Comm_split(MPI_COMM_WORLD, nodehash, rank, &nodeComm);
   MPI_Comm_rank(nodeComm,&nodeRank);
   //create communicator for inter-node communication
   MPI_Comm_split(MPI_COMM_WORLD, nodeRank, rank, &interComm);
   MPI_Comm_rank(interComm, &interRank);
   MPI_Comm_size(interComm, &nNodes);

#ifdef PAPI_MEM
   /*If we have PAPI, we can report the resident usage of the process*/
   if (PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT) {
      PAPI_dmem_info_t dmem;  
      PAPI_get_dmem_info(&dmem);
      double mem_papi[4] = {};
      double node_mem_papi[4] = {};
      double sum_mem_papi[4];
      double min_mem_papi[4];
      double max_mem_papi[4];
      /*PAPI returns memory in KB units, transform to bytes*/
      mem_papi[0] = dmem.high_water_mark * 1024;
      mem_papi[1] = dmem.high_water_mark * 1024 + extra_bytes;
      mem_papi[2] = dmem.resident * 1024;
      mem_papi[3] = extra_bytes;
      //sum node mem
      MPI_Reduce(mem_papi, node_mem_papi, 4, MPI_DOUBLE, MPI_SUM, 0, nodeComm);
      
      //rank 0 on all nodes do total reduces
      if(nodeRank == 0) {
         MPI_Reduce(node_mem_papi, sum_mem_papi, 4, MPI_DOUBLE, MPI_SUM, 0, interComm);
         MPI_Reduce(node_mem_papi, min_mem_papi, 4, MPI_DOUBLE, MPI_MIN, 0, interComm);
         MPI_Reduce(node_mem_papi, max_mem_papi, 4, MPI_DOUBLE, MPI_MAX, 0, interComm);
         if (max_mem_papi[3] != 0.0) {
            logFile << "(MEM) Estimating increased high water mark from refinement" << endl;
         }
         logFile << "(MEM) tstep " << Parameters::tstep << " t " << Parameters::t << " Resident per node (avg, min, max): " << sum_mem_papi[2]/nNodes/GiB << " " << min_mem_papi[2]/GiB << " "  << max_mem_papi[2]/GiB << endl;
         logFile << "(MEM) tstep " << Parameters::tstep << " t " << Parameters::t << " High water mark per node                 (GiB) avg: " << sum_mem_papi[0]/nNodes/GiB << " min: " << min_mem_papi[0]/GiB << " max: "  << max_mem_papi[0]/GiB <<
            " sum (TiB): " << sum_mem_papi[0]/TiB << " on "<< nNodes << " nodes" << endl;
         if(max_mem_papi[3] != 0.0) {
            logFile << "(MEM) tstep " << Parameters::tstep << " t " << Parameters::t << " High water mark per node with refinement (GiB) avg: " << sum_mem_papi[1]/nNodes/GiB << " min: " << min_mem_papi[1]/GiB << " max: "  << max_mem_papi[1]/GiB <<
               " sum (TiB): " << sum_mem_papi[1]/TiB << " on "<< nNodes << " nodes" << endl;
         }
      }
      if(rank == MASTER_RANK) {
         bailout(max_mem_papi[1]/GiB > Parameters::bailout_max_memory, "Memory high water mark per node exceeds bailout threshold", __FILE__, __LINE__);
      }
   }
   
#endif


   /*
   // Report /proc/meminfo memory consumption.      
   double mem_proc_free = (double)get_node_free_memory();
   double total_mem_proc = 0;
   double min_free,max_free;
   const int root = 0;
   const int numberOfParameters = 1;
   MPI_Reduce( &mem_proc_free, &total_mem_proc, numberOfParameters, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD );
   MPI_Reduce( &mem_proc_free, &min_free, numberOfParameters, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD );
   MPI_Reduce( &mem_proc_free, &max_free, numberOfParameters, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
   logFile << "(MEM) Node free memory (avg, min, max): " << total_mem_proc/n_procs << " " << min_free << " " << max_free << endl;
   logFile << writeVerbose;
   */

   MPI_Comm_free(&interComm);
   MPI_Comm_free(&nodeComm);

}



