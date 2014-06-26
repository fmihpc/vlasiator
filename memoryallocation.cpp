/*
This file is part of Vlasiator.

Copyright 2014 Finnish Meteorological Institute
*/
#include <cstdlib>
#include <string.h>
#include <iostream>
#include "logger.h"
#include "memoryallocation.h"
#ifdef PAPI_MEM
#include "papi.h" 
#endif 

extern Logger logFile, diagnostic;
using namespace std;

#ifdef USE_JEMALLOC
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

// Global delete using jemalloc
void operator delete(void *p)
{
   je_free(p);
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

// Global delete[] using jemalloc
void operator delete[](void *p)
{
   je_free(p);
}

#endif 


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

/*! Measures memory consumption and writes it into logfile. Collective operation on MPI_COMM_WORLD
 */
void report_process_memory_consumption(){
   /*now report memory consumption into logfile*/
   int rank,n_procs;
   MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef PAPI_MEM
   /*If we have PAPI, we can report the resident usage of the process*/
   if (PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT) {
      PAPI_dmem_info_t dmem;  
      PAPI_get_dmem_info(&dmem);
      double mem_papi[2] = {};
      double sum_mem_papi[2];
      double min_mem_papi[2];
      double max_mem_papi[2];
      /*PAPI returns memory in KB units, transform to bytes*/
      mem_papi[0] = dmem.size * 1024; 
      mem_papi[1] = dmem.resident * 1024;
      MPI_Reduce(mem_papi, sum_mem_papi, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(mem_papi, min_mem_papi, 2, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(mem_papi, max_mem_papi, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      logFile << "(MEM) PAPI Resident (avg, min, max): " << sum_mem_papi[1]/n_procs << " " << min_mem_papi[1] << " "  << max_mem_papi[1] << endl;
   }   
#endif

   // Report /proc/meminfo memory consumption:
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
}



