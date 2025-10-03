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

#include "common.h"
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>

#ifdef _OPENMP
  #include <omp.h>
#endif
#include "grid.h"
#include "logger.h"
#include "object_wrapper.h"

#include "memory_report.h"

#ifdef PAPI_MEM
#include "papi.h"
#endif

#ifdef USE_GPU
#include "arch/gpu_base.hpp"
#endif
extern Logger logFile;

/*! Report spatial cell counts per refinement level as well as velocity cell counts per population into logfile
 */
void report_cell_and_block_counts(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid){
   cint maxRefLevel = mpiGrid.get_maximum_refinement_level();
   const std::vector<CellID> localCells = getLocalCells();
   cint popCount = getObjectWrapper().particleSpecies.size();

   // popCount+1 as we store the spatial cell counts and then the populations' v_cell counts.
   // maxRefLevel+1 as e.g. there's 2 levels at maxRefLevel == 1
   std::vector<int64_t> localCounts((popCount+1)*(maxRefLevel+1), 0), globalCounts((popCount+1)*(maxRefLevel+1), 0);

   for (const auto cellid : localCells) {
      cint level = mpiGrid.get_refinement_level(cellid);
      localCounts[level]++;
      for(int pop=0; pop<popCount; pop++) {
         localCounts[maxRefLevel+1 + level*popCount + pop] += mpiGrid[cellid]->get_number_of_velocity_blocks(pop);
      }
   }

   MPI_Reduce(localCounts.data(), globalCounts.data(), (popCount+1)*(maxRefLevel+1), MPI_INT64_T, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);

   logFile << "(CELLS) tstep = " << P::tstep << " time = " << P::t << " spatial cells [ ";
   for(int level = 0; level <= maxRefLevel; level++) {
      logFile << globalCounts[level] << " ";
   }
   logFile << "] blocks ";
   for(int pop=0; pop<popCount; pop++) {
      logFile << getObjectWrapper().particleSpecies[pop].name << " [ ";
      for(int level = 0; level <= maxRefLevel; level++) {
         logFile << globalCounts[maxRefLevel+1 + level*popCount + pop] << " ";
      }
      logFile << "] ";
   }
   logFile << std::endl << std::flush;

}

/*! Return the amount of free memory on the node in bytes*/
uint64_t get_node_free_memory(){
   uint64_t mem_proc_free = 0;
   FILE * in_file = fopen("/proc/meminfo", "r");
   char attribute_name[200];
   uint64_t memory=0;
   char memory_unit[10];
   const char * memfree_attribute_name = "MemFree:";
   if (in_file) {
      while (!feof(in_file)) {
         int retval = fscanf(in_file, "%199s %lu %9s", attribute_name, &memory, memory_unit);
         if (retval >= 2) {
            if (strcmp(attribute_name, memfree_attribute_name) == 0) {
               // free memory in KB, transform to B
               mem_proc_free = memory * 1024;
               break;
            }
         } else {
            // Skip to neww line
            fscanf(in_file, "%*[^\n]\n");
         }
      }
      fclose(in_file);
   }

   return mem_proc_free;
}

/*! Measures memory consumption and writes it into logfile.
 *  Collective operation on MPI_COMM_WORLD
 *  extra_bytes is used for additional buffer for the high water mark,
 *  for example when estimating refinement memory usage
 */
void report_memory_consumption(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   double extra_bytes/*=0*/
){
   /*Report memory consumption into logfile*/

   char nodename[MPI_MAX_PROCESSOR_NAME];
   int namelength, nodehash;
   int rank, nProcs, nodeRank, interRank;
   int nNodes;
   const double GiB = pow(2,30);

   std::hash<std::string> hasher;
   MPI_Comm nodeComm;
   MPI_Comm interComm;


   MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   //get name of this node
   MPI_Get_processor_name(nodename,&namelength);
   nodehash=(int)(hasher(std::string(nodename)) % std::numeric_limits<int>::max());

   //intra-node communicator
   MPI_Comm_split(MPI_COMM_WORLD, nodehash, rank, &nodeComm);
   MPI_Comm_rank(nodeComm,&nodeRank);
   //create communicator for inter-node communication
   MPI_Comm_split(MPI_COMM_WORLD, nodeRank, rank, &interComm);
   MPI_Comm_rank(interComm, &interRank);
   MPI_Comm_size(interComm, &nNodes);

   // Report /proc/meminfo memory consumption first so we then get resident and HWM from Papi below, easier to compare by eye in logfile in this order.
   double mem_proc_free = (double)get_node_free_memory();
   double total_mem_proc = 0;
   double min_free,max_free;
   const int numberOfParameters = 1;
   MPI_Reduce( &mem_proc_free, &total_mem_proc, numberOfParameters, MPI_DOUBLE, MPI_SUM, 0, interComm );
   MPI_Reduce( &mem_proc_free, &min_free, numberOfParameters, MPI_DOUBLE, MPI_MIN, MASTER_RANK, MPI_COMM_WORLD );
   MPI_Reduce( &mem_proc_free, &max_free, numberOfParameters, MPI_DOUBLE, MPI_MAX, MASTER_RANK, MPI_COMM_WORLD );

   char reportstring[512];
   snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/node; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g on %i nodes\n",
         P::tstep, P::t, "Free", total_mem_proc/nNodes/GiB, min_free/GiB, max_free/GiB, total_mem_proc/GiB, nNodes);
   logFile << reportstring;

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
      mem_papi[1] = dmem.resident * 1024 + extra_bytes;
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
            logFile << "(MEM) Estimating increased high water mark from refinement" << std::endl;
         }
         snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/node; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g on %i nodes\n",
                P::tstep, P::t, "Resident",            sum_mem_papi[2]/nNodes/GiB, min_mem_papi[2]/GiB, max_mem_papi[2]/GiB, sum_mem_papi[2]/GiB, nNodes);
         logFile << reportstring;
         snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/node; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g on %i nodes\n",
               P::tstep, P::t, "High water mark \U0001F30A   ",     sum_mem_papi[0]/nNodes/GiB, min_mem_papi[0]/GiB, max_mem_papi[0]/GiB, sum_mem_papi[0]/GiB, nNodes);
         logFile << reportstring;
         if (max_mem_papi[3] != 0.0) {
            snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/node; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g on %i nodes\n",
               P::tstep, P::t, "Resident with refines", sum_mem_papi[1]/nNodes/GiB, min_mem_papi[1]/GiB, max_mem_papi[1]/GiB, sum_mem_papi[1]/GiB, nNodes);
            logFile << reportstring;
         }
      }
      if(rank == MASTER_RANK) {
         bailout(max_mem_papi[0]/GiB > P::bailout_max_memory, "Memory high water mark per node exceeds bailout threshold", __FILE__, __LINE__);
         bailout(max_mem_papi[1]/GiB > P::bailout_max_memory, "Estimated resident per node exceeds bailout threshold", __FILE__, __LINE__);
      }
   }
#endif

   /*now report memory consumption of mpiGrid specifically into logfile*/
   const std::vector<CellID>& cells = getLocalCells();
   const std::vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary();

   /* Compute memory statistics of the memory consumption of the spatial cells.
    * Internally we use double as MPI does
    * not define proper uint64_t datatypes for MAXLOCNot Real, as we
    * want double here not to loose accuracy.
    */

   /*report data for memory needed by blocks*/
   double mem[6] = {0};
   double sum_mem[6];

   for(unsigned int i=0;i<cells.size();i++){
      #ifdef USE_GPU
      // GPU version keeps track of largest attained velocity mesh throughout acceleration cycles
      mem[0] += mpiGrid[cells[i]]->largestvmesh * WID3 * sizeof(Realf);
      #else
      mem[0] += mpiGrid[cells[i]]->get_cell_memory_size();
      #endif
      mem[3] += mpiGrid[cells[i]]->get_cell_memory_capacity();
   }

   for(unsigned int i=0;i<remote_cells.size();i++){
      if(mpiGrid[remote_cells[i]] != NULL) {
         mem[1] += mpiGrid[remote_cells[i]]->get_cell_memory_size();
         mem[4] += mpiGrid[remote_cells[i]]->get_cell_memory_capacity();
      }
   }

   mem[2] = mem[0] + mem[1];//total memory according to size()
   mem[5] = mem[3] + mem[4];//total memory according to capacity()


   MPI_Reduce(mem, sum_mem, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


   struct {
      double val;
      int   rank;
   } max_mem[3],mem_usage_loc[3],min_mem[3];
   for(uint i = 0; i<3; i++){
      mem_usage_loc[i].val = mem[i + 3]; //report on capacity numbers (6: local cells, 7: remote cells, 8: all cells)
      mem_usage_loc[i].rank = rank;
   }

   MPI_Reduce(mem_usage_loc, max_mem, 3, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   MPI_Reduce(mem_usage_loc, min_mem, 3, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);


   snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/rank; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g min rank %i max rank %i\n",
         P::tstep, P::t, "Local cells capacity",  sum_mem[3]/nProcs/GiB, min_mem[0].val/GiB, max_mem[0].val/GiB, sum_mem[3]/GiB, min_mem[0].rank, max_mem[0].rank);
   logFile << reportstring;
   snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/rank; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g min rank %i max rank %i\n",
         P::tstep, P::t, "Remote cells capacity", sum_mem[4]/nProcs/GiB, min_mem[1].val/GiB, max_mem[1].val/GiB, sum_mem[4]/GiB, min_mem[1].rank, max_mem[1].rank);
   logFile << reportstring;
   snprintf(reportstring,512, "(MEM) tstep %i t %.3g %-21s (GiB/rank; avg, min, max, sum): %-8.3g %-8.3g %-8.3g %-8.3g min rank %i max rank %i\n",
         P::tstep, P::t, "Total cells capacity",  sum_mem[5]/nProcs/GiB, min_mem[2].val/GiB, max_mem[2].val/GiB, sum_mem[5]/GiB, min_mem[2].rank, max_mem[2].rank);
   logFile << reportstring;

   snprintf(reportstring,512, "(MEM) tstep %i t %.3g Total size and capacity of SpatialCells        (GiB): %-8.3g %-8.3g\n",
         P::tstep, P::t, sum_mem[2]/GiB, sum_mem[5]/GiB);
   logFile << reportstring;

   logFile << writeVerbose;

   MPI_Comm_free(&interComm);
   MPI_Comm_free(&nodeComm);

   #ifdef USE_GPU
   // TODO: Clear duplicate output
   // (local_cells_capacity, ghost_cells_capacity, local_cells_size, ghost_cells_size)
   gpu_reportMemory(mem[3], mem[4], mem[0], mem[1]);

   #endif
}
