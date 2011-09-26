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

#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <stdint.h>
#include <queue>
#include <vector>

namespace priorityqueue {
   /** Struct used to store spatial cell IDs, along with their priorities, 
    * into a priority queue. 
    */
   struct CellPriority {
      unsigned int priority; /**< Priority of the spatial cell.*/
      uint64_t cellID;       /**< Global ID of the spatial cell.*/
      
      /** Constructor for CellPriority. The constructor just assings given priority and 
       * cellID values to the corresponding member variables.
       * @param priority The priority of the spatial cell.
       * @param cellID Global ID of the spatial cell.
       */
      CellPriority(const unsigned int& priority,const uint64_t& cellID): priority(priority),cellID(cellID) { }
   };

   /** Comparator struct used to determine if spatial cell a has a lower priority 
    * than cell b (weak ordering).
    */
   struct Comparator {
      /** Compare the given CellPrioritys. If a has a lower priority than b, return true. 
       * Otherwise return false.
       * @param a First spatial cell to compare.
       * @param b Second spatial cell to compare.
       * @return If true, a < b.
       */
      bool operator()(const CellPriority& a,const CellPriority& b) {return a.priority < b.priority;}
   };
}

/** Priority queue for computing spatial cells in order of importance. 
 * The "importance" is defined by the user, priority queues just gives out 
 * the cell with the highest priority. This class is just a wrapper for 
 * STL priority_queue, which requires one to supply a comparator object.
 */
template<typename T> class PriorityQueue {
 public:
   void clear();
   bool empty() const;
   void insert(const T& cellID,const unsigned int& priority);
   bool pop(T& cellID,uint& priority);
   size_t size() {return priorityQueue.size();}
   
 private:
 
   std::priority_queue<priorityqueue::CellPriority,std::vector<priorityqueue::CellPriority>,priorityqueue::Comparator> priorityQueue;
};

/** Clear the contents of the PriorityQueue.
 */
template<typename T> void PriorityQueue<T>::clear() {
   while (priorityQueue.empty() == false) priorityQueue.pop();
}

/** Test if PriorityQueue is empty.
 * @return If true, the PriorityQueue is empty.
 */
template<typename T> bool PriorityQueue<T>::empty() const {return priorityQueue.empty();}

/** Insert a spatial cell into PriorityQueue with the given priority.
 * @param cellID Global ID of the spatial cell.
 * @param priority Computation priority of the cell.
 */
template<typename T> void PriorityQueue<T>::insert(const T& cellID,const unsigned int& priority) {
   priorityQueue.push(priorityqueue::CellPriority(priority,cellID));
}

/** Pop the cell with the highest priority from the PriorityQueue.
 * @param cellID A variable in which the global ID of the popped cell is written into. 
 * @param priority A variable in which the priority of the popped cell is written into.
 * @return If true, cellID and priority contain valid values.
 */
template<typename T> bool PriorityQueue<T>::pop(T& cellID,uint& priority) {
   if (priorityQueue.empty() == true) return false;
   cellID = priorityQueue.top().cellID;
   priority = priorityQueue.top().priority;
   priorityQueue.pop();
   return true;
}

#endif
