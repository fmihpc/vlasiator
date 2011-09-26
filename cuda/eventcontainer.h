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

#ifndef EVENTCONTAINER_H
#define EVENTCONTAINER_H

#include <list>
#include <stdint.h>
#include <cuda.h>
#include <cuda_runtime.h>

/** Container class for CUDA events. CUDA events can be used 
 * to track the completion of streams' tasks.
 */
class EventContainer {
 public:
   EventContainer();
   ~EventContainer();

   bool empty() const;
   bool getCompleted(unsigned int& streamID,uint64_t& cellID);
   bool insert(const unsigned int& streamID,cudaStream_t& stream,const uint64_t& cellID);
   
 private:
   
   /** Struct EventToken is used to associate CUDA streams with 
    * CUDA events and with the spatial cell the stream is computing.
    */
   struct EventToken {
      unsigned int streamID; /**< ID number of the stream associated with the event.*/
      cudaEvent_t event;     /**< CUDA event.*/
      uint64_t cellID;       /**< ID number of the spatial cell which is computed (if applicable).*/
      
      EventToken();
      EventToken(const unsigned int& streamID,const cudaEvent_t& event,const uint64_t& cellID);
      ~EventToken();
   };
   
   std::list<EventContainer::EventToken> events;       /**< List of recorded CUDA events currently in EventContainer.*/
   std::list<EventContainer::EventToken>::iterator it; /**< Iterator to EventContainer::events.*/
};

#endif
