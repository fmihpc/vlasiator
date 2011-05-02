#include <cstdlib>
#include <iostream>
#include <limits>

#include <mpilogger.h>

#include "eventcontainer.h"

extern MPILogger mpilogger;

using namespace std;

/** Constructor for EventToken. Writes invalid values 
 * to streamID and cellID member variables.
 */
EventContainer::EventToken::EventToken() {
   streamID = numeric_limits<unsigned int>::max();
   cellID = numeric_limits<unsigned int>::max();
}

/** Overloaded constructor for EventToken. Copies the 
 * given variables into internal variables.
 * @param streamID ID number of the CUDA stream associated with the CUDA event.
 * @param event CUDA event.
 * @param cellID ID number of the spatial cell the stream is computing.
 */
EventContainer::EventToken::EventToken(const unsigned int& streamID,const cudaEvent_t& event,const uint64_t& cellID) {
   this->streamID = streamID;
   this->event = event;
   this->cellID = cellID;
}

/** Destructor for EventToken. The destructor is empty.
 */
EventContainer::EventToken::~EventToken() {
   
}

/** Constructor for EventContainer. The constructor is empty.
 */
EventContainer::EventContainer() {
   
}

/** Destructor for EventContainer. Destroys all CUDA events in 
 * EventContainer.
 */
EventContainer::~EventContainer() {
   for (it = events.begin(); it != events.end(); ++it) {
      cudaEventDestroy(it->event);
   }
}

/** Query if the EventContainer is empty.
 * @return If true, the EventContainer is empty.
 */
bool EventContainer::empty() const {return events.empty();}

/**< Test if any of the events in EventContainer has completed. This function 
 * returns the first completed event, if possible.
 * @param streamID ID number of the stream associated with the completed event.
 * @param cellID Global ID of the spatial cell the stream was computing.
 * @return If true, a completed event was found and streamID, cellID contain valid values.
 */
bool EventContainer::getCompleted(unsigned int& streamID,uint64_t& cellID) {
   cudaError_t result;
   for (it = events.begin(); it != events.end(); ++it) {
      result = cudaEventQuery(it->event);
      if (result == cudaSuccess) {
	 if (cudaEventDestroy(it->event) != cudaSuccess) {
	    mpilogger << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred in event destruction!" << endl << write;
	 }
	 streamID = it->streamID;
	 cellID   = it->cellID;
	 events.erase(it);
	 return true;
      } else if (result != cudaErrorNotReady) {
	 mpilogger << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred in event query!" << endl << write;
      }
   }
   return false;
}

/** Insert a new event into EventContainer and associate it with the given CUDA stream.
 * @param streamID ID number of the CUDA stream associated with the CUDA event.
 * @param stream CUDA stream.
 * @param cellID Global ID of the cell the stream is computing.
 * @return If true, a new event was recorded for the given stream.
 */
bool EventContainer::insert(const unsigned int& streamID,cudaStream_t& stream,const uint64_t& cellID) {
   cudaEvent_t event;
   if (cudaEventCreate(&event) != cudaSuccess) {
      mpilogger << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while creating an event!" << endl << write;
   }
   
   if (cudaEventRecord(event,stream) != cudaSuccess) {
      mpilogger << "ERROR '" << cudaGetErrorString(cudaGetLastError()) << "' occurred while recording event with stream #" << streamID << endl << write;
   }
   events.push_back(EventContainer::EventToken(streamID,event,cellID));
   return true;
}











