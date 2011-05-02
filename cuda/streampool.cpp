#include <cstdlib>
#include <iostream>

#include "streampool.h"

using namespace std;

/** Constructor for class StreamPool. The constructor is empty - 
 * one must call StreamPool::initialze in order to initialize 
 * the StreamPool.
 */
StreamPool::StreamPool(): streams(NULL) {

}

/** Destructor for class StreamPool. Deallocates CUDA streams.
 */
StreamPool::~StreamPool() {
   if (streams != NULL) {
      for (unsigned int i=0; i<N_streams; ++i) {
	 cudaStreamDestroy(streams[i]);
      }
   }
   delete streams;
   streams = NULL;
}

/** Get the number of available streams.
 * @return Number of streams currently in StreamPool.
 */
unsigned int StreamPool::available() const {return availableStreams.size();}

/** Initialize StreamPool.
 * @param N_streams Number of CUDA streams to create for StreamPool.
 * @param threshold Threshold number of streams. If available streams is less than or 
 * equal this value, StreamPool::get does not return a stream ID if parameter mayEmptyPool
 * is false.
 * @return If true, StreamPool initialized correctly and is ready for use.
 */
bool StreamPool::initialize(const unsigned int& N_streams,const unsigned int& threshold) {
   this->N_streams = N_streams;
   this->threshold = threshold;
   
   initialized = true;
   streams = new cudaStream_t[N_streams];
   for (unsigned int i=0; i<N_streams; ++i) {
      if (cudaStreamCreate(&(streams[i])) != cudaSuccess) {
	 initialized = false;
	 cerr << "ERROR " << cudaGetErrorString(cudaGetLastError()) << " occurred while creating streams!" << endl;
      } else {
	 availableStreams.push_back(i);
      }
   }
   return initialized;
}

/** Get the CUDA stream with the given stream ID.
 * @param streamID ID number of the requested stream.
 * @return CUDA stream.
 */
cudaStream_t& StreamPool::get(unsigned int& streamID) {return streams[streamID];}

/** Get an available CUDA stream from StreamPool. 
 * @param streamID ID number of the returned stream is written into this variable.
 * @param stream CUDA stream is copied into this variable.
 * @param mayEmptyPool If true, StreamPool::get is allowed to empty the StreamPool completely.
 * @return If true, a valid stream was returned.
 */
bool StreamPool::get(unsigned int& streamID,cudaStream_t& stream,const bool& mayEmptyPool) {
   if (availableStreams.size() == 0) return false;
   if (mayEmptyPool == false && availableStreams.size() < threshold) return false;
   
   streamID = availableStreams.front();
   availableStreams.pop_front();
   stream = streams[streamID];
   return true;
}

/** Return a CUDA stream to StreamPool.
 * @param streamID ID number of the CUDA stream. The ID number was 
 * given in StreamPool::get member function.
 * @return If true, the stream was returned correctly.
 */
bool StreamPool::insert(const unsigned int& streamID) {
   if (streamID >= N_streams) return false;
   availableStreams.push_back(streamID);
   return true;
}
