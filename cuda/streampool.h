#ifndef STREAMPOOL_H
#define STREAMPOOL_H

#include <list>
#include <cuda.h>
#include <cuda_runtime.h>

/** A class for maintaining CUDA streams. Using this class a user can 
 * request available streams for operations. After streams have 
 * completed their computations the streams need to be returned to the 
 * StreamPool.
 */
class StreamPool {
 public:
   StreamPool();
   ~StreamPool();
   
   unsigned int available() const;
   bool get(unsigned int& streamID,cudaStream_t& stream,const bool& mayEmptyPool = true);
   cudaStream_t& get(unsigned int& streamID);
   bool insert(const unsigned int& streamID);
   bool initialize(const unsigned int& streams,const unsigned int& threshold);
   
 private:
   std::list<unsigned int> availableStreams; /**< List of available stream IDs.*/
   bool initialized;                         /**< If true, StreamPool has initialized correctly and is ready for use.*/
   unsigned int N_streams;                   /**< Total number of streams in StreamPool.*/
   cudaStream_t* streams;                    /**< Pointers to CUDA streams. Array size is given by N_streams.*/
   unsigned int threshold;                   /**< Sometimes a user might not want to empty StreamPool completely.
					      * If number of available streams <= threshold, a stream ID is not 
					      * returned by StreamPool::get if mayEmptyPool parameter is false.*/
};

#endif

