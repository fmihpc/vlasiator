#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>
#include <map>

#ifndef CPU_SORT_IDS_HPP
#define CPU_SORT_IDS_HPP


// Comparator function for sorting vector of pairs
template<typename ID> inline bool paircomparator( const std::pair<ID, uint> & l, const std::pair<ID, uint> & r ) {
   return l.first < r.first;
}



template<typename ID, typename LENGTH> inline void sortIds(const uint dimension, 
                                                           const LENGTH meshSize, 
                                                           const std::vector<ID>& ids,
							   std::map<ID, ID>& mapping){

  //sortedIds.resize(ids.size());
  //TODO conditionally parallel version?
  //#pragma omp parallel for
   for (uint i = 0; i < ids.size() ; ++i ) {
     const ID id = ids[i];
     
     if (id > meshSize[0] * meshSize[1] * meshSize[2])
       continue;
     
     switch( dimension ) {
     case 0: {
       const ID idMapped = id; // Mapping the block id to different coordinate system if dimension is not zero:
       mapping[idMapped] = id;
     }
       break;
     case 1: {
       // Do operation: 
       //   id = x + y*x_max + z*y_max*x_max 
       //=> id' = id - (x + y*x_max) + y + x*y_max = x + y*x_max + z*y_max*x_max - (x + y*x_max) + y + x*y_max
       //          = y + x*y_max + z*y_max*x_max
       const ID x_index = (id-1) % meshSize[0];
       const ID y_index = ((id-1) / meshSize[0]) % meshSize[1];
       
       // Mapping the block id to different coordinate system if dimension is not zero:
       const ID idMapped = id - (x_index + y_index*meshSize[0]) + y_index + x_index * meshSize[1];
       
       mapping[idMapped] = id;
     }
     break;
   case 2: {
     // Do operation: 
     //   id  = x + y*x_max + z*y_max*x_max 
     //=> id' = z + y*z_max + x*z_max*y_max
     const ID x_index = (id-1) % meshSize[0];
     const ID y_index = ((id-1) / meshSize[0]) % meshSize[1];
     const ID z_index = ((id-1) / (meshSize[0] * meshSize[1]));
     
     // Mapping the id id to different coordinate system if dimension is not zero:
     //const uint idMapped 
     //  = z_indice 
     //  + y_indice * meshSize[2] 
     //  + x_indice*meshSize[1]*meshSize[2];
     const ID idMapped = 1 + z_index + y_index*meshSize[2] + x_index*meshSize[1]*meshSize[2];
     mapping[idMapped] = id;
   }
     break;
   }
}
// Finally, sort the list of pairs
// std::sort( sortedIds.begin(), sortedIds.end(), paircomparator<ID> );

//for (auto j = 1; j <= mapping.size(); j++)
//  std::cout << j << ", " << mapping[j] << "\n";

}






template<typename ID, typename LENGTH>  void sortIdlistByDimension(const uint dimension, const LENGTH meshSize, 
								   std::vector<ID> & ids,
								   std::vector<uint> & columnIdOffsets,
								   std::vector<uint> & columnNumIds,
								   std::vector<uint> & setColumnOffsets,
								   std::vector<uint> & setNumColumns) {
  
  const uint nIds = ids.size();

   //sort Ids
   std::vector<std::pair<ID, ID> > sortedIdPairs;
   sortIds<ID, LENGTH>(dimension, meshSize, ids, sortedIdPairs);
   
   
   // Put in the sorted ids, and also compute column offsets and lengths:
   columnIdOffsets.push_back(0); //first offset
   setColumnOffsets.push_back(0); //first offset   
   uint prev_column_id, prev_dimension_id;

   for (uint i=0; i<nIds; ++i) {
      // identifies a particular column
      uint column_id = sortedIdPairs[i].first / meshSize[dimension];     
       
      // identifies a particular id in a column (along the dimension)
      uint dimension_id = sortedIdPairs[i].first % meshSize[dimension];
      
      //sorted list
      ids[i] = sortedIdPairs[i].second;
       
      if ( i > 0 &&  ( column_id != prev_column_id || dimension_id != (prev_dimension_id + 1) )){
         //encountered new column! For i=0, we already entered the correct offset (0).
         //We also identify it as a new column if there is a break in the column (e.g., gap between two populations)
         /*add offset where the next column will begin*/
         columnIdOffsets.push_back(i); 
         /*add length of the current column that now ended*/
         columnNumIds.push_back(columnIdOffsets[columnIdOffsets.size()-1] - columnIdOffsets[columnIdOffsets.size()-2]);

         if (column_id != prev_column_id ){
            //encountered new set of columns, add offset to new set starting at present column
            setColumnOffsets.push_back(columnIdOffsets.size() - 1);
            /*add length of the previous column set that ended*/
            setNumColumns.push_back(setColumnOffsets[setColumnOffsets.size()-1] - setColumnOffsets[setColumnOffsets.size()-2]);
         }
      }      
      prev_column_id = column_id;
      prev_dimension_id = dimension_id;                        
   }
   
   columnNumIds.push_back(nIds - columnIdOffsets[columnIdOffsets.size()-1]);
   setNumColumns.push_back(columnNumIds.size() - setColumnOffsets[setColumnOffsets.size()-1]);
}






#endif
