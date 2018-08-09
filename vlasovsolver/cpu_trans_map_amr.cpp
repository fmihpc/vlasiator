//#include "vlasovsolver/cpu_1d_ppm_nonuniform.hpp"
#include "vlasovsolver/cpu_1d_ppm_nonuniform_conserving.hpp"

struct grid_data {

   int value = 0;

   std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
   {
      return std::make_tuple(this, 0, MPI_BYTE);
   }
    
};

struct setOfPencils {

   uint N; // Number of pencils in the set
   uint sumOfLengths;
   std::vector<uint> lengthOfPencils; // Lengths of pencils
   std::vector<CellID> ids; // List of cells
   std::vector<Real> x,y; // x,y - position 

   setOfPencils() {
      N = 0;
      sumOfLengths = 0;
   }

   void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn, vector<Real> zIn) {

      N += 1;
      sumOfLengths += idsIn.size();
      lengthOfPencils.push_back(idsIn.size());
      ids.insert(ids.end(),idsIn.begin(),idsIn.end());
      x.push_back(xIn);
      y.push_back(yIn);
      z.insert(z.end(),zIn.begin(),zIn.end());
  
   }

   std::vector<CellID> getIds(uint pencilId) {

      if (pencilId > N) {
         return;
      }

      CellID ibeg = 0;
      for (uint i = 0; i < pencilId; i++) {
         ibeg += lengthOfPencils[i];
      }
      CellID iend = ibeg + lengthOfPencils[pencilId];

      vector<CellID> idsOut;
    
      for (uint i = ibeg; i <= iend; i++) {
         idsOut.push_back(ids[i]);
      }

      return idsOut;
   }

};

CellID selectNeighbor(dccrg::Dccrg<grid_data> &grid, CellID id, int dimension = 0, uint path = 0) {

   const auto neighbors = grid.get_face_neighbors_of(id);

   vector < CellID > myNeighbors;
   // Collect neighbor ids in the positive direction of the chosen dimension.
   // Note that dimension indexing starts from 1 (of course it does)
   for (const auto cell : neighbors) {
      if (cell.second == dimension + 1)
         myNeighbors.push_back(cell.first);
   }

   CellID neighbor;
  
   switch( myNeighbors.size() ) {
      // Since refinement can only increase by 1 level the only possibilities
      // Should be 0 neighbors, 1 neighbor or 4 neighbors.
   case 0 : {
      // did not find neighbors
      neighbor = INVALID_CELLID;
      break;
   }
   case 1 : {
      neighbor = myNeighbors[0];
      break;
   }
   case 4 : {
      neighbor = myNeighbors[path];
      break;
   }
   default: {
      // something is wrong
      neighbor = INVALID_CELLID;
      throw "Invalid neighbor count!";
      break;
   }
   }

   return neighbor;
  
}

setOfPencils buildPencilsWithNeighbors( dccrg::Dccrg<grid_data> &grid, 
					setOfPencils &pencils, CellID startingId,
					vector<CellID> ids, uint dimension, 
					vector<uint> path) {

   const bool debug = false;
   CellID nextNeighbor;
   uint id = startingId;
   uint startingRefLvl = grid.get_refinement_level(id);

   if( ids.size() == 0 )
      ids.push_back(startingId);

   // If the cell where we start is refined, we need to figure out which path
   // to follow in future refined cells. This is a bit hacky but we have to
   // use the order or the children of the parent cell to figure out which
   // corner we are in.
   // Maybe you could use physical coordinates here?
   if( startingRefLvl > path.size() ) {
      for ( uint i = path.size(); i < startingRefLvl; i++) {
         auto parent = grid.get_parent(id);
         auto children = grid.get_all_children(parent);
         auto it = std::find(children.begin(),children.end(),id);
         auto index = std::distance(children.begin(),it);
         auto index2 = index;
      
         switch( dimension ) {
         case 0: {
            index2 = index / 2;
            break;
         }
         case 1: {
            index2 = index - index / 2;
            break;
         }
         case 2: {
            index2 = index % 4;
            break;
         }
         }
         path.insert(path.begin(),index2);
         id = parent;
      }
   }

   id = startingId;
  
   while (id > 0) {

      // Find the refinement level in the neighboring cell. Any neighbor will do
      // since refinement level can only increase by 1 between neighbors.
      nextNeighbor = selectNeighbor(grid,id,dimension);

      // If there are no neighbors, we can stop.
      if (nextNeighbor == 0)
         break;
    
      uint refLvl = grid.get_refinement_level(nextNeighbor);

      if (refLvl > 0) {
    
         // If we have encountered this refinement level before and stored
         // the path this builder follows, we will just take the same path
         // again.
         if ( path.size() >= refLvl ) {
      
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << " ";
               std::cout << std::endl;
            }
	
            nextNeighbor = selectNeighbor(grid,id,dimension,path[refLvl-1]);      
	
         } else {
	
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have NOT seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << ' ';
               std::cout << std::endl;
            }
	
            // New refinement level, create a path through each neighbor cell
            for ( uint i : {0,1,2,3} ) {
	  
               vector < uint > myPath = path;
               myPath.push_back(i);
	  
               nextNeighbor = selectNeighbor(grid,id,dimension,myPath.back());
	  
               if ( i == 3 ) {
	    
                  // This builder continues with neighbor 3
                  ids.push_back(nextNeighbor);
                  path = myPath;
	    
               } else {
	    
                  // Spawn new builders for neighbors 0,1,2
                  buildPencilsWithNeighbors(grid,pencils,id,ids,dimension,myPath);
	    
               }
	  
            }
	
         }

      } else {
         if(debug) {
            std::cout << "I am cell " << id << ". ";
            std::cout << " I am on refinement level 0." << std::endl;
         }
      }// Closes if (refLvl == 0)

      // If we found a neighbor, add it to the list of ids for this pencil.
      if(nextNeighbor != INVALID_CELLID) {
         if (debug) {
            std::cout << " Next neighbor is " << nextNeighbor << "." << std::endl;
         }
         ids.push_back(nextNeighbor);
      }

      // Move to the next cell.
      id = nextNeighbor;
    
   } // Closes while loop

   // Get the x,y - coordinates of the pencil (in the direction perpendicular to the pencil)
   const auto coordinates = grid.get_center(ids[0]);
   double x,y;
   uint ix,iy,iz
      switch(dimension) {
      case 0: {
         ix = 1;
         iy = 2;
         iz = 0;
         break;
      }
      case 1: {
         ix = 2;
         iy = 0;
         iz = 1;
         break;
      }
      case 2: {
         ix = 0;
         iy = 1;
         iz = 2;
         break;
      }
      default: {
         ix = 0;
         iy = 1;
         iz = 2;
         break;
      }
      }

   x = coordinates[ix];
   y = coordinates[iy];
  
   pencils.addPencil(ids,x,y);
   return pencils;
  
}

void propagatePencil(Vec dr[], Vec values[], Vec z_translation, uint blocks_per_dim ) {

   // Determine direction of translation
   // part of density goes here (cell index change along spatial direcion)
   const int target_scell_index = (z_translation > 0) ? 1: -1; 

   // Vector buffer where we write data, initialized to 0*/
   Vec targetValues[(blocks_per_dim + 2) * WID];
  
   for (uint k_block = 0; k_block < blocks_per_dim; k_block++){

      for (uint k_cell=0; k_cell < WID; ++k_cell) {

         uint gid = k_block * WID + k_cell + WID;
         // init target_values
         targetValues[gid] = 0.0;
      
      }
   }
   for (uint k_block = 0; k_block < blocks_per_dim; k_block++){
    
      for (uint k_cell=0; k_cell < WID; ++k_cell){

         uint gid = k_block * WID + k_cell + WID;
         //uint gid = (blocks_per_dim + 2) * WID - (k_block * WID + k_cell + WID);
      
         // Calculate normalized coordinates in current cell.
         // The coordinates (scaled units from 0 to 1) between which we will
         // integrate to put mass in the target  neighboring cell.
         // Normalize the coordinates to the origin cell. Then we scale with the difference
         // in volume between target and origin later when adding the integrated value.
         Realv z_1,z_2;
         if ( z_translation < 0 ) {
            z_1 = 0;
            z_2 = -z_translation / dr[gid]; 
         } else {
            z_1 = 1.0 - z_translation / dr[gid]; 
            z_2 = 1.0;
         }

         if( abs(z_1) > 1.0 || abs(z_2) > 1.0 ) {
            std::cout << "Error, CFL condition violated\n";
            std::cout << "Exiting\n";
            std::exit(1);
         }
      
         // Compute polynomial coefficients
         Vec a[3];
         //compute_ppm_coeff_nonuniform(dr, values, h4, gid + target_scell_index, a);
         compute_ppm_coeff_nonuniform(dr, values, h4, gid, a);

         // Compute integral
         const Vec ngbr_target_density =
            z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
            z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );

         // Store mapped density in two target cells
         // in the neighbor cell we will put this density        
         targetValues[gid + target_scell_index] +=  ngbr_target_density * dr[gid] / dr[gid + target_scell_index];
         // in the current original cells we will put the rest of the original density
         targetValues[gid]                      +=  values[gid] - ngbr_target_density;
      }
   }

   // Store target data into source data
   for (uint k_block = 0; k_block<blocks_per_dim;k_block++){
    
      for (uint k_cell=0; k_cell<WID; ++k_cell){

         uint gid = k_block * WID + k_cell + WID;
         //uint gid = (blocks_per_dim + 2) * WID - (k_block * WID + k_cell + WID);
         values[gid] = targetValues[gid];
      
      }
    
   }  
  
}


bool trans_map_1d(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const vector<CellID>& localPropagatedCells,
                  const vector<CellID>& remoteTargetCells,
                  const uint dimension,
                  const Realv dt,
                  const uint popID) {
   
   vector<Realv> dz; /*< cell size in the dimension of the pencil */
   Realv dvz,vz_min;  
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned char  cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/

   // return if there's no cells to propagate
   if(localPropagatedCells.size() == 0) 
      return true;
  
   // Vector with all cell ids
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());
  
   const uint nSourceNeighborsPerCell = 1 + 2 * VLASOV_STENCIL_WIDTH;

   // Vectors of pointers to the cell structs
   std::vector<SpatialCell*> allCellsPointer(allCells.size());
   std::vector<SpatialCell*> sourceNeighbors(localPropagatedCells.size() * nSourceNeighborsPerCell);
   std::vector<SpatialCell*> targetNeighbors(3 * localPropagatedCells.size() );
  
   // Initialize allCellsPointer
#pragma omp parallel for
   for(uint celli = 0; celli < allCells.size(); celli++){         
      allCellsPointer[celli] = mpiGrid[allCells[celli]];
   }
   
   // ****************************************************************************
   
   // compute pencils => set of pencils (shared datastructure)
   vector<CellID> seedIds;  

#pragma omp parallel for
   for(auto celli: localPropagatedCells){             
      // Collect a list of cell ids that do not have a neighbor in the negative direction
      // These are the seed ids for the pencils.
      vector<CellID> negativeNeighbors;
      // Returns all neighbors as (id, direction-dimension) pairs.
      for ( const auto neighbor : grid.get_face_neighbors_of(allCellsPointer[celli]) ) {      

         // select the neighbor in the negative dimension of the propagation
         if (neighbor.second == - (dimension + 1))

            // add the id of the neighbor to a list
            negativeNeighbors.push_back(neighbor.first);
      }
      // if no neighbors were found in the negative direction, add this cell id to the seed cells
      if (negativeNeighbors.size() == 0)
         seedIds.push_back(celli);    
   }

   // Empty vectors for internal use of buildPencilsWithNeighbors. Could be default values but
   // default vectors are complicated. Should overload buildPencilsWithNeighbors like suggested here
   // https://stackoverflow.com/questions/3147274/c-default-argument-for-vectorint
   vector<CellID> ids;
   vector<uint> path;

   // Output vectors for ready pencils
   setOfPencils pencils;
   vector<setOfPencils> pencilSets;
   
   for (const auto seedId : seedIds) {
      // Construct pencils from the seedIds into a set of pencils.
      pencils = buildPencilsWithNeighbors(mpiGrid, pencils, seedId, ids, dimension, path);
   }
   // Add the final set of pencils to the pencilSets - vector.
   // Only one set is created for now but we retain support for multiple sets
   pencilSets.push_back(pencils);
   // ****************************************************************************
  
   // Fiddle indices x,y,z
   switch (dimension) {
   case 0:
      // set values in array that is used to convert block indices 
      // to global ID using a dot product.
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
   case 1:
      // set values in array that is used to convert block indices 
      // to global ID using a dot product
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
   case 2:
      // set values in array that is used to convert block indices
      // to global id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   default:
      cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
      abort();
      break;
   }
           
   // init cellid_transpose
   for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
         for (uint i=0; i<WID; ++i) {
            const uint cell =
               i * cell_indices_to_id[0] +
               j * cell_indices_to_id[1] +
               k * cell_indices_to_id[2];
            cellid_transpose[ i + j * WID + k * WID2] = cell;
         }
      }
   }
      // ****************************************************************************

      const uint8_t VMESH_REFLEVEL = 0;
      
      // Get a pointer to the velocity mesh of the first spatial cell
      const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = allCellsPointer[0]->get_velocity_mesh(popID);

      // set cell size in dimension direction
      dvz = vmesh.getCellSize(VMESH_REFLEVEL)[dimension];
      vz_min = vmesh.getMeshMinLimits()[dimension];
  
      // Get a unique sorted list of blockids that are in any of the
      // propagated cells. First use set for this, then add to vector (may not
      // be the most nice way to do this and in any case we could do it along
      // dimension for data locality reasons => copy acc map column code, TODO: FIXME
      // TODO: Do this separately for each pencil?
      std::unordered_set<vmesh::GlobalID> unionOfBlocksSet;    
  
      for(auto cell : allCellsPointer) {
         vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
         for (vmesh::LocalID block_i=0; block_i< vmesh.size(); ++block_i) {
            unionOfBlocksSet.insert(vmesh.getGlobalID(block_i));
         }
      }
  
      std::vector<vmesh::GlobalID> unionOfBlocks;
      unionOfBlocks.reserve(unionOfBlocksSet.size());
      for(const auto blockGID:  unionOfBlocksSet) {
         unionOfBlocks.push_back(blockGID);
      }
      // ****************************************************************************
  
      int t1 = phiprof::initializeTimer("mappingAndStore");
      
#pragma omp parallel 
      {      
         //std::vector<bool> targetsValid(localPropagatedCells.size());
         //std::vector<vmesh::LocalID> allCellsBlockLocalID(allCells.size());
             
#pragma omp for schedule(guided)
         // Loop over velocity space blocks. Thread this loop (over vspace blocks) with OpenMP.    
         for(uint blocki = 0; blocki < unionOfBlocks.size(); blocki++){
            
            phiprof::start(t1);

            // Get global id of the velocity block
            vmesh::GlobalID blockGID = unionOfBlocks[blocki];

            velocity_block_indices_t block_indices;
            uint8_t vRefLevel;
            vmesh.getIndices(blockGID,vRefLevel, block_indices[0],
                             block_indices[1], block_indices[2]);      
      
            // Loop over sets of pencils
            // This loop only has one iteration for now
            for ( auto pencils: pencilSets ) {

               // Allocate targetdata sum(lengths of pencils)*WID3)
               Vec targetData[pencils.sumOfLengths * WID3];

               // Initialize targetdata to 0
               for( uint i = 0; i<targetData.size(); i++ ) {
                  targetData[i] = Vec(0.0);
               }

               // TODO: There's probably a smarter way to keep track of where we are writing
               //       in the target data structure.
               uint targetDataIndex = 0;

               // Compute spatial neighbors for target cells.
               // For targets we only have actual cells as we do not
               // want to propagate boundary cells (array may contain
               // INVALID_CELLIDs at boundaries).
               for ( auto celli: pencils.ids ) {
                  compute_spatial_target_neighbors(mpiGrid, localPropagatedCells[celli], dimension,
                                                   targetNeighbors.data() + celli * 3);
               }
	
               // Loop over pencils	
               for(uint pencili = 0; pencili < pencils.N; pencili++){

                  // Allocate source data: sourcedata<length of pencil * WID3)
                  Vec sourceData[pencils.lengthOfPencils[pencili] * WID3];

                  // Compute spatial neighbors for source cells. In
                  // source cells we have a wider stencil and take into account
                  // boundaries. 
                  for( auto celli: pencils.getIds(pencili)) {
                     compute_spatial_source_neighbors(mpiGrid, localPropagatedCells[celli],
                                                      dimension, sourceNeighbors.data()
                                                      + celli * nSourceNeighborsPerCell);

                     // load data(=> sourcedata) / (proper xy reconstruction in future)
                     // copied from regular code, should work?

                     // TODO: Does the index of sourceData need adjustments for vector length?
                     copy_trans_block_data(sourceNeighbors.data() + celli * nSourceNeighborsPerCell,
                                           blockGID, sourceData[celli], cellid_transpose, popID);
	    
                     // At the same time, calculate dz's and store them in a vector.
                     dz.push_back(dz_ini / 2.0 ** mpiGrid.get_refinement_level(celli));

                  }	 
	    	  	  
                  // Calculate cell centered velocity
                  const Realv cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min;
	  
                  const Realv z_translation = dt * cell_vz / dz[celli];
                  // propagate pencil(blockid = velocities, pencil-ids = dzs ),
                  propagatePencil(dz, sourceData, z_translation, blocks_per_dim);

                  // sourcedata => targetdata[this pencil])
                  for (auto value: sourceData) {
                     targetData[targetDataIndex] = value;
                     targetDataindex++;
                  }
                  
                  // dealloc source data -- Should be automatic since it's declared in this iteration?
	  
               }
      
               // Loop over pencils again
               for(uint pencili = 0; pencili < pencils.N; pencili++){

                  // store_data(target_data =>)  :Aggregate data for blockid to original location 
                  
                  //store values from target_values array to the actual blocks
                  for(auto celli: pencils.ids) {
                     //TODO: Figure out validity check later
                     //if(targetsValid[celli]) {
                        for(uint ti = 0; ti < 3; ti++) {
                           SpatialCell* spatial_cell = targetNeighbors[celli * 3 + ti];
                           if(spatial_cell ==NULL) {
                              //invalid target spatial cell
                              continue;
                           }
                           
                           // Get local ID of the velocity block
                           const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blockGID, popID);
                           
                           if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                              // block does not exist. If so, we do not create it and add stuff to it here.
                              // We have already created blocks around blocks with content in
                              // spatial sense, so we have no need to create even more blocks here
                              // TODO add loss counter
                              continue;
                           }
                           // Pointer to the data field of the velocity block
                           Realf* blockData = spatial_cell->get_data(blockLID, popID);
                           for(int i = 0; i < WID3 ; i++) {

                              // Write data into target block
                              blockData[i] += targetData[(celli * 3 + ti) * WID3 + i];
                           }
                        }
                        //}
                     
                  }
                  
                  // dealloc target data -- Should be automatic again?
               }
            }
         }
      }

      return true;
   }
