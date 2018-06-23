bool trans_map_1d_amr(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
		      const vector<CellID>& localPropagatedCells,
		      const vector<CellID>& remoteTargetCells,
		      const uint dimension,
		      const Realv dt,
		      const uint popID) {


  // values used with an stencil in 1 dimension, initialized to 0. 
  // Contains a block, and its spatial neighbours in one dimension.
  vector<Realv> dz;
  Realv z_min, dvz,vz_min;  
  uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
  unsigned char  cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/
  
  if(localPropagatedCells.size() == 0) 
    return true; 
  //vector with all cells
  vector<CellID> allCells(localPropagatedCells);
  allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());
  
  const uint nSourceNeighborsPerCell = 1 + 2 * VLASOV_STENCIL_WIDTH;
  std::vector<SpatialCell*> allCellsPointer(allCells.size());
  std::vector<SpatialCell*> sourceNeighbors(localPropagatedCells.size() * nSourceNeighborsPerCell);
  std::vector<SpatialCell*> targetNeighbors(3 * localPropagatedCells.size() );
  
  
#pragma omp parallel for
  for(uint celli = 0; celli < allCells.size(); celli++){         
    allCellsPointer[celli] = mpiGrid[allCells[celli]];
  }
  

  vector<CellID> seedIds;  

#pragma omp parallel for
  for(uint celli = 0; celli < localPropagatedCells.size(); celli++){         
    // compute spatial neighbors, separately for targets and source. In
    // source cells we have a wider stencil and take into account
    // boundaries. For targets we only have actual cells as we do not
    // want to propagate boundary cells (array may contain
    // INVALID_CELLIDs at boundaries).
    compute_spatial_source_neighbors(mpiGrid, localPropagatedCells[celli], dimension,
				     sourceNeighbors.data() + celli * nSourceNeighborsPerCell);
    compute_spatial_target_neighbors(mpiGrid, localPropagatedCells[celli], dimension,
				     targetNeighbors.data() + celli * 3);
    
    // Collect a list of cell ids that do not have a neighbor in the negative direction
    // These are the seed ids for the pencils.
    vector<CellID> negativeNeighbors;
    for (const auto neighbor : grid.get_face_neighbors_of(celli)) {      
      
      if (neighbor.second == - (dimension + 1))
	negativeNeighbors.push_back(neighbor.first);
    }
    if (negativeNeighbors.size() == 0)
      seedIds.push_back(celli);    
  }

  // compute pencils => set of pencils (shared datastructure)
  vector<CellID> ids;
  vector<uint> path;  
  setOfPencils pencils;

  for (const auto id : seedIds) {
    pencils = buildPencilsWithNeighbors(mpiGrid, pencils, id, ids, dimension, path);
  }
  
  // Get a unique sorted list of blockids that are in any of the
  // propagated cells. First use set for this, then add to vector (may not
  // be the most nice way to do this and in any case we could do it along
  // dimension for data locality reasons => copy acc map column code, TODO: FIXME
  std::unordered_set<vmesh::GlobalID> unionOfBlocksSet;    
  
  for(uint celli = 0; celli < allCellsPointer.size(); celli++) {
    vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = allCellsPointer[celli]->get_velocity_mesh(popID);
    for (vmesh::LocalID block_i=0; block_i< vmesh.size(); ++block_i) {
      unionOfBlocksSet.insert(vmesh.getGlobalID(block_i));
    }
  }
  
  std::vector<vmesh::GlobalID> unionOfBlocks;
  unionOfBlocks.reserve(unionOfBlocksSet.size());
  for(const auto blockGID:  unionOfBlocksSet) {
    unionOfBlocks.push_back(blockGID);
  }
  
  // Fiddle indices x,y,z

  const uint8_t VMESH_REFLEVEL=0;
  const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = allCellsPointer[0]->get_velocity_mesh(popID);
  // set cell size in dimension direction
  dvz = vmesh.getCellSize(VMESH_REFLEVEL)[dimension];
  vz_min = vmesh.getMeshMinLimits()[dimension];
  switch (dimension) {
  case 0:
    // TODO: Query the pencil for the dz's? Do this inside the loop over pencils
    dz = ;
    z_min = P::xmin;      
    // set values in array that is used to convert block indices 
    // to global ID using a dot product.
    cell_indices_to_id[0]=WID2;
    cell_indices_to_id[1]=WID;
    cell_indices_to_id[2]=1;
    break;
  case 1:
    z_min = P::ymin;
    // set values in array that is used to convert block indices 
    // to global ID using a dot product
    cell_indices_to_id[0]=1;
    cell_indices_to_id[1]=WID2;
    cell_indices_to_id[2]=WID;
    break;
  case 2:
    z_min = P::zmin;
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
           
  // init plane_index_to_id
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
  
  int t1 = phiprof::initializeTimer("mapping");
  int t2 = phiprof::initializeTimer("store");

#pragma omp parallel 
  {      
    std::vector<Realf> targetBlockData(3 * localPropagatedCells.size() * WID3);
    std::vector<bool> targetsValid(localPropagatedCells.size());
    std::vector<vmesh::LocalID> allCellsBlockLocalID(allCells.size());
             
#pragma omp for schedule(guided)
    // Loop over unionOfBlocks
    for(uint blocki = 0; blocki < unionOfBlocks.size(); blocki++){
      vmesh::GlobalID blockGID = unionOfBlocks[blocki];
      phiprof::start(t1);

      // Loop over sets of pencils

      // Loop over pencils
      
    }
  }
}
