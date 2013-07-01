


//template<typename T> inline T cell_ib_id(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> inline T cell_ib_id(const T& i,const T& j,const T& k) {return k+j*6+i*36;} //internal get ib cell id. This takes in indices in natural 6x6x6 frame
template<typename T> inline T cell_id(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}

enum BlockInterpolationType { CONSTANT, HINGED_HYPERPLANE};
enum HingedHyperplaneParams {X0,Y0,Z0,DFDX,DFDY,DFDZ,NUM_HH_PARAMS};

class interpolated_block {

  public:
   
   interpolated_block (BlockInterpolationType type) : interpolationType(type) {}
   
   inline double get_value(double x,double y,double z){
      switch(this->interpolationType) {
          case CONSTANT:
             return eval_nointerpolation(x,y,z);
          case HINGED_HYPERPLANE:
             return eval_hinged_hyperplane(x,y,z);

      }
      return 0.0;
   }
   
   inline double get_value(unsigned int cellid,double x,double y,double z){
      switch(this->interpolationType) {
          case CONSTANT:
             return eval_nointerpolation(x,y,z);
          case HINGED_HYPERPLANE:
             return eval_hinged_hyperplane(cellid,x,y,z);
      }
      return 0.0;
   }

   void load_block(Velocity_Block* block){
      block_ptr = block;
      this->load_data();
      switch(this->interpolationType) {
          case CONSTANT:
             break; 
          case HINGED_HYPERPLANE:     
             this->prepare_hinged_hyperplane();
             break;
      }
   }
   inline unsigned int get_cell_id(const unsigned int i,const unsigned int j,const unsigned int k){
      return cell_ib_id(i+1,j+1,k+1);
   }

private:
   /*private parameters*/
   const BlockInterpolationType interpolationType;
   Velocity_Block*  block_ptr;
   //TODO FIX 216 now hardcoded
   double avgs[216];
   //TODO, extra stuff stored for all, at least hh_parameters should perhaps be dynamic
   double hh_parameters[216][NUM_HH_PARAMS]; /*< Parameters for hinged hyperplanes*/

   /*private functions*/
   double eval_nointerpolation(double x,double y,double z) {
      const unsigned int cell_x=(x-block_ptr->parameters[BlockParams::VXCRD])/block_ptr->parameters[BlockParams::DVX];
      const unsigned int cell_y=(y-block_ptr->parameters[BlockParams::VYCRD])/block_ptr->parameters[BlockParams::DVY];
      const unsigned int cell_z=(z-block_ptr->parameters[BlockParams::VZCRD])/block_ptr->parameters[BlockParams::DVZ];
      const unsigned int cell_id=cell_ib_id(cell_x+1,cell_y+1,cell_z+1);
      return avgs[cell_id];
   }
   
   double eval_hinged_hyperplane(double x,double y,double z) {
      //TODO, we could store x0,y0,z0 and dfdxyz for each cell when reading the block and then reuse them when evaluating. Coudl save for multiple subcells.
      //TODO, derivatives need to be limited to avoid negative regions. Do we also need similar limiters like in fvm...?
      const unsigned int cell_x=(x-block_ptr->parameters[BlockParams::VXCRD])/block_ptr->parameters[BlockParams::DVX];
      const unsigned int cell_y=(y-block_ptr->parameters[BlockParams::VYCRD])/block_ptr->parameters[BlockParams::DVY];
      const unsigned int cell_z=(z-block_ptr->parameters[BlockParams::VZCRD])/block_ptr->parameters[BlockParams::DVZ];
      const unsigned int cell_id=cell_ib_id(cell_x+1,cell_y+1,cell_z+1);
      return avgs[cell_id]+
         (x - hh_parameters[cell_id][X0]) * hh_parameters[cell_id][DFDX]+
         (y - hh_parameters[cell_id][Y0]) * hh_parameters[cell_id][DFDY]+
         (z - hh_parameters[cell_id][Z0]) * hh_parameters[cell_id][DFDZ];
   }
   /*get value in a particular velocity cell, at x,y,z*/
   double eval_hinged_hyperplane(unsigned int cell_id,double x,double y,double z){
      //TODO, we could store x0,y0,z0 and dfdxyz for each cell when reading the block and then reuse them when evaluating. Coudl save for multiple subcells.
      //TODO, derivatives need to be limited to avoid negative regions. Do we also need similar limiters like in fvm...?
      return avgs[cell_id]+
         (x - hh_parameters[cell_id][X0]) * hh_parameters[cell_id][DFDX]+
         (y - hh_parameters[cell_id][Y0]) * hh_parameters[cell_id][DFDY]+
         (z - hh_parameters[cell_id][Z0]) * hh_parameters[cell_id][DFDZ];
   }
   

   void prepare_hinged_hyperplane(){
      for(unsigned int cell_x=0;cell_x<WID;cell_x++){
         for(unsigned int cell_y=0;cell_y<WID;cell_y++){
            for(unsigned int cell_z=0;cell_z<WID;cell_z++){
               const unsigned int cell_id=cell_ib_id(cell_x+1,cell_y+1,cell_z+1);
               hh_parameters[cell_id][X0]=block_ptr->parameters[BlockParams::VXCRD]+(cell_x+0.5)*block_ptr->parameters[BlockParams::DVX];
               hh_parameters[cell_id][Y0]=block_ptr->parameters[BlockParams::VYCRD]+(cell_y+0.5)*block_ptr->parameters[BlockParams::DVY];
               hh_parameters[cell_id][Z0]=block_ptr->parameters[BlockParams::VZCRD]+(cell_z+0.5)*block_ptr->parameters[BlockParams::DVZ];

               if(avgs[cell_id]<0) {
                  //Not-so-nice situation of negative dist function produced by limiters of Leveque, or something else*/
                  //Let's  just have a flat hyperplane here to avoid more serious problems
                  hh_parameters[cell_id][DFDX]=0.0;
                  hh_parameters[cell_id][DFDY]=0.0;
                  hh_parameters[cell_id][DFDZ]=0.0;
                  continue;
               }
               hh_parameters[cell_id][DFDX]=(avgs[cell_ib_id(cell_x+2,cell_y+1,cell_z+1)]-
                                  avgs[cell_ib_id(cell_x+0,cell_y+1,cell_z+1)])/(2.0*block_ptr->parameters[BlockParams::DVX]);
               hh_parameters[cell_id][DFDY]=(avgs[cell_ib_id(cell_x+1,cell_y+2,cell_z+1)]-
                                             avgs[cell_ib_id(cell_x+1,cell_y+0,cell_z+1)])/(2.0*block_ptr->parameters[BlockParams::DVY]);
               hh_parameters[cell_id][DFDZ]=(avgs[cell_ib_id(cell_x+1,cell_y+1,cell_z+2)]-
                                  avgs[cell_ib_id(cell_x+1,cell_y+1,cell_z+0)])/(2.0*block_ptr->parameters[BlockParams::DVZ]);
               
               /*Get the minimum value of the interpolation from one
                * of the 8 corners of the cell. We know based on the
                * derivatives which it is*/
               double minVal=eval_hinged_hyperplane(cell_id,
                                                    block_ptr->parameters[BlockParams::VXCRD]+
                                                    (cell_x+(hh_parameters[cell_id][DFDX]<0))*block_ptr->parameters[BlockParams::DVX],
                                                    block_ptr->parameters[BlockParams::VYCRD]+
                                                    (cell_y+(hh_parameters[cell_id][DFDY]<0))*block_ptr->parameters[BlockParams::DVY],
                                                    block_ptr->parameters[BlockParams::VZCRD]+
                                                    (cell_z+(hh_parameters[cell_id][DFDZ]<0))*block_ptr->parameters[BlockParams::DVZ]);
               if(minVal<0){
                  /*scale derivatives so that the minimum value of the interpolation is zero*/
                  double alpha=avgs[cell_id]/(avgs[cell_id]-minVal);
                  hh_parameters[cell_id][DFDX]*=alpha;
                  hh_parameters[cell_id][DFDY]*=alpha;
                  hh_parameters[cell_id][DFDZ]*=alpha;
               }
               
               
            } 
         } 
      } 
   }
   


   void load_data(){
      Real* nbrAvgs;
      //typedef Parameters P;
      
      // NOTE: velocity block can have up to 26 neighbours, and we construct
      // a 6x6x6 sized block with this blocks values, and the stencil 1 values from neighboring blocs.
      //If block does not exist, its values are 0.0.
      //we copy from fx as we assume dist function has been copied there in SL algorithm.
      // NOTE2: avgs is a (6,6,6) sized block, nbrAvgs is (4,4,4).
  
      const unsigned int MIN = 0;
      const unsigned int MAX = 5;
      const unsigned int BLOCKMIN = 0;
      const unsigned int BLOCKMAX = 3;

      for (unsigned int i=0; i<216;i++)
         avgs[i]=0.0;

      
      // ***** NEIGHBOURS TO NEGATIVE VZ DIRECTION *****           
      // Neighbour (iv-1,jv-1,kv-1):
   
      if (block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZM1] != NULL){
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZM1]->fx;
         avgs[cell_ib_id(MIN,MIN,MIN)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMAX,BLOCKMAX)];
      }
   
      // Neighbour (iv  ,jv-1,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZM1] != NULL) {
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZM1]->fx;
         for (unsigned int i=0; i<WID; ++i)
            avgs[cell_ib_id(i+1,MIN,MIN)]=nbrAvgs[cell_id(i,BLOCKMAX,BLOCKMAX)];
      }
      
      // Neighbour (iv+1,jv-1,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YM1_ZM1] != NULL) {
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XP1_YM1_ZM1]->fx;
         avgs[cell_ib_id(MAX,MIN,MIN)]=nbrAvgs[cell_id(BLOCKMIN,BLOCKMAX,BLOCKMAX)]; 
      }
   
      // Neighbour (iv-1,jv  ,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZM1] != NULL) {
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZM1]->fx;
         for (unsigned int j=0; j<WID; ++j)
            avgs[cell_ib_id(MIN,j+1,MIN)] = nbrAvgs[cell_id(BLOCKMAX,j,BLOCKMAX)];         
      }

      
      // Neighbour (iv  ,jv  ,kv-1):       
      if (block_ptr->neighbors[velocity_neighbor::XCC_YCC_ZM1] != NULL) {
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XCC_YCC_ZM1]->fx;
         for (unsigned int j=0; j<WID; ++j)
            for (unsigned int i=0; i<WID; ++i)
               avgs[cell_ib_id(i+1,j+1,MIN)]=nbrAvgs[cell_id(i,j,BLOCKMAX)];
      }

      // Neighbour (iv+1,jv  ,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZM1] != NULL) {
         nbrAvgs  = block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZM1]->fx;
         for (unsigned int j=0; j<WID; ++j)
            avgs[cell_ib_id(MAX,j+1,MIN)]=nbrAvgs[cell_id(BLOCKMIN,j,BLOCKMAX)];
      }
         
      // Neighbour (iv-1,jv+1,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YP1_ZM1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YP1_ZM1]->fx;
         avgs[cell_ib_id(MIN,MAX,MIN)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMIN,BLOCKMAX)];
      }
      // Neighbour (iv  ,jv+1,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZM1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZM1]->fx;
         for (unsigned int i=0; i<WID; ++i)
            avgs[cell_ib_id(i+1,MAX,MIN)]=nbrAvgs[cell_id(i,BLOCKMIN,BLOCKMAX)];
      }
      // Neighbour (iv+1,jv+1,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YP1_ZM1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XP1_YP1_ZM1]->fx;
         avgs[cell_ib_id(MAX,MAX,MIN)]=nbrAvgs[cell_id(BLOCKMIN,BLOCKMIN,BLOCKMAX)];
      }
   
      // ***** NEIGHBOURS IN SAME VZ PLANE *****
      // Neighbour (iv-1,jv-1,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            avgs[cell_ib_id(MIN,MIN,k+1)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMAX,k)];
      }
      // Neighbour (iv  ,jv-1,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            for (unsigned int i=0; i<WID; ++i)
               avgs[cell_ib_id(i+1,MIN,k+1)]=nbrAvgs[cell_id(i,BLOCKMAX,k)];
      }
      // Neighbour (iv+1,jv-1,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YM1_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XP1_YM1_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            avgs[cell_ib_id(MAX,MIN,k+1)]=nbrAvgs[cell_id(BLOCKMIN,BLOCKMAX,k)];
      }
      // Neighbour (iv-1,jv  ,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            for (unsigned int j=0; j<WID; ++j)
               avgs[cell_ib_id(MIN,j+1,k+1)]=nbrAvgs[cell_id(BLOCKMAX,j,k)];
      }
   
      // This block (iv  ,jv  ,kv  ):
      for (unsigned int k=0; k<WID; ++k)
         for (unsigned int j=0; j<WID; ++j)
            for (unsigned int i=0; i<WID; ++i) {
               avgs[cell_ib_id(i+1,j+1,k+1)]=block_ptr->fx[cell_id(i,j,k)];
            }
   
      // Neighbour (iv+1,jv  ,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            for (unsigned int j=0; j<WID; ++j)
               avgs[cell_ib_id(MAX,j+1,k+1)]=nbrAvgs[cell_id(BLOCKMIN,j,k)];
      }
      // Neighbour (iv-1,jv+1,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YP1_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YP1_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            avgs[cell_ib_id(MIN,MAX,k+1)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMIN,k)];
      }
      // Neighbour (iv  ,jv+1,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            for (unsigned int i=0; i<WID; ++i)
               avgs[cell_ib_id(i+1,MAX,k+1)]=nbrAvgs[cell_id(i,BLOCKMIN,k)];
      }   
      // Neighbour (iv+1,jv+1,kv  ):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YP1_ZCC] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XP1_YP1_ZCC]->fx;
         for (unsigned int k=0; k<WID; ++k)
            avgs[cell_ib_id(MAX,MAX,k+1)]=nbrAvgs[cell_id(BLOCKMIN,BLOCKMIN,k)];
      }

      // ***** NEIGHBOURS TO POSITIVE VZ DIRECTION *****
      // Neighbour (iv-1,jv-1,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZP1]->fx;
         avgs[cell_ib_id(MIN,MIN,MAX)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMAX,BLOCKMIN)] ;
      }   
      // Neighbour (iv  ,jv-1,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZP1]->fx;
         for (unsigned int i=0; i<WID; ++i)
            avgs[cell_ib_id(i+1,MIN,MAX)]=nbrAvgs[cell_id(i,BLOCKMAX,BLOCKMIN)];
      }
      // Neighbour (iv+1,jv-1,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YM1_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XP1_YM1_ZP1]->fx;
         avgs[cell_ib_id(MAX,MIN,MAX)]=nbrAvgs[cell_id(BLOCKMIN,BLOCKMAX,BLOCKMIN)];
      }
      // Neighbour (iv-1,jv  ,kv+1):       
      if (block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZP1]->fx;
         for (unsigned int j=0; j<WID; ++j)
            avgs[cell_ib_id(MIN,j+1,MAX)]=nbrAvgs[cell_id(BLOCKMAX,j,BLOCKMIN)];
      }
      // Neighbour (iv  ,jv  ,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YCC_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XCC_YCC_ZP1]->fx;
         for (unsigned int j=0; j<WID; ++j)
            for (unsigned int i=0; i<WID; ++i)
               avgs[cell_ib_id(i+1,j+1,MAX)]=nbrAvgs[cell_id(i,j,BLOCKMIN)];
      }
      // Neighbour (iv+1,jv  ,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZP1]->fx;
         for (unsigned int j=0; j<WID; ++j)
            avgs[cell_ib_id(MAX,j+1,MAX)]=nbrAvgs[cell_id(BLOCKMIN,j,BLOCKMIN)];
      }
      // Neighbour (iv-1,jv+1,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XM1_YP1_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XM1_YP1_ZP1]->fx;
         avgs[cell_ib_id(MIN,MAX,MAX)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMIN,BLOCKMIN)];
      }
      // Neighbour (iv  ,jv+1,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZP1] != NULL) {
         nbrAvgs =  block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZP1]->fx;
         for (unsigned int i=0; i<WID; ++i)
            avgs[cell_ib_id(i+1,MAX,MAX)]=nbrAvgs[cell_id(i,BLOCKMIN,BLOCKMIN)];
      }
      // Neighbour (iv+1,jv+1,kv+1):
      if (block_ptr->neighbors[velocity_neighbor::XP1_YP1_ZP1] != NULL) {
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XP1_YP1_ZP1]->fx;
         avgs[cell_ib_id(MAX,MAX,MAX)]=nbrAvgs[cell_id(BLOCKMIN,BLOCKMIN,BLOCKMIN)];
      }
   
   }
};

