#include <einspline/bspline.h>

//template<typename T> inline T cell_ib_id(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> inline T cell_ib_id(const T& i,const T& j,const T& k) {return k*36+j*6+i;}
template<typename T> inline T cell_id(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}

class interpolated_block {

  public:
   
   interpolated_block () {}
   
   inline double get_value(double x,double y,double z){
      double value;
      eval_UBspline_3d_d(spline,x,y,z,&value);
      return value
   }

   void set_block(Velocity_Block*  block){
      block_ptr=block;
      this->load_data();
      this->prepare_einspline();
   }
   
   
private:
   
   void prepare_einspline(){
      /* Here we set the grid, used for the computation */
      Ugrid x_grid;
      Ugrid y_grid;
      Ugrid z_grid;
      x_grid.start = block_ptr->parameters[BlockParams::VXCRD] - 0.5*block_ptr->parameters[BlockParams::DVX];
      x_grid.end = block_ptr->parameters[BlockParams::VXCRD] + (WID+0.5)*block_ptr->parameters[BlockParams::DVX];
      x_grid.num = WID+2; //< how long is our array of points ?
      y_grid.start = block_ptr->parameters[BlockParams::VYCRD] - 0.5*block_ptr->parameters[BlockParams::DVY];
      y_grid.end = block_ptr->parameters[BlockParams::VYCRD] + (WID+0.5)*block_ptr->parameters[BlockParams::DVY];
      y_grid.num = WID+2; //< how long is our array of points ? 
      z_grid.start = block_ptr->parameters[BlockParams::VZCRD] - 0.5*block_ptr->parameters[BlockParams::DVZ];
      z_grid.end = block_ptr->parameters[BlockParams::VZCRD] + (WID+0.5)*block_ptr->parameters[BlockParams::DVZ];
      z_grid.num = WID+2; //< how long is our array of points ?


      /** Boundary Conditions **/
      BCtype_d xBC = {NATURAL, NATURAL , 0.,0.};
      BCtype_d yBC = {NATURAL, NATURAL , 0.,0.};
      BCtype_d zBC = {NATURAL, NATURAL , 0.,0.};

      spline=create_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,avgs);


      
   }
   


   
   
   load_data(){
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
   
      // ***** NEIGHBOURS TO NEGATIVE VZ DIRECTION *****           
      // Neighbour (iv-1,jv-1,kv-1):
   
      if (block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZM1] != NULL){
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XM1_YM1_ZM1]->fx;
         avgs[cell_ib_id(MIN,MIN,MIN)]=nbrAvgs[cell_id(BLOCKMAX,BLOCKMAX,BLOCKMAX)];
      }
      else {
         avgs[cell_ib_id(MIN,MIN,MIN)]=0.0;
      }
   
   
      // Neighbour (iv  ,jv-1,kv-1):
      if (block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZM1] != NULL) {
         nbrAvgs = block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZM1]->fx;
         for (unsigned int i=0; i<WID; ++i)
            avgs[cell_ib_id(i+1,MIN,MIN)]=nbrAvgs[cell_id(i,BLOCKMAX,BLOCKMAX)];
      }
      else{
         for (unsigned int i=0; i<WID; ++i)
            avgs[cell_ib_id(i+1,MIN,MIN)]=0.0;
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
   
   Velocity_Block*  block_ptr;
   double avgs[512];
   UBspline_3d_d *spline;
}

