#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <mpi.h>
#include <vector>
#include <string>
#include <limits>

#include "definitions.h"

const uint MAX_SPAT_CELLS = 8100;
const uint MAX_VEL_BLOCKS = 250000;
//const uint MAX_VEL_BLOCKS = 2500000;

const uint CUDA_WIDTH = 65536; // Width of CUDA array (for 2D textures)
const uint CUDA_HEIGHT = 3000; // Height of CUDA array
                          // Make sure width*height / 64 >= MAX_VEL_BLOCKS

#ifdef PARGRID
   const uint INVALID_CELLID = std::numeric_limits<uint>::max();
#else
   const uint64_t INVALID_CELLID = 0;
#endif

namespace Transmit {
   const uint CELL_PARAMS  = 1;
   const uint BLOCK_PARAMS = CELL_PARAMS << 1;
   const uint AVGS         = CELL_PARAMS << 2;
   const uint FLUXES       = CELL_PARAMS << 3;
   const uint DERIV1       = CELL_PARAMS << 4;
   const uint DERIV2       = CELL_PARAMS << 5;
   const uint NBRSVEL      = CELL_PARAMS << 6;
}

struct Parameters {
   static Real xmin;  /**< X-coordinate of the lower left corner of the spatial grid. */
   //static Real xmax;  /**< Y-coordinate of the lower left corner of the spatial grid. */
   static Real ymin;  /**< Z-coordinate of the lower left corner of the spatial grid. */
   //static Real ymax;  /**< X-coordinate of the upper right corner of the spatial grid. */
   static Real zmin;  /**< Y-coordinate of the upper right corner of the spatial grid. */
   //static Real zmax;  /**< Z-coordinate of the upper right corner of the spatial grid. */
   static Real dx_ini; /**< Initial size of spatial cell in x-direction. */
   static Real dy_ini; /**< Initial size of spatial cell in y-direction. */
   static Real dz_ini; /**< Initial size of spatial cell in z-direction. */
   
   static Real vxmin; /**< VX-coordinate of the lower left corner of velocity grid. */
   static Real vxmax; /**< VY-coordinate of the lower left corner of velocity grid. */
   static Real vymin; /**< VZ-coordinate of the lower left corner of velocity grid. */
   static Real vymax; /**< VX-coordinate of the upper right corner of velocity grid. */
   static Real vzmin; /**< VY-coordinate of the upper right corner of velocity grid. */
   static Real vzmax; /**< VZ-coordinate of the upper right corner of velocity grid. */
   
   static uint xcells_ini; /**< Initial number of spatial cells in x-direction. */
   static uint ycells_ini; /**< Initial number of spatial cells in y-direction. */
   static uint zcells_ini; /**< Initial number of spatial cells in z-direction. */
   static uint vxblocks_ini; /**< Initial number of velocity grid blocks in vx-direction. */
   static uint vyblocks_ini; /**< Initial number of velocity grid blocks in vy-direction. */
   static uint vzblocks_ini; /**< Initial number of velocity grid blocks in vz-direction. */

   static Real q;                    /**< Charge of simulated particle species.*/
   static Real m;                    /**< Mass of simulated particle species.*/
   static Real q_per_m;              /**< Charge-to-mass ratio of simulated particle species,
				      * calculated from Parameters::q and Parameters::m.*/
   static Real t;                    /**< Current simulation time. */
   static Real dt;                   /**< The value of the timestep to use in propagation. */
   static luint tstep_min;           /**< Timestep when simulation starts, needed for restarts.*/
   static luint tstep;               /**< The number of the current timestep. 0=initial state. */
   static luint tsteps;              /**< Total number of timesteps to calculate. */
   static luint saveRestartInterval;
   static luint diagnInterval;

   static std::string solar_wind_file;	/**< Read solar wind data from this file. */

   static bool save_spatial_grid;	/**< Save spatial cell averages for the whole simulation. */
   static bool save_velocity_grid;	/**< Save the velocity grid of every spatial cell in the simulation. */
   static std::vector<Real> save_spatial_cells_x;	/**< Save the velocity grid of spatial cells at these locations. */
   static std::vector<Real> save_spatial_cells_y;
   static std::vector<Real> save_spatial_cells_z;
   
   static uint transmit; /**< Indicates the data that needs to be transmitted to remote nodes.
			  * This is created with bitwise or from the values defined in 
			  * namespace Transmit.*/

   static bool recalculateStencils; /**< If true, MPI stencils should be recalculated because of 
				     * load balancing.*/
   
   static bool propagateField;      /**< If true, magnetic field is propagated during the simulation.*/
   static bool propagateVlasov;     /**< If true, distribution function is propagated during the simulation.*/
   static bool periodic_x, periodic_y, periodic_z; /**< Whether spatial vlasov grid is periodic */

};

struct Readparameters {
    Readparameters(int argc, char* argv[],MPI_Comm comm);
    static bool add(const std::string& name,const std::string& desc,const std::string& defValue);
    static bool add(const std::string& name,const std::string& desc,const bool& defValue);
    static bool add(const std::string& name,const std::string& desc,const int& defValue);
    static bool add(const std::string& name,const std::string& desc,const unsigned int& defValue);
    static bool add(const std::string& name,const std::string& desc,const float& defValue);
    static bool add(const std::string& name,const std::string& desc,const double& defValue);

    static bool get(const std::string& name,std::string& value);
    static bool get(const std::string& name,bool& value);
    static bool get(const std::string& name,int& value);
    static bool get(const std::string& name,unsigned int& value);
    static bool get(const std::string& name,unsigned long& value);
    static bool get(const std::string& name,float& value);
    static bool get(const std::string& name,double& value);

//Functions for composing options (can be defined multiple times and are all returned as a vector)
    static bool addComposing(const std::string& name,const std::string& desc);
    static bool get(const std::string& name,std::vector<std::string>& value);
    static bool get(const std::string& name,std::vector<int>& value);
    static bool get(const std::string& name,std::vector<float>& value);
    static bool get(const std::string& name,std::vector<double>& value);

    
    static bool finalize();
    static bool helpMessage();
    static bool isInitialized();
    static bool parse();
   
private:
    static int argc;                  /**< How many entries argv contains.*/
    static char** argv;              /**< Pointer to char* array containing command line parameters.*/
    static int rank;
    static MPI_Comm comm;

  
    /** Private default constructor to prevent incorrect initialization.*/
    Readparameters();
    static bool addDefaultParameters();
};


#endif
