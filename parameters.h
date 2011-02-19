#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "definitions.h"
#include <vector>
#include <string>

//cuint MAX_VEL_BLOCKS = 300000;
cuint MAX_VEL_BLOCKS = 300000;
namespace Transmit {
   cuint CELL_PARAMS  = 1;
   cuint BLOCK_PARAMS = CELL_PARAMS << 1;
   cuint AVGS         = CELL_PARAMS << 2;
   cuint FLUXES       = CELL_PARAMS << 3;
   cuint DERIV1       = CELL_PARAMS << 4;
   cuint DERIV2       = CELL_PARAMS << 5;
   cuint NBRSVEL      = CELL_PARAMS << 6;
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

   Parameters(int argc, char* argv[]);
   
   static bool add(const std::string& name,const std::string& desc,char& var,const char& defValue);
   static bool add(const std::string& name,const std::string& desc,int& var,const int& defValue);
   static bool add(const std::string& name,const std::string& desc,unsigned int& var,const unsigned int& defValue);
   static bool add(const std::string& name,const std::string& desc,long unsigned int& var,const long unsigned int& defValue);
   static bool add(const std::string& name,const std::string& desc,float& var,const float& defValue);
   static bool add(const std::string& name,const std::string& desc,double& var,const double& defValue);
   static bool add(const std::string& name,const std::string& desc,std::string& var,const std::string& defValue,const bool& multitoken=false);

   static bool finalize();
   
   static bool get(const std::string& name,char& value);
   static bool get(const std::string& name,int& value);
   static bool get(const std::string& name,unsigned int& value);
   static bool get(const std::string& name,long unsigned int& value);
   static bool get(const std::string& name,float& value);
   static bool get(const std::string& name,double& value);
   static bool get(const std::string& name,std::string& value);
   
   static bool helpMessage();
   static bool isInitialized();
   static bool parse();
   
 private:
   static int argc;                  /**< How many entries argv contains.*/
   static char** argv;              /**< Pointer to char* array containing command line parameters.*/
   
   /** Private default constructor to prevent incorrect initialization.*/
   Parameters();
   
   static bool addDefaultParameters();
};

#endif
