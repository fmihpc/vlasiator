/*
Background magnetic field test of GUMICS.

Copyright 1997, 1998, 1999, 2000, 2001, 2010,
2011, 2012 Finnish Meteorological Institute
*/

#include "boost/program_options.hpp"
#include "cmath"
#include "cstdlib"
#include "limits"

#include "../B0.hpp"

#define R_E 6.3712e6

using namespace std;

void get_dipole_B(
   const double x, const double y, const double z,
   double& Bx, double& By, double& Bz
) {
   const double k0 = 8e15;
   Bx = -3 * k0 * x * z / pow(x*x + y*y + z*z, 5.0/2.0);
   By = -3 * k0 * y * z / pow(x*x + y*y + z*z, 5.0/2.0);
   Bz = -3 * k0 * (z*z - (x*x + y*y + z*z) / 3.0) / pow(x*x + y*y + z*z, 5.0/2.0);
}


int main(int argc, char* argv[])
{
   boost::program_options::options_description options(
      "Usage: test1 [options (options given on the command line "
      "override options given everywhere else)], where options are:"
   );
   options.add_options()("help", "print this help message")("verbose", "print run time information");

   TB0 background_B(&options);

   /*
   Option parsing
   */
   boost::program_options::variables_map option_variables;

   // read options from command line
   boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
   boost::program_options::notify(option_variables);
   // read options from environment variables
   boost::program_options::store(boost::program_options::parse_environment(options, "GUMICS_"), option_variables);
   boost::program_options::notify(option_variables);

   // print a help message if asked
   if (option_variables.count("help") > 0) {
      return EXIT_SUCCESS;
   }

   bool verbose = false;
   if (option_variables.count("verbose") > 0) {
      verbose = true;
   }

   background_B.initialize();

   if (verbose) {
      cout << "Testing background magnetic field calculation..." << endl;
   }

   // calculate volume average for smaller and smaller cube,
   // result should converge to the value at the center
   double
      start[3] = {0, 0, 0},
      end[3] = {0, 0, 0},
      center[3] = {0, 0, 0};

   // center in positive x direction
   center[0] = 5 * R_E;
   center[1] = 0 * R_E;
   center[2] = 0 * R_E;

   if (verbose) {
      cout << "Cube center at "
         << center[0] / R_E << ", " << center[1] / R_E << ", " << center[2] / R_E << " (R_E):"
         << endl;
   }

   double
      previous_difference = std::numeric_limits<double>::max(),
      difference = std::numeric_limits<double>::max();

   for (double size = 1 * R_E; size >= 0.01 * R_E; size /= 2) {
      start[0] = center[0] - size / 2;
      start[1] = center[1] - size / 2;
      start[2] = center[2] - size / 2;
      end[0] = center[0] + size / 2;
      end[1] = center[1] + size / 2;
      end[2] = center[2] + size / 2;

      double
         result_x = 0, result_y = 0, result_z = 0,
         reference_x = 0, reference_y = 0, reference_z = 0;

      get_dipole_B(center[0], center[1], center[2], reference_x, reference_y, reference_z);

      background_B.BackgroundVolumeAverageFast(start, end, result_x, result_y, result_z);

      const double
         dx = result_x - reference_x,
         dy = result_y - reference_y,
         dz = result_z - reference_z,

      difference = sqrt(dx*dx + dy*dy + dz*dz);

      if (verbose) {
         cout << "Difference for size " << size / R_E << " (R_E): " << difference << endl;
      }

      if (difference >= previous_difference) {
         cerr << __FILE__ << ":" << __LINE__
            << " New difference (" << difference
            << ") is larger than previous (" << previous_difference
            << ")"
            << endl;
         abort();
      }

      previous_difference = difference;
   }


   // center in -y direction
   center[0] = 0 * R_E;
   center[1] = -3 * R_E;
   center[2] = 0 * R_E;

   if (verbose) {
      cout << "Cube center at "
         << center[0] / R_E << ", " << center[1] / R_E << ", " << center[2] / R_E << " (R_E):"
         << endl;
   }

   previous_difference = std::numeric_limits<double>::max(),
   difference = std::numeric_limits<double>::max();

   for (double size = 1 * R_E; size >= 0.01 * R_E; size /= 2) {
      start[0] = center[0] - size / 2;
      start[1] = center[1] - size / 2;
      start[2] = center[2] - size / 2;
      end[0] = center[0] + size / 2;
      end[1] = center[1] + size / 2;
      end[2] = center[2] + size / 2;

      double
         result_x = 0, result_y = 0, result_z = 0,
         reference_x = 0, reference_y = 0, reference_z = 0;

      get_dipole_B(center[0], center[1], center[2], reference_x, reference_y, reference_z);

      background_B.BackgroundVolumeAverageFast(start, end, result_x, result_y, result_z);

      const double
         dx = result_x - reference_x,
         dy = result_y - reference_y,
         dz = result_z - reference_z,

      difference = sqrt(dx*dx + dy*dy + dz*dz);

      if (verbose) {
         cout << "Difference for size " << size / R_E << " (R_E): " << difference << endl;
      }

      if (difference >= previous_difference) {
         cerr << __FILE__ << ":" << __LINE__
            << " New difference (" << difference
            << ") is larger than previous (" << previous_difference
            << ")"
            << endl;
         abort();
      }

      previous_difference = difference;
   }


   // +z
   center[0] = 0 * R_E;
   center[1] = 0 * R_E;
   center[2] = 2 * R_E;

   if (verbose) {
      cout << "Cube center at "
         << center[0] / R_E << ", " << center[1] / R_E << ", " << center[2] / R_E << " (R_E):"
         << endl;
   }

   previous_difference = std::numeric_limits<double>::max(),
   difference = std::numeric_limits<double>::max();

   for (double size = 1 * R_E; size >= 0.01 * R_E; size /= 2) {
      start[0] = center[0] - size / 2;
      start[1] = center[1] - size / 2;
      start[2] = center[2] - size / 2;
      end[0] = center[0] + size / 2;
      end[1] = center[1] + size / 2;
      end[2] = center[2] + size / 2;

      double
         result_x = 0, result_y = 0, result_z = 0,
         reference_x = 0, reference_y = 0, reference_z = 0;

      get_dipole_B(center[0], center[1], center[2], reference_x, reference_y, reference_z);

      background_B.BackgroundVolumeAverageFast(start, end, result_x, result_y, result_z);

      const double
         dx = result_x - reference_x,
         dy = result_y - reference_y,
         dz = result_z - reference_z,

      difference = sqrt(dx*dx + dy*dy + dz*dz);

      if (verbose) {
         cout << "Difference for size " << size / R_E << " (R_E): " << difference << endl;
      }

      if (difference >= previous_difference) {
         cerr << __FILE__ << ":" << __LINE__
            << " New difference (" << difference
            << ") is larger than previous (" << previous_difference
            << ")"
            << endl;
         abort();
      }

      previous_difference = difference;
   }


   // +x, -y, -z
   center[0] = 1 * R_E;
   center[1] = -2 * R_E;
   center[2] = 3 * R_E;

   if (verbose) {
      cout << "Cube center at "
         << center[0] / R_E << ", " << center[1] / R_E << ", " << center[2] / R_E << " (R_E):"
         << endl;
   }

   previous_difference = std::numeric_limits<double>::max(),
   difference = std::numeric_limits<double>::max();

   for (double size = 1 * R_E; size >= 0.01 * R_E; size /= 2) {
      start[0] = center[0] - size / 2;
      start[1] = center[1] - size / 2;
      start[2] = center[2] - size / 2;
      end[0] = center[0] + size / 2;
      end[1] = center[1] + size / 2;
      end[2] = center[2] + size / 2;

      double
         result_x = 0, result_y = 0, result_z = 0,
         reference_x = 0, reference_y = 0, reference_z = 0;

      get_dipole_B(center[0], center[1], center[2], reference_x, reference_y, reference_z);

      background_B.BackgroundVolumeAverageFast(start, end, result_x, result_y, result_z);

      const double
         dx = result_x - reference_x,
         dy = result_y - reference_y,
         dz = result_z - reference_z,

      difference = sqrt(dx*dx + dy*dy + dz*dz);

      if (verbose) {
         cout << "Difference for size " << size / R_E << " (R_E): " << difference << endl;
      }

      if (difference >= previous_difference) {
         cerr << __FILE__ << ":" << __LINE__
            << " New difference (" << difference
            << ") is larger than previous (" << previous_difference
            << ")"
            << endl;
         abort();
      }

      previous_difference = difference;
   }


   // cube overlaps the origin
   center[0] = 0.26 * R_E;
   center[1] = 0.26 * R_E;
   center[2] = 0.26 * R_E;

   if (verbose) {
      cout << "Cube center at "
         << center[0] / R_E << ", " << center[1] / R_E << ", " << center[2] / R_E << " (R_E):"
         << endl;
   }

   previous_difference = std::numeric_limits<double>::max(),
   difference = std::numeric_limits<double>::max();

   for (double size = 1 * R_E; size >= 0.01 * R_E; size /= 2) {
      start[0] = center[0] - size / 2;
      start[1] = center[1] - size / 2;
      start[2] = center[2] - size / 2;
      end[0] = center[0] + size / 2;
      end[1] = center[1] + size / 2;
      end[2] = center[2] + size / 2;

      double
         result_x = 0, result_y = 0, result_z = 0,
         reference_x = 0, reference_y = 0, reference_z = 0;

      get_dipole_B(center[0], center[1], center[2], reference_x, reference_y, reference_z);

      background_B.BackgroundVolumeAverageFast(start, end, result_x, result_y, result_z);

      const double
         dx = result_x - reference_x,
         dy = result_y - reference_y,
         dz = result_z - reference_z,

      difference = sqrt(dx*dx + dy*dy + dz*dz);

      if (verbose) {
         cout << "Difference for size " << size / R_E << " (R_E): " << difference << endl;
      }

      if (difference >= previous_difference) {
         cerr << __FILE__ << ":" << __LINE__
            << " New difference (" << difference
            << ") is larger than previous (" << previous_difference
            << ")"
            << endl;
         abort();
      }

      previous_difference = difference;
   }

   cout << "PASSED" << endl;

   return EXIT_SUCCESS;
}

