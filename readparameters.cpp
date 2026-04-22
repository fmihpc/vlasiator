/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "readparameters.h"
#include "CLI11.hpp"
#include "common.h"
#include "object_wrapper.h"
#include "projects/project.h"
#include <algorithm>

using namespace std;
// namespace PO = boost::program_options;

// Initialize static member of class ReadParameters
bool Readparameters::helpRequested = false;
bool Readparameters::versionRequested = false;
vector<string> Readparameters::populations = {};
CLI::App app_new{"Usage: main [options (options given on the command line override "
                 "options given everywhere else)], where options are:"};
CLI::App* Readparameters::app = &app_new;
// PO::options_description* Readparameters::descriptions = NULL;
// PO::variables_map* Readparameters::variables = NULL;

string Readparameters::global_config_file_name = "";
string Readparameters::user_config_file_name = "";
string Readparameters::run_config_file_name = "";

int Readparameters::argc;
char** Readparameters::argv;

map<string, string> Readparameters::options;
map<string, bool> Readparameters::isOptionParsed;
map<string, vector<string>> Readparameters::vectorOptions;
map<string, bool> Readparameters::isVectorOptionParsed;

/** Constructor for class ReadParameters.
 * The constructor defines some default parameters and parses the input files.
 * @param cmdargc Command line argc.
 * @param cmdargv Command line argv.
 */
Readparameters::Readparameters(int cmdargc, char* cmdargv[]) {
   argc = cmdargc;
   argv = cmdargv;
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == MASTER_RANK) {
      // variables = new PO::variables_map;
      addDefaultParameters();
      app->set_config("--run_config");

      // Read options from command line, first time for help message parsing, second time in parse() below.
      // PO::store(PO::command_line_parser(argc, argv).options(*descriptions).allow_unregistered().run(), *variables);
      // PO::notify(*variables);
   }
   MPI_Bcast(&Readparameters::helpRequested, sizeof(bool), MPI_BYTE, 0, MPI_COMM_WORLD);
}

Readparameters::~Readparameters() {
   // delete app;
}

/** Add a new composing input parameter.
 * Note that parse must be called in order for the input file(s) to be re-read.
 * Only needs to be called by root process.
 * It can be defined multiple times and are all returned as a vector.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 */
// void Readparameters::addComposing(const string& name, const string& desc) {
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if (rank == MASTER_RANK) {
//       isVectorOptionParsed[name] = false;
//       descriptions->add_options()(name.c_str(), PO::value<vector<string>>(&(vectorOptions[name]))->composing(),
//                                   desc.c_str());
//    }
// }

/** Write the descriptions of known input options to standard output if
 * an option called "help" has been read, and exit in that case.
 */
void Readparameters::helpMessage() {
   if (helpRequested) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == MASTER_RANK) {
         cout << app->help() << endl;
      }
      MPI_Finalize();
      exit(0);
   }
}

/** Write version information to standard output if
 * an option called "version" has been read.
 * @return If true, option called "version" was found and descriptions were
 * written to standard output.
 */
bool Readparameters::versionMessage() {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == MASTER_RANK) {
      if (Readparameters::versionRequested) {
         printVersion();
         return true;
      }
      return false;
   }
   return true;
}

/** Helper wrapper function to get version info
   @return std string with the version information
 */
std::string Readparameters::versionInfo() { return getVersion(); }

/** Helper wrapper function to get the config info
   @return std string with the config information
 */
std::string Readparameters::configInfo() { return getConfig(run_config_file_name.c_str()); }

/** Request Parameters to reparse input file(s). This function needs to be
 * called after new options have been added via Parameters:add functions.
 * Otherwise the values of the new options are not read. This is a collective
 * function, all processes have to see it.
 * @param needsRunConfig Whether or not this program can run without a runconfig
 * file (esm: vlasiator can't, but particle pusher can).
 * @param allowUnknown true if unregistered options are parsed without error.
 * @return True if input file(s) were parsed successfully.
 */
void Readparameters::parse(bool main) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == MASTER_RANK) {
      try {

         // app->allow_config_extras();
         app->parse(argc, argv);


         string conf1 = app->config_to_str();
         int confsize1 = conf1.size();

      } catch (const CLI::ParseError& e) {
         app->exit(e);
      }
   }
  
   // if (main) {

std::string conf;
int confsize;

if (rank == MASTER_RANK) {
    conf = app->config_to_str();
    confsize = conf.size();

}

/* One thread per rank performs MPI calls */
// #pragma omp master
// {
    /* Broadcast size */
    MPI_Bcast(&confsize, 1, MPI_INT,
              MASTER_RANK, MPI_COMM_WORLD);

    /* Non-root ranks allocate space */
    if (rank != MASTER_RANK) {
        conf.resize(confsize);
    }

    /* Broadcast data */
    MPI_Bcast(conf.data(), confsize, MPI_CHAR,
              MASTER_RANK, MPI_COMM_WORLD);
// }
// #pragma omp barrier
      if (rank != MASTER_RANK) {
         stringstream strs(conf);
         std::istream_iterator<string> it(strs);
         std::istream_iterator<string> end;
         std::string parsed_conf = "";
         for (it = it; it != end; ++it) {
            // lists/vectors in the config are parsed to have a space between the items
            // since strings are parsed with quotes, we can prevent adding " --" to items inside the list
            // by checking if the first character is an alphabet.
            if (std::isalpha((*it).at(0))) {
               parsed_conf.append(" --" + *it);
            } else {
               parsed_conf.append(*it);
            }
         }

         // std::replace(parsed_conf.begin(), parsed_conf.end(), '"','\0' );
         parsed_conf.erase(remove(parsed_conf.begin(), parsed_conf.end(), '"'), parsed_conf.end());
         std::cout << parsed_conf << std::endl;
         // app->allow_extras();
         // string part = " --ParticlePopulations=[proton]";
         // app->parse(parsed_conf + part);
         // app->remove_option(app->get_option("--ParticlePopulations"));
         app->allow_extras(!main);
         // getObjectWrapper().project->addParameters();
         //

         // projects::createProject();




         // std::cout << "CSTR=" << parsed_conf << std::endl;
         app->parse(parsed_conf);

         // // getObjectWrapper().project->getParameters();
         // auto subs = app->get_subcommands();
         // auto opts1 = app->get_options();
         // for (auto opt : opts1) {
         //    std::cout << "OPTION=" << opt->get_name() << std::endl;
         // }
         // std::cout << "SUBSSIZE=" << subs.size() << std::endl;
         // string last_name;
         // for (auto sub : subs) {
         //    if (last_name == sub->get_name()) {
         //       continue;
         //    }
         //    auto opts = sub->get_options();
         //    std::cout << "OPTSIZE=" << opts.size() << std::endl;
         //    std::cout << sub->get_name() << std::endl;
         //    for (auto opt : opts) {
         //       std::cout << "SUBCOM=" << sub->get_name() << "      OPTION=" << opt->get_name() << std::endl;
         //    }
         //    last_name = sub->get_name();
         // }
         // abort();
      }
   // }
}
// bool Readparameters::parse(const bool needsRunConfig, const bool allowUnknown) {
//
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    if (rank == MASTER_RANK) {
//       // Read options from command line:
//       PO::store(PO::parse_command_line(argc, argv, *descriptions), *variables);
//       PO::notify(*variables);
//       // Read options from environment variables:
//       PO::store(PO::parse_environment(*descriptions, "MAIN_"), *variables);
//       PO::notify(*variables);
//       // Read options from run config file:
//       if (run_config_file_name.size() > 0) {
//          ifstream run_config_file(run_config_file_name.c_str(), fstream::in);
//          if (run_config_file.good()) {
//             PO::store(PO::parse_config_file(run_config_file, *descriptions, allowUnknown), *variables);
//             PO::notify(*variables);
//             run_config_file.close();
//          } else {
//             cerr << __FILE__ << ":" << __LINE__ << "Couldn't open or read run config file " + run_config_file_name
//                  << endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//          }
//       }
//       // Read options from user config file:
//       if (user_config_file_name.size() > 0) {
//          ifstream user_config_file(user_config_file_name.c_str(), fstream::in);
//          if (user_config_file.good()) {
//             PO::store(PO::parse_config_file(user_config_file, *descriptions, allowUnknown), *variables);
//             PO::notify(*variables);
//             user_config_file.close();
//          } else {
//             cerr << __FILE__ << ":" << __LINE__ << "Couldn't open or read user config file " + user_config_file_name
//                  << endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//          }
//       }
//       // Read options from global config file:
//       if (global_config_file_name.size() > 0) {
//          ifstream global_config_file(global_config_file_name.c_str(), fstream::in);
//          if (global_config_file.good()) {
//             PO::store(PO::parse_config_file(global_config_file, *descriptions, allowUnknown), *variables);
//             PO::notify(*variables);
//             global_config_file.close();
//          } else {
//             cerr << __FILE__ << ":" << __LINE__ << "Couldn't open or read global config file " +
//             global_config_file_name
//                  << endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//          }
//       }
//    }
//
//    // Check if the user has specified --version
//    bool hasVersionOption = versionMessage();
//    MPI_Bcast(&hasVersionOption, sizeof(bool), MPI_BYTE, 0, MPI_COMM_WORLD);
//    if (hasVersionOption) {
//       MPI_Finalize();
//       exit(0);
//    }
//
//    // Require that there is a run config file.
//    // There are so many options so it is unlikely that one would like to define
//    // all on the command line, or only in global/user run files.
//    bool hasRunConfigFile = (run_config_file_name.size() > 0);
//    MPI_Bcast(&hasRunConfigFile, sizeof(bool), MPI_BYTE, 0, MPI_COMM_WORLD);
//    if (needsRunConfig && !hasRunConfigFile && !helpRequested) {
//       if (rank == MASTER_RANK) {
//          cout << "Run config file required. Use --help to list all options" << endl;
//       }
//       MPI_Finalize();
//       exit(0);
//    }
//
//    int nOptionsToBroadcast;
//    int vectorSize;
//    const int maxStringLength = 1024;
//    char value[maxStringLength];
//    char name[maxStringLength];
//    value[maxStringLength - 1] = '\0';
//    name[maxStringLength - 1] = '\0';
//
//    /*
//      loop through options and broadcast all options from root rank to the others
//    Separate bcasts not optimal from performance point of view, but parse is
//    normally just called a few times so it should not matter.
//    */
//
//    // count number of options not parsed/broadcasted previously
//    if (rank == MASTER_RANK) {
//       nOptionsToBroadcast = 0;
//       for (map<string, bool>::iterator ip = isOptionParsed.begin(); ip != isOptionParsed.end(); ++ip) {
//          if (!ip->second)
//             nOptionsToBroadcast++;
//       }
//    }
//
//    MPI_Bcast(&nOptionsToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    if (rank == MASTER_RANK) {
//       // iterate through map and bcast cstrings of key/value pairs not parsed before
//       for (map<string, string>::iterator p = options.begin(); p != options.end(); ++p) {
//
//          if (!isOptionParsed[p->first]) {
//             strncpy(name, p->first.c_str(), maxStringLength - 1);
//             strncpy(value, p->second.c_str(), maxStringLength - 1);
//             MPI_Bcast(name, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//             MPI_Bcast(value, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//             isOptionParsed[p->first] = true;
//          }
//       }
//    } else {
//       // receive new options
//       for (int p = 0; p < nOptionsToBroadcast; p++) {
//          MPI_Bcast(name, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//          MPI_Bcast(value, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//          string sName(name);
//          string sValue(value);
//          options[sName] = sValue;
//       }
//    }
//
//    /*
//      loop through vector options and broadcast all vector options from root rank
//    to the others. Separate bcasts not optimal from performance point of view,
//    but parse is normally just called a few times so it should not matter.
//    */
//
//    // count number of vector options not parsed/broarcasted previously
//    if (rank == MASTER_RANK) {
//       nOptionsToBroadcast = 0;
//       for (map<string, bool>::iterator ip = isVectorOptionParsed.begin(); ip != isVectorOptionParsed.end(); ++ip) {
//          if (!ip->second)
//             nOptionsToBroadcast++;
//       }
//    }
//
//    // root broadcasts its new vector values
//    MPI_Bcast(&nOptionsToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    if (rank == MASTER_RANK) {
//       // iterate through map and bcast cstrings of key/value pairs not parsed before
//       for (map<string, vector<string>>::iterator p = vectorOptions.begin(); p != vectorOptions.end(); ++p) {
//          if (!isVectorOptionParsed[p->first]) {
//             strncpy(name, p->first.c_str(), maxStringLength - 1);
//             MPI_Bcast(name, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//             vectorSize = vectorOptions[p->first].size();
//             MPI_Bcast(&vectorSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//             for (vector<string>::iterator v = vectorOptions[p->first].begin(); v != vectorOptions[p->first].end();
//                  ++v) {
//                strncpy(value, v->c_str(), maxStringLength - 1);
//                MPI_Bcast(value, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//             }
//             isVectorOptionParsed[p->first] = true;
//          }
//       }
//    } else {
//       // others receive new options
//       for (int p = 0; p < nOptionsToBroadcast; p++) {
//          MPI_Bcast(name, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//          string sName(name);
//          MPI_Bcast(&vectorSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//          for (int i = 0; i < vectorSize; i++) {
//             MPI_Bcast(value, maxStringLength, MPI_CHAR, 0, MPI_COMM_WORLD);
//             string sValue(value);
//             vectorOptions[sName].push_back(sValue);
//          }
//       }
//    }
//
//    return true;
// }

// add names of input files
void Readparameters::addDefaultParameters() {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == MASTER_RANK) {
      app->remove_option(app->get_help_ptr());
      Readparameters::add_flag("--help", "print this help message", Readparameters::helpRequested);
      // std::cout << "INSIDE DEFAULT PARAM ADD" << std::endl;
      // Readparameters::app->get_option("--help")->each([](const string){
      //   std::cout << "test" << std::endl;
      //   Readparameters::helpRequested=true;
      // });
      Readparameters::add_flag("--version", "print version information", Readparameters::versionRequested);

      // // Parameters which set the names of the configuration file(s):
      // descriptions->add_options()(
      //     "global_config", PO::value<string>(&global_config_file_name)->default_value(""),
      //     "read options from the global configuration file arg (relative to the current working directory). Options "
      //     "given in this file are overridden by options given in the user's and run's configuration files and by "
      //     "options given in environment variables (prefixed with MAIN_) and the command line")(
      //     "user_config", PO::value<string>(&user_config_file_name)->default_value(""),
      //     "read options from the user's configuration file arg (relative to the current working directory). Options "
      //     "given in this file override options given in the global configuration file. Options given in this file "
      //     "are "
      //     "overridden by options given in the run's configuration file and by options given in environment "
      //     "variables "
      //     "(prefixed with MAIN_) and the command line")("run_config",
      //                                                   PO::value<string>(&run_config_file_name)->default_value(""),
      //                                                   "read options from the run's configuration file arg "
      //                                                   "(relative to the current working directory). Options "
      //                                                   "given in this file override options given in the user's "
      //                                                   "and global configuration files. Options given in "
      //                                                   "this override options given in the user's and global "
      //                                                   "configuration files. Options given in this file are "
      //                                                   "overridden by options given in environment variables "
      //                                                   "(prefixed with MAIN_) and the command line");
   }
}
