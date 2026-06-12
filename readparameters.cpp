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
#include <CLI11.hpp>
#include "common.h"
#include <algorithm>

using namespace std;

// Initialize static member of class ReadParameters
bool Readparameters::helpRequested = false;
bool Readparameters::fullHelp = false;
bool Readparameters::legacyHelp = false;
bool Readparameters::checkCfg = false;
bool Readparameters::versionRequested = false;
vector<string> Readparameters::populations = {};
CLI::App app_new{"Usage: main [options (options given on the command line override "
                 "options given everywhere else)], where options are:","vlasiator"};
CLI::App* Readparameters::app = &app_new;

int Readparameters::argc;
char** Readparameters::argv;

map<string, string> Readparameters::options;
map<string, string> Readparameters::optionsComposing;
map<string, bool> Readparameters::isOptionParsed;
map<string, bool> Readparameters::isSubComParsed;
map<string, string> Readparameters::subcommandDescriptions;

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
      addDefaultParameters();
      subcommandDescriptions["io"]="I/O options";
      subcommandDescriptions["gridbuilder"]="Spatial grid options";
      app->set_config("--run_config","config.cfg","Configuration file, when passing multiple configuration files, with precedence last to first, so configs passed last will be overridden by the configs passed to it first.");
   }
   MPI_Bcast(&Readparameters::helpRequested, sizeof(bool), MPI_BYTE, 0, MPI_COMM_WORLD);
}

Readparameters::~Readparameters() {
   // delete app;
}


/** Write the descriptions of known input options to standard output if
 * an option called "help" has been read, and exit in that case.
 */
void Readparameters::helpMessage() {
   if (helpRequested) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == MASTER_RANK) {
         if (fullHelp) {
            cout << app->help("",CLI::AppFormatMode::All) << endl;
         } else if (legacyHelp) {
            for(std::map<std::string,std::string>::iterator iter = options.begin(); iter != options.end(); ++iter) {
               cout << iter->first << endl;
            }
            for(std::map<std::string,std::string>::iterator iter = optionsComposing.begin(); iter != optionsComposing.end(); ++iter) {
               cout << iter->first << endl;
            }
         } else {
            cout << app->help() << endl;
         }
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
std::string Readparameters::configInfo() { return getConfig(app->get_config_ptr()->as<std::string>().c_str()); }

/** Request Parameters to reparse input file(s). This function needs to be
 * called after new options have been added via Parameters:add functions.
 * Otherwise the values of the new options are not read. This is a collective
 * function, all processes have to see it.
 * @param extras true if unregistered options are parsed without error.
 * @return True if input file(s) were parsed successfully.
 */
void Readparameters::parse(bool extras) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == MASTER_RANK) {
      try {
         app->allow_extras(extras);
         app->allow_config_extras(extras);
         app->parse(argc, argv);
      } catch (const CLI::ConfigError& err) {
         auto config_file=app->get_config_ptr()->as<std::string>();
         std::vector<CLI::ConfigItem> givenOptions = app->get_config_formatter()->from_file(config_file);
         std::string invalidOptions="";
         bool lastSubcomValid=true;
         for (auto &opt : givenOptions) {
            auto optName=opt.fullname();
            std::string subcom="";
            if (optName.back()=='+') {
              subcom=optName.substr(0,optName.size()-3);
              if (!isSubComParsed[subcom]) {
                invalidOptions+="["+subcom+"]\n";
                lastSubcomValid=false;
              }
              continue;
            }
            else if (optName.back()=='-') { 
              lastSubcomValid=true;
              continue;
            };
            if( (options.find(optName)==options.end()) && (optionsComposing.find(optName)==optionsComposing.end()) && lastSubcomValid) {
                invalidOptions+=" "+opt.fullname()+'\n';
            } 
          }
         std::cerr << "Error parsing config, following options are invalid:\n"<<invalidOptions << std::endl;
         MPI_Finalize();
         exit(1);
      }
   }
    std::string conf;
    int confsize;

    if (rank == MASTER_RANK) {
        conf = app->config_to_str();
        confsize = conf.size();
    }

    MPI_Bcast(&confsize, 1, MPI_INT,
              MASTER_RANK, MPI_COMM_WORLD);
    if (rank != MASTER_RANK) {
        conf.resize(confsize);
    }
    MPI_Bcast(conf.data(), confsize, MPI_CHAR,
              MASTER_RANK, MPI_COMM_WORLD);

    //send the parsed configuration file as string to other ranks
    if (rank != MASTER_RANK) {
      stringstream strs(conf);
      std::istream_iterator<string> it(strs);
      std::istream_iterator<string> end;
      std::string parsed_conf = "";
      for (it = it; it != end; ++it) {
          // lists/vectors in the config are parsed to have a space between the items
          // since strings are parsed with quotes, we can prevent adding " --" to items inside the list
          // by checking if the first character is an alphabet.
          //special handling incase we have parameter with just single letter or a flag
          if (((*it).at(1)=='=') || ((*it).size()==1) ) {
            parsed_conf.append(" -" + *it);
          }
          else if (std::isalpha((*it).at(0))) {
            parsed_conf.append(" --" + *it);
          } else {
            parsed_conf.append(*it);
          }
      }

      parsed_conf.erase(remove(parsed_conf.begin(), parsed_conf.end(), '"'), parsed_conf.end());
      app->allow_extras(extras);
      app->parse(parsed_conf);
    }
}


/** Add basic program parameters **/
void Readparameters::addDefaultParameters() {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == MASTER_RANK) {
      //Remove the original CLI11-made flag, since we want to replace it with helpRequested.
      app->remove_option(app->get_help_ptr());
      Readparameters::addFlag("--help", "print this help message", Readparameters::helpRequested);
      Readparameters::addFlag("--full_help","print the full help without subcommand categorization", Readparameters::fullHelp);
      Readparameters::addFlag("--legacy_help","print all the options as a long list", Readparameters::legacyHelp);
      Readparameters::addFlag("--version", "print version information", Readparameters::versionRequested);
      Readparameters::addFlag("--check_cfg","flag whether to validate the config file",Readparameters::checkCfg); 
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
