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

#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include "CLI11.hpp"
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <optional>
#include <stdint.h>
#include <string>
#include <string_view>
#include <typeinfo>
#include "projects/project.h"
#include <vector>

#include "common.h"
#include "version.h"

class Readparameters {
public:
   Readparameters(int cmdargc, char* cmdargv[]);
   ~Readparameters();

   /** Add a new input parameter with lambda function that does something with each parsed value.
    * Note that parse must be called in order for the input file(s) to be re-read.
    * Only called by the root process.
    * @param name The name of the parameter, as given in the input file(s).
    * @param desc Description for the parameter.
    * @param var variable to which the parsed value is assigned.
    * @param lambda function<void(string)> for doing operations on values as they are parsed. 
    */

    template <typename T> static CLI::Option* add_each_lambda(const std::string& name, const std::string& desc, T& var,
        std::function<void(const std::string)> lambda
                                                    ) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      CLI::Option* opt;
      if (name.find('.') != std::string::npos) {

        auto indx = name.find('.');
        auto subcom = name.substr(0, indx);
        auto namein = name.substr(indx + 1, name.size());
        CLI::App* sub = nullptr;

        if (!isOptionParsed[subcom]){
          
          sub = app->add_subcommand(subcom, "");

        } else {
          sub = app->get_subcommand(subcom);
        };
        if (sub!=nullptr)
        {
          std::string dashes="-";
          if (namein.size() != 1){
            dashes+="-";
          }
          opt=sub->add_option((dashes+namein).c_str(), var, desc.c_str())->each(lambda);//->callback_priority(priority)->force_callback();
          isOptionParsed[subcom]=true;
        } else {
        std::cerr << "Something went wrong with adding subcommand "+subcom+"!" << std::endl;
        abort();
        };
      } else {
        opt=app->add_option(("--"+name).c_str(), var, desc.c_str())->each(lambda);
      }  
      return opt;
   }

   static CLI::App* get_app(){
      return app;
   }
   /** Add a new input parameter 
    * Note that parse must be called in order for the input file(s) to be re-read.
    * Only called by the root process.
    * @param name The name of the parameter, as given in the input file(s).
    * @param desc Description for the parameter.
    * @param var variable to which the parsed value is assigned.
    * @param defval Default value for variable.
    */

   template <typename T> static CLI::Option* add(const std::string& name, const std::string& desc,
       T& var, std::optional<T> defval=std::nullopt
      ) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      CLI::Option* opt;
      if (name.find('.') != std::string::npos) {

        auto indx = name.find('.');
        auto subcom = name.substr(0, indx);
        auto namein = name.substr(indx + 1, name.size());
        CLI::App* sub = nullptr;
        
        if (!isOptionParsed[subcom]){
          sub = app->add_subcommand(subcom, "");
        } else {
          sub = app->get_subcommand(subcom);
        };
        if (sub!=nullptr)
        {
          std::string dashes="-";
          if (namein.size() != 1){
            dashes+="-";
          }
          opt = sub->add_option((dashes+namein).c_str(), var, desc.c_str())->capture_default_str(); //->each(lambda);
          isOptionParsed[subcom]=true;
        } else {
        std::cerr << "Something went wrong with adding subcommand "+subcom+"!" << std::endl;
        abort();
        };
      } else {
        opt=app->add_option(("--"+name).c_str(), var, desc.c_str())->capture_default_str();//->expected(0,-1); //->each(lambda);
      }
      if (defval && opt != nullptr){
          opt->default_val(*defval);
      } 

      return opt;
   }
   static void parseComposing(){
    int rank;
    ifstream configFile;
    std::string line;
    std::string subcom;
    std::string commandName;
    bool hasConfig=true;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==MASTER_RANK && app->get_config_ptr()->count()==0) {
      hasConfig=false;
    }
    MPI_Bcast(&hasConfig,sizeof(bool),MPI_BYTE,MASTER_RANK,MPI_COMM_WORLD); 
    if (!hasConfig) {
      return;
    }
    if (rank==MASTER_RANK) {
      
      configFile.open(app->get_config_ptr()->as<std::string>());
      if (configFile.is_open()) {
        while ( getline(configFile,line) ){
          //remove white spaces 
          auto it = std::remove_if(line.begin(),line.end(),[](unsigned char x){return std::isspace(x);});
          line.erase(it,line.end());

          //skip if empty
          if (line.empty()) {
            continue;
          }

          //check if subcommand, if found make it current subcommand to be parsed
          if (line.front()=='[') {
              if (line.back()!=']') {
                std::cerr << "Invalid configuration line, found a line starting with '[' but it did not end with ']':\n"+line << std::endl;
              }
              subcom=line.substr(1,line.length() - 2);
              continue;
          } 
          
          //Find eqIndx, skip if not found
          auto eqIndx=line.find('=');
          if (eqIndx == std::string::npos ) {
            continue;
          }
          
          //Find '=' to capture option name and value as string_view for now 
          std::string_view optName=std::string_view(line).substr(0,eqIndx);
          std::string_view optValue=std::string_view(line).substr(eqIndx + 1,line.size());
          

          //Are we inside subcommand or not?
          if (!subcom.empty()) {
            commandName=subcom+'.'+std::string(optName);
          } else {
            commandName=std::string(optName);
          }

          //Skip vector values and let CLI11 handle them
          if (optValue.front()=='[') {
              //Do we need to do something?
              isOptionParsed.erase(commandName);
              continue;
          } 

          //Skip if the command is not addComposing/not known
          if (auto search=options.find(commandName); search == options.end()) {
            continue;
          }
          //this is used for checking whether to pass the handling to CLI11 later(or rather whether to rerun callback)
          if (auto search=isOptionParsed.find(commandName); search == isOptionParsed.end()) {
            continue;
          }

          if (options[commandName].empty()) {
            options[commandName]=std::string(optValue);
          } else {
            options[commandName]=options[commandName]+','+std::string(optValue);
          }
          isOptionParsed[commandName]=false;
    

        }
        configFile.close();
      }
    }
    //Add brackets and throw the values to CLI option
    for(std::map<std::string,std::string>::iterator iter = options.begin(); iter != options.end(); ++iter) {
      std::string commandName=iter->first;
      if (rank == MASTER_RANK) {
        if (auto search=isOptionParsed.find(commandName); search == isOptionParsed.end()) {
          isOptionParsed[commandName]=true;  
        } else {
          isOptionParsed[commandName]=false;
        }
      }

      MPI_Bcast(&isOptionParsed[commandName], sizeof(bool) ,MPI_BYTE, MASTER_RANK, MPI_COMM_WORLD);
      if (!isOptionParsed[commandName]) {
        iter->second='['+iter->second+']'; 
        std::string optVal=iter->second;
        int strSize=optVal.size();
        MPI_Bcast(&strSize, 1, MPI_INT,
              MASTER_RANK, MPI_COMM_WORLD);

        if (rank != MASTER_RANK) {
            optVal.resize(strSize);
        }
        MPI_Bcast(optVal.data(), strSize, MPI_CHAR,
                  MASTER_RANK, MPI_COMM_WORLD);

        auto opt= getOption(commandName);
        opt->clear();
        opt->add_result(optVal);
        opt->run_callback(); 
      }
    }


   }
    /** Add a new composing input parameter.
    * Note that parse must be called in order for the input file(s) to be re-read.
    * Only needs to be called by root process.
    * It can be defined multiple times and are all returned as a vector.
    * @param name The name of the parameter, as given in the input file(s).
    * @param desc Description for the parameter.
    * @param var variable to which the parsed value is assigned.
    * @param defval Default value for variable.
  */
   template <typename T> static CLI::Option* addComposing(const std::string& name, const std::string& desc,
       T& value, std::optional<T> defval=std::nullopt
      ) {

      //Options are added her to be handled by the config parser defined in this file NOT CLI11
      options[name] = "";
      isOptionParsed[name]=false;

      //Add an option for help print and CLI
      //Will also be used to pass the parameter to other MPI ranks 
      return add(name,desc,value,defval); 
   }
   static string getPops(int i){
     return populations.at(i);
   };
   template <typename T> static void addFlag(const std::string& name, const std::string& desc,  T& defValue) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (name.find('.') != std::string::npos) {

        auto indx = name.find('.');
        auto subcom = name.substr(0, indx);
        auto namein = name.substr(indx + 1, name.size());
        CLI::App* sub = nullptr;
        
        if (!isOptionParsed[subcom]){
          sub = app->add_subcommand(subcom, "");
        } else {
          sub = app->get_subcommand(subcom);
        };
        if (sub!=nullptr)
        {
          sub->add_flag(namein.c_str(), defValue, desc.c_str());
          isOptionParsed[subcom]=true;
        } else {
        std::cerr << "Something went wrong with adding subcommand "+subcom+"!" << std::endl;
        abort();
        };
      } else {
        app->add_flag(name.c_str(), defValue, desc.c_str());
      }

   }

   // Determine whether a given variable has been set.
   static bool isSet(const std::string& name) {
      return(app->get_option_no_throw(name)!=nullptr);
   }

   static CLI::Option* getOption(const std::string& name) {
      if (auto indx=name.find('.'); indx!=std::string::npos) {
         std::string subcomName=name.substr(0,indx);
         std::string optName=name.substr(indx+1,name.size());
         CLI::App* subcom = app->get_subcommand(subcomName);
         std::string dashes="-";
         if (optName.size() != 1){
            dashes+="-";
          }
         return subcom->get_option(dashes+optName);
      }

      std::string dashes="";
      if (name.front() != '-') {
        dashes="-";
        if (name.size() != 1){
            dashes+="-";
        }
      }
      return app->get_option(dashes+name);
   }
   static bool removeOption(const std::string& name, bool master_rank) {
      if (master_rank) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank!=MASTER_RANK){
          return false;
        }
      }
      if (auto indx=name.find('.'); indx!=std::string::npos) {
         std::string subcomName=name.substr(0,indx);
         std::string optName=name.substr(indx+1,name.size());
         CLI::App* subcom = app->get_subcommand(subcomName);
         subcom->remove_option(subcom->get_option(optName));

      } else {
        app->remove_option(app->get_option(name));
      }
      return true;
   }
   static void helpMessage();

   static bool versionMessage();

   static std::string versionInfo();
   
   static std::string configInfo();

   static void parse(bool main=false);

   static bool helpRequested;
   static bool fullHelp;
   static bool versionRequested;

   static std::vector<std::string> populations;

private:
   static int argc;    /**< How many entries argv contains.*/
   static char** argv; /**< Pointer to char* array containing command line parameters.*/
   static CLI::App* app; 

   // static boost::program_options::options_description* descriptions;
   // static boost::program_options::variables_map* variables;
   static std::map<std::string, std::string> options;
   static std::map<std::string, bool> isOptionParsed;
   static std::map<std::string, std::vector<std::string>> vectorOptions;
   static std::map<std::string, bool> isVectorOptionParsed;

   static std::string global_config_file_name;
   static std::string user_config_file_name;
   static std::string run_config_file_name;
   static void addDefaultParameters();
};

#endif
