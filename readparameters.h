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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <optional>
#include <stdint.h>
#include <string>
#include <typeinfo>
#include "projects/project.h"
#include <vector>

#include "common.h"
#include "version.h"

class Readparameters {
public:
   Readparameters(int cmdargc, char* cmdargv[]);
   ~Readparameters();

   /** Add a new input parameter.
    * Note that parse must be called in order for the input file(s) to be re-read.
    * Only called by the root process.
    * @param name The name of the parameter, as given in the input file(s).
    * @param desc Description for the parameter.
    * @param defValue Default value for variable.
    */
   // static void add(const std::string& name, const std::string& desc, const std::string& defValue) {
   //    int rank;
   //    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //    if (rank == MASTER_RANK) {
   //       options[name] = "";
   //       isOptionParsed[name] = false;
   //       app->add_option(
   //           name.c_str(), boost::program_options::value<std::string>(&(options[name]))->default_value(defValue),
   //           desc.c_str());
   //    }
   // }
    template <typename T> static CLI::Option* add_each_lambda(const std::string& name, const std::string& desc, T& defValue,
        std::function<void(const std::string)> lambda
                                                    ) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         // std::stringstream ss;
         //
         // static constexpr bool n = (std::is_floating_point<T>::value);
         // if (n) {
         //    ss << std::setprecision(std::numeric_limits<double>::digits10 + 1) << defValue;
         // } else {
         //    ss << defValue;
         // }
      options[name] = "";
      isOptionParsed[name] = false;
      // std::cout << name << std::endl;

      CLI::Option* opt;
      if (name.find('.') != std::string::npos) {

        auto indx = name.find('.');
        auto subcom = name.substr(0, indx);
        auto namein = name.substr(indx + 1, name.size());
        CLI::App* sub = nullptr;
        // CLI::CallbackPriority priority = CLI::CallbackPriority::First;
        //
        //
        // if (name=="proton_properties.mass") {
        //   std::cout << "priority" << std::endl;
        //   CLI::CallbackPriority priority = CLI::CallbackPriority::Last;
        // };
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
          opt=sub->add_option((dashes+namein).c_str(), defValue, desc.c_str())->each(lambda);//->callback_priority(priority)->force_callback();
          isOptionParsed[subcom]=true;
        } else {
        std::cerr << "Something went wrong with adding subcommand "+subcom+"!" << std::endl;
        abort();
        };
      } else {

        opt=app->add_option(("--"+name).c_str(), defValue, desc.c_str())->each(lambda);
        // app->callback([](){projects::Project::addParameters();});
      }  
      return opt;
      // options[name] = "";
      // isOptionParsed[name] = false;
      // app->add_option(
      //     name.c_str(), defValue,
          // desc.c_str())->each(lambda);
   }

   static CLI::App* get_app(){
      return app;
   }
   template <typename T> static CLI::Option* add(const std::string& name, const std::string& desc,
       T& value,
       std::optional<T> defval=std::nullopt, bool join=false, bool required=false

              ) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         
         // std::stringstream ss;
         //
         // static constexpr bool n = (std::is_floating_point<T>::value);
         // if (n) {
         //    ss << std::setprecision(std::numeric_limits<double>::digits10 + 1) << defValue;
         // } else {
         //    ss << defValue;
         // }
      //
      options[name] = "";
      isOptionParsed[name] = false;
      // std::cout << name << std::endl;
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
          opt = sub->add_option((dashes+namein).c_str(), value, desc.c_str())->capture_default_str(); //->each(lambda);
          isOptionParsed[subcom]=true;
        } else {
        std::cerr << "Something went wrong with adding subcommand "+subcom+"!" << std::endl;
        abort();
        };
      } else {
        opt=app->add_option(("--"+name).c_str(), value, desc.c_str())->capture_default_str();//->expected(0,-1); //->each(lambda);
      }
      if (defval && opt != nullptr){
          opt->default_val(*defval);
      } 

      return opt;

      // std::cerr << "Something went wrong with adding option " << name << std::endl;
      // abort();
      // return nullptr;
   }
   static string getPops(int i){
     return populations.at(i);
   };
   template <typename T> static void add_flag(const std::string& name, const std::string& desc,  T& defValue) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      options[name] = "";
      isOptionParsed[name] = false;
      // std::cout << name << std::endl;
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
          sub->add_flag(namein.c_str(), defValue, desc.c_str()); //->each(lambda);
          isOptionParsed[subcom]=true;
        } else {
        std::cerr << "Something went wrong with adding subcommand "+subcom+"!" << std::endl;
        abort();
        };
      } else {
        app->add_flag(name.c_str(), defValue, desc.c_str()); //->each(lambda);
      }
      // app->add_option(
      //     name.c_str(), defValue,
          // desc.c_str());
   }

   // Determine whether a given variable has been set.
   static bool isSet(const std::string& name) {
      return(app->get_option_no_throw(name)!=nullptr);
   }
   static CLI::Option* get_option(const std::string& name) {
      return app->get_option(name);
   }

   static void helpMessage();

   static bool versionMessage();

   static std::string versionInfo();
   
   static std::string configInfo();

   static void parse(bool main=false);

   static bool helpRequested;
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
