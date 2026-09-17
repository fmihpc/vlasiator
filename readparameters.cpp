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
#include "common.h"
using namespace std;
bool Readparameters::helpRequested = false;
bool Readparameters::versionRequested = false;
bool Readparameters::checkCfg = false;
vector<string> Readparameters::populations = {};
map<string, string> Readparameters::subcommandDescriptions;

int Readparameters::argc;
char** Readparameters::argv;
string Readparameters::configFileName = "config.cfg";

Readparameters::Readparameters(int cmdargc, char* cmdargv[]) {
   argc = cmdargc;
   argv = cmdargv;
   addDefaultParameters();
   subcommandDescriptions["io"] = "I/O options";
   subcommandDescriptions["gridbuilder"] = "Spatial grid options";
}

Readparameters::~Readparameters() {}
Readparameters::Option* Readparameters::registerOption(const std::string& key, Option&& opt) {
   auto result = registry().insert_or_assign(key, std::move(opt));
   if (result.second) {
      registryOrder().push_back(key);
   }
   return &(result.first->second);
}

void Readparameters::addDefaultParameters() {
   add<std::string>("run_config",
                    "Configuration file, overridden by options given on the command line.",
                    configFileName, std::string("config.cfg"));
   addFlag("help", "print help message", Readparameters::helpRequested);
   addFlag("version", "print version ", Readparameters::versionRequested);
   addFlag("check_cfg", "flag whether to validate the config file", Readparameters::checkCfg);
}

void Readparameters::resetAll() {
   for (const auto& key : registryOrder()) {
      Option& opt = registry().at(key);
      opt.resetToDefault();
      opt.wasSet = false;
   }
}

std::vector<std::string> Readparameters::make_tokens(int argcIn, char** argvIn) {
   std::vector<std::string> out;
   if (argcIn <= 0) {
      return out;
   }
   QdArgParser<' '> parser(argcIn, argvIn);
   for (auto it = parser.begin(); it != parser.end(); ++it) {
      out.emplace_back(*it);
   }
   return out;
}

std::vector<std::string> Readparameters::make_tokens(const std::string& buffer) {
   std::vector<std::string> out;
   if (buffer.empty()) {
      return out;
   }
   QdArgParser<' '> parser(buffer.data(), buffer.size());
   for (auto it = parser.begin(); it != parser.end(); ++it) {
      out.emplace_back(*it);
   }
   return out;
}

void Readparameters::applyAssignment(Option& opt, const std::string& rawValue, std::set<std::string>& touched) {
   const bool first = touched.insert(opt.name).second;
   if (first) {
      opt.clearValue();
   }
   opt.wasSet = true;

   if (rawValue.size() >= 2 && rawValue.front() == '[' && rawValue.back() == ']') {
      const std::string inner = rawValue.substr(1, rawValue.size() - 2);
      std::size_t start = 0;
      while (!inner.empty() && start <= inner.size()) {
         const auto comma = inner.find(',', start);
         const std::string tok = (comma == std::string::npos) ? inner.substr(start) : inner.substr(start, comma - start);
         if (!tok.empty()) {
            opt.assignOne(tok);
         }
         if (comma == std::string::npos) {
            break;
         }
         start = comma + 1;
      }
   } else {
      opt.assignOne(rawValue);
   }
}

void Readparameters::applyArgTokens(const std::vector<std::string>& tokens, bool extras, std::vector<std::string>& invalid) {
   std::set<std::string> touched;
   std::size_t i = 0;
   while (i < tokens.size()) {
      const std::string& tok = tokens[i];
      if (tok.size() < 3 || tok[0] != '-' || tok[1] != '-') {
         ++i;
         continue;
      }
      const std::string body = tok.substr(2);
      std::string name = body;
      std::string value;
      bool hasValue = false;
      if (const auto eq = body.find('='); eq != std::string::npos) {
         name = body.substr(0, eq);
         value = body.substr(eq + 1);
         hasValue = true;
      }

      const std::string key = normalizeName(name);
      auto it = registry().find(key);
      if (it == registry().end()) {
         if (!extras) {
            invalid.push_back(name);
         }
         ++i;
         continue;
      }
      Option& opt = it->second;

      if (opt.isFlag) {
         applyAssignment(opt, hasValue ? value : std::string(""), touched);
         ++i;
         continue;
      }

      if (!hasValue) {
         if (i + 1 < tokens.size()) {
            value = tokens[i + 1];
            i += 2;
         } else {
            invalid.push_back(name + " (missing value)");
            ++i;
            continue;
         }
      } else {
         ++i;
      }
      applyAssignment(opt, value, touched);
   }
}

void Readparameters::applyConfigFile(const std::string& filename, bool extras, std::vector<std::string>& invalid) {
   std::ifstream in(filename);
   if (!in.is_open()) {
      return;
   }

   std::string section;
   std::set<std::string> touched;
   std::string rawLine;
   while (std::getline(in, rawLine)) {
      //allow # comments inline
      const auto commentPos = rawLine.find('#');
      if (commentPos != std::string::npos) {
         rawLine.erase(commentPos);
      }

      const auto firstNonSpace = rawLine.find_first_not_of(" \t\r\n");
      if (firstNonSpace == std::string::npos) {
         continue;
      }
      std::string line = rawLine;
      line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char ch) { return std::isspace(ch); }), line.end());
      if (line.empty()) {
         continue;
      }

      if (line.front() == '[') {
         if (line.back() == ']') {
            section = line.substr(1, line.size() - 2);
         } else {
            std::cerr << "Invalid configuration line, found a line starting with '[' which does not end with ']':\n"
                      << rawLine << std::endl;
         }
         continue;
      }

      const auto eq = line.find('=');
      if (eq == std::string::npos) {
         continue;
      }
      const std::string name = line.substr(0, eq);
      const std::string value = line.substr(eq + 1);
      const std::string fullName = section.empty() ? name : section + '.' + name;
      const std::string key = normalizeName(fullName);

      auto it = registry().find(key);
      if (it == registry().end()) {
         if (!extras) {
            invalid.push_back(fullName);
         }
         continue;
      }
      applyAssignment(it->second, value, touched);
   }
}

std::string Readparameters::serializeAll() {
   std::ostringstream oss;
   bool first = true;
   for (const auto& key : registryOrder()) {
      const Option& opt = registry().at(key);
      if (!opt.wasSet) {
         continue;
      }
      if (!first) {
         oss << ' ';
      }
      first = false;
      oss << "--" << key << '=' << opt.serializeValue();
   }
   return oss.str();
}

std::string Readparameters::finalizeFileName(const std::vector<std::string>& tokens) {
   static const std::string prefix = "--run_config=";
   for (std::size_t i = 0; i < tokens.size(); ++i) {
      const std::string& tok = tokens[i];
      if (tok.rfind(prefix, 0) == 0) {
         return tok.substr(prefix.size());
      }
      if (tok == "--run_config" && i + 1 < tokens.size()) {
         return tokens[i + 1];
      }
   }
   return configFileName;
}

void Readparameters::helpMessage() {
   if (!helpRequested) {
      return;
   }
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank == MASTER_RANK) {
      cout << "Usage: main [options (options given on the command line override options given "
              "everywhere else)], where options are:\n"
           << endl;

      map<string, vector<string>> bySection;
      for (const auto& key : registryOrder()) {
         const auto dot = key.find('.');
         const string section = (dot == string::npos) ? string() : key.substr(0, dot);
         bySection[section].push_back(key);
      }

      auto printGroup = [](const string& section, const vector<string>& keys) {
         if (section.empty()) {
            cout << "General options:" << endl;
         } else {
            auto it = subcommandDescriptions.find(section);
            cout << section;
            if (it != subcommandDescriptions.end() && !it->second.empty()) {
               cout << ": " << it->second;
            }
            cout << ":" << endl;
         }
         for (const auto& key : keys) {
            const Option& opt = registry().at(key);
            cout << "  --" << key;
            if (!opt.isFlag) {
               cout << "=<" << (opt.isVector ? "list" : "value") << ">";
            }
            cout << "\n      " << opt.desc;
            if (!opt.isFlag) {
               cout << " (default: " << opt.defaultStr << ")";
            }
            cout << endl;
         }
         cout << endl;
      };

      if (const auto it = bySection.find(string()); it != bySection.end()) {
         printGroup(string(), it->second);
      }
      for (const auto& entry : bySection) {
         if (entry.first.empty()) {
            continue;
         }
         printGroup(entry.first, entry.second);
      }
   }
   MPI_Finalize();
   exit(0);
}

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

std::string Readparameters::versionInfo() { return getVersion(); }

std::string Readparameters::configInfo() { return getConfig(configFileName.c_str()); }

std::vector<std::string> Readparameters::parse(bool extras) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   std::string conf;
   int confsize = 0;
   std::vector<std::string> invalid;

   if (rank == MASTER_RANK) {
      resetAll();
      const std::vector<std::string> tokens = make_tokens(argc, argv);
      configFileName = finalizeFileName(tokens);

      applyConfigFile(configFileName, extras, invalid);
      applyArgTokens(tokens, extras, invalid);
      if (!extras && !invalid.empty()) {
         std::cerr << "Error parsing config, following options are invalid:\n";
         for (const auto& name : invalid) {
            std::cerr << " " << name << "\n";
         }
         std::cerr << std::endl;
         MPI_Finalize();
         exit(1);
      }

      conf = serializeAll();
      confsize = static_cast<int>(conf.size());
   }

   MPI_Bcast(&confsize, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
   if (rank != MASTER_RANK) {
      resetAll();
      conf.resize(confsize);
   }
   MPI_Bcast(conf.data(), confsize, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

   if (rank != MASTER_RANK) {
      applyArgTokens(make_tokens(conf), true, invalid);
   }

   return invalid;
}
