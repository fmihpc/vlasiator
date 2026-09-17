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
#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mpi.h>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
#include "common.h"
#include "version.h"
#include "qdparser.h"

class Readparameters {
public:
   Readparameters(int cmdargc, char* cmdargv[]);
   ~Readparameters();

   template <typename T> struct is_vector : public std::false_type {};

   template <typename T, typename A> struct is_vector<std::vector<T, A>> : public std::true_type {};
   struct Option {
      Option* expected(long, long) { return this; }
      Option* required(bool = true) { return this; }
      std::string name;
      std::string desc;
      std::string defaultStr;
      bool isVector = false;
      bool isFlag = false;
      bool wasSet = false;
      std::function<void()> resetToDefault;
      std::function<void()> clearValue;
      std::function<void(const std::string&)> assignOne;
      std::function<std::string()> serializeValue;
   };

   template <typename T>
   static Option* add(const std::string& name, const std::string& desc, T& var,
                       std::optional<T> defval = std::nullopt) {
      const std::string key = normalizeName(name);
      T baseline = defval.has_value() ? *defval : var;
      var = baseline;

      Option opt;
      opt.name = key;
      opt.desc = desc;
      opt.isVector = is_vector<T>::value;
      opt.defaultStr = formatValue(baseline);
      opt.resetToDefault = [&var, baseline]() { var = baseline; };
      if constexpr (is_vector<T>::value) {
         using ElemT = typename T::value_type;
         opt.clearValue = [&var]() { var.clear(); };
         opt.assignOne = [&var](const std::string& tok) { var.push_back(parseScalar<ElemT>(tok)); };
      } else {
         opt.clearValue = [&var]() { var = T{}; };
         opt.assignOne = [&var](const std::string& tok) { var = parseScalar<T>(tok); };
      }
      opt.serializeValue = [&var]() { return formatValue(var); };
      return registerOption(key, std::move(opt));
   }

   template <typename T>
   static Option* addComposing(const std::string& name, const std::string& desc, T& var,
                                std::optional<T> defval = std::nullopt) {
      return add(name, desc, var, defval);
   }

   static Option* addFlag(const std::string& name, const std::string& desc, bool& var) {
      const std::string key = normalizeName(name);
      var = false;

      Option opt;
      opt.name = key;
      opt.desc = desc;
      opt.isFlag = true;
      opt.defaultStr = "false";
      opt.resetToDefault = [&var]() { var = false; };
      opt.clearValue = [&var]() { var = false; };
      opt.assignOne = [&var](const std::string& tok) { var = tok.empty() ? true : parseScalar<bool>(tok); };
      opt.serializeValue = [&var]() { return std::string(var ? "true" : "false"); };
      return registerOption(key, std::move(opt));
   }

   static bool isSet(const std::string& name) {
      auto it = registry().find(normalizeName(name));
      return it != registry().end() && it->second.wasSet;
   }

   static void helpMessage();

   static bool versionMessage();

   static std::string versionInfo();

   static std::string configInfo();

   static std::vector<std::string> parse(bool extras = false);

   static void parseComposing() {}
   static bool helpRequested;
   static bool versionRequested;
   static bool checkCfg;
   static std::vector<std::string> populations;
   static std::map<std::string, std::string> subcommandDescriptions;

private:
   static int argc; 
   static char** argv;
   static std::string configFileName;

   static std::string normalizeName(const std::string& name) {
      std::size_t i = 0;
      while (i < name.size() && name[i] == '-') {
         ++i;
      }
      return name.substr(i);
   }

   static std::map<std::string, Option>& registry() {
      static std::map<std::string, Option> reg;
      return reg;
   }
   static std::vector<std::string>& registryOrder() {
      static std::vector<std::string> order;
      return order;
   }
   static Option* registerOption(const std::string& key, Option&& opt);

   static void addDefaultParameters();
   static void resetAll();
   static std::string serializeAll();
   static std::string finalizeFileName(const std::vector<std::string>& tokens);
   static std::vector<std::string> make_tokens(int argcIn, char** argvIn);
   static std::vector<std::string> make_tokens(const std::string& buffer);
   static void applyAssignment(Option& opt, const std::string& rawValue, std::set<std::string>& touched);
   static void applyArgTokens(const std::vector<std::string>& tokens, bool extras, std::vector<std::string>& invalid);
   static void applyConfigFile(const std::string& filename, bool extras, std::vector<std::string>& invalid);

   template <typename T> static T parseScalar(const std::string& tok) {
      if constexpr (std::is_same_v<T, std::string>) {
         return tok;
      } else if constexpr (std::is_same_v<T, bool>) {
         std::string low = tok;
         std::transform(low.begin(), low.end(), low.begin(), [](unsigned char c) { return std::tolower(c); });
         if (low == "1" || low == "true" || low == "yes" || low == "on") {
            return true;
         }
         if (low == "0" || low == "false" || low == "no" || low == "off") {
            return false;
         }
         throw std::runtime_error("cannot parse '" + tok + "' as a boolean");
      } else if constexpr (std::is_floating_point_v<T>) {
         return static_cast<T>(std::stold(tok));
      } else if constexpr (std::is_integral_v<T>) {
         if constexpr (std::is_unsigned_v<T>) {
            return static_cast<T>(std::stoull(tok));
         } else {
            return static_cast<T>(std::stoll(tok));
         }
      } else if constexpr (std::is_enum_v<T>) {
         return static_cast<T>(std::stoll(tok));
      } else {
         static_assert(!sizeof(T), "Readparameters::what the hell is this type!!!");
      }
   }

   template <typename T> static std::string formatScalar(const T& val) {
      if constexpr (std::is_same_v<T, std::string>) {
         return val;
      } else if constexpr (std::is_same_v<T, bool>) {
         return val ? "true" : "false";
      } else if constexpr (std::is_floating_point_v<T>) {
         std::ostringstream ss;
         ss << std::setprecision(std::numeric_limits<T>::digits10 + 1) << val;
         return ss.str();
      } else if constexpr (std::is_enum_v<T>) {
         return std::to_string(static_cast<std::underlying_type_t<T>>(val));
      } else {
         return std::to_string(val);
      }
   }

   template <typename T> static std::string formatValue(const T& val) {
      if constexpr (is_vector<T>::value) {
         std::string out = "[";
         for (std::size_t i = 0; i < val.size(); ++i) {
            if (i) {
               out += ',';
            }
            out += formatScalar(val[i]);
         }
         out += ']';
         return out;
      } else {
         return formatScalar(val);
      }
   }
};

#endif
