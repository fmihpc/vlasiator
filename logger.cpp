#include <cstdlib>
#include <iostream>
#include <fstream>

#include "logger.h"

using namespace std;

Logger::Logger() {
   initialized = false;
}

Logger::Logger(const std::string& fname) {
   setOutputFile(fname);
}

Logger::~Logger() {
   if (initialized == true) out->close();
   initialized = false;
}

bool Logger::setOutputFile(const std::string& fname) {
   out = new fstream(fname.c_str(), fstream::out);
   if (out->good() == false) {
      std::cerr << "(LOGGER) ERROR: Failed to open output file!" << std::endl;
      initialized = false;
      out = NULL;
   }
   else
     initialized = true;
   return initialized;
}

Logger& Logger::operator<<(std::ostream& os) {
   if (initialized == true) (*out) << os;
   return *this;
}

Logger& Logger::operator<<(std::ostream& (*pf)(std::ostream& )) {
   if (initialized == true) (*out) << pf;
   return *this;
}

Logger& Logger::operator<<(std::ios_base& (*pf)(std::ios_base& )) {
   if (initialized == true) (*out) << pf;
   return *this;
}

