#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>

class Logger {
 public:
   Logger(const std::string& fname);
   ~Logger();
   
   Logger& endl(Logger& lgr);
   Logger& operator<<(std::ostream& os);
   Logger& operator<<(std::ostream& (*pf)(std::ostream& ));
   Logger& operator<<(std::ios_base& (*pf)(std::ios_base& ));

   template<typename T> Logger& operator<<(const T& val);
   
 private:
   bool initialized;
   std::fstream* out;
};

template<typename T> Logger& Logger::operator<<(const T& val) {
   if (initialized == true) (*out) << val;
   return *this;
}

#endif
