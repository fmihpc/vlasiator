#ifndef PROFILE_H
#define PROFILE_H

namespace profile 
{
  
    bool start(const std::string &label);
    bool stop (const std::string &label);
    bool print(MPI_Comm comm);
}
#endif
