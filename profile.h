#ifndef PROFILE_H
#define PROFILE_H

namespace profile 
{
    //start a profiling block with a certain label
    bool start(const std::string &label);
    //stop a profiling block with a certain label

    
    bool stop (const std::string &label,
               double workUnits=-1.0,
               const std::string &workUnitLabel="");
    bool print(MPI_Comm comm);
}
#endif
