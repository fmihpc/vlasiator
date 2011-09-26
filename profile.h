/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

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
