#!/bin/bash
#
# This file is part of Vlasiator.
# Copyright 2010-2016 Finnish Meteorological Institute
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

curdir=$(pwd)

. /etc/profile

# if you use client/server visit and have "cd /lustre/tmp/..." in your ~/.bashrc this workaround is needed
cd $curdir

cat > version.cpp <<EOF
#include <iostream>
#include "mpi.h"

using namespace std;

bool printVersion() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0){ 
EOF

echo "    cout << endl << \"----------- Compilation --------- \"<<endl;" >>version.cpp
echo "    cout <<  \"date:       $(date)\" <<endl;" >>version.cpp
echo "    cout <<  \"folder:     $PWD \"<<endl;" >>version.cpp
echo "    cout <<  \"CMP:        $1 \"<<endl;" >>version.cpp
echo "    cout <<  \"CXXFLAGS:   $2 \"<<endl;" >>version.cpp
echo "    cout <<  \"FLAGS:      $3 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_MPI:    $4 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_DCCRG:  $5 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_ZOLTAN: $6 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_BOOST:  $7 \"<<endl;" >>version.cpp


echo "    cout << endl << \"----------- git branch --------- \"<<endl;" >>version.cpp
git branch  | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> version.cpp


echo "    cout << endl << \"----------- git log (last 10 commits) --------- \"<<endl;" >>version.cpp
git log --pretty=oneline | head | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> version.cpp


echo "    cout << endl << \"----------- module list --------- \"<<endl;" >>version.cpp
module list 2>&1 | gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> version.cpp


echo "    cout << endl << \"----------- git status --------- \"<<endl;" >>version.cpp
git status | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g'  |gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> version.cpp

echo "    cout << endl << \"----------- git diff ---------- \"<<endl;" >>version.cpp

echo "    const char diff_data[] = {" >> version.cpp
DIFF=$(git diff `git diff --name-only |grep -v generate_version.sh` | xxd -i)
if [[ -n $DIFF ]]; then
   echo -n $DIFF >> version.cpp
   echo "    ,0 };" >> version.cpp
else
   echo "    0 };" >> version.cpp
fi
echo "    cout << diff_data << endl;" >> version.cpp

cat >> version.cpp <<EOF
  }
  return true;
}
EOF

