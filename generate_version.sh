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
#include <fstream>


using namespace std;

bool printVersion() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0){ 
EOF

echo "    cout << endl << \"----------- Compilation --------- \"<<endl;" >>version.cpp
echo "    cout <<  \"date:            $(date)\" <<endl;" >>version.cpp
echo "    cout <<  \"folder:          $PWD \"<<endl;" >>version.cpp
echo "    cout <<  \"CMP:             $1 \"<<endl;" >>version.cpp
echo "    cout <<  \"CXXFLAGS:        $2 \"<<endl;" >>version.cpp
echo "    cout <<  \"FLAGS:           $3 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_MPI:         $4 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_ZOLTAN:      $5 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_BOOST:       $6 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_DCCRG:       $7 \"<<endl;" >>version.cpp
echo "    cout <<  \"                 commit: $8 \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_FSGRID:      $9 \"<<endl;" >>version.cpp
echo "    cout <<  \"                 commit: ${10} \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_VLSV:        ${11} \"<<endl;" >>version.cpp
echo "    cout <<  \"                 commit: ${12} \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_HASHINATOR:  ${13} \"<<endl;" >>version.cpp
echo "    cout <<  \"                 commit: ${14} \"<<endl;" >>version.cpp
echo "    cout <<  \"INC_PHIPROF:     ${15} \"<<endl;" >>version.cpp
echo "    cout <<  \"                 commit: ${16} \"<<endl;" >>version.cpp

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
DIFF=$(git diff `git diff --name-only |grep -v generate_version.sh` | xxd -i | sed "s/0x\([0-9a-f]\{2\}\)/'\\\\x\1'/g")
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



cat >> version.cpp <<EOF

std::string getVersion() {
  std::string  versionInfo;

EOF

echo "  versionInfo+=\"----------- Compilation --------- \n\";" >>version.cpp
echo "  versionInfo+=\"date:       $(date)\n\";" >>version.cpp
echo "  versionInfo+=\"CMP:        $1 \n\";" >>version.cpp
echo "  versionInfo+=\"CXXFLAGS:   $2 \n\";" >>version.cpp
echo "  versionInfo+=\"FLAGS:      $3 \n\";" >>version.cpp
echo "  versionInfo+=\"INC_MPI:    $4 \n\";" >>version.cpp
echo "  versionInfo+=\"INC_DCCRG:  $5 \n\";" >>version.cpp
echo "  versionInfo+=\"INC_ZOLTAN: $6 \n\";" >>version.cpp
echo "  versionInfo+=\"INC_BOOST:  $7 \n\";" >>version.cpp


echo "     versionInfo+= \"----------- git branch ---------n\";" >>version.cpp
git branch  | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s""\"%s\n","  versionInfo+=",$0"\\n"," ;")}' >> version.cpp


echo "   versionInfo+= \"----------- git log (last 10 commits) --------- \";" >>version.cpp
git log --pretty=oneline | head | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s\"%s\n","     versionInfo+= ",$0"\\n"," ;")}' >> version.cpp


echo "     versionInfo+=\"----------- module list --------- \";" >>version.cpp
module list 2>&1 | gawk '{printf("%s\"%s\"%s\n","   versionInfo+= ",$0"\\n"," ;")}' >> version.cpp


echo "     versionInfo+=\"----------- git status --------- \";" >>version.cpp
git status | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g'  |gawk '{printf("%s\"%s\"%s\n","   versionInfo+= ",$0"\\n"," ;")}' >> version.cpp


echo "   versionInfo+=\"----------- git diff ---------- \";" >>version.cpp

echo "    const char diff_data[] = {" >> version.cpp
DIFF=$(git diff `git diff --name-only |grep -v generate_version.sh` | xxd -i | sed "s/0x\([0-9a-f]\{2\}\)/'\\\\x\1'/g")
if [[ -n $DIFF ]]; then
   echo -n $DIFF >> version.cpp
   echo "    ,0 };" >> version.cpp
else
   echo "    0 };" >> version.cpp
fi

echo "std::string buffer;" >> version.cpp
echo "buffer+=diff_data;" >> version.cpp
echo "versionInfo+=buffer;" >> version.cpp

cat >> version.cpp <<EOF
  
  return versionInfo;
}
EOF



cat >> version.cpp <<EOF

std::string getConfig(const char* filename) {
  std::string  configInfo;


configInfo+="\n";
configInfo+="*------------Configuration File------------*\n";
configInfo+="\n";



std::ifstream file(filename);
if (file.is_open()) {
  std::string line;
  while (std::getline(file, line)) {
      configInfo+=line.c_str();
      configInfo+="\n";

    }
  }
  file.close();


  return configInfo;
}
EOF
