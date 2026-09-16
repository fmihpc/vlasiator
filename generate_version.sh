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

# if you use client/server visit and have "cd /lustre/tmp/..." in your ~/.bashrc this workaround is needed
cd $curdir

tmpfilename=$(mktemp -u -p . version.XXXXXXXXXXXXXXXX.cpp)

cat > $tmpfilename <<EOF
#include <iostream>
#include "mpi.h"
#include <fstream>


using namespace std;

bool printVersion() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0){ 
EOF

echo "    cout << endl << \"----------- Compilation --------- \"<<endl;" >>$tmpfilename
echo "    cout <<  \"date:            $(date)\" <<endl;" >>$tmpfilename
echo "    cout <<  \"folder:          $PWD \"<<endl;" >>$tmpfilename
echo "    cout <<  \"CMP:             $1 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"CXXFLAGS:        $2 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"FLAGS:           $3 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_MPI:         $4 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_ZOLTAN:      $5 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_BOOST:       $6 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_DCCRG:       $7 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"                 commit: $8 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_FSGRID:      $9 \"<<endl;" >>$tmpfilename
echo "    cout <<  \"                 commit: ${10} \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_VLSV:        ${11} \"<<endl;" >>$tmpfilename
echo "    cout <<  \"                 commit: ${12} \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_HASHINATOR:  ${13} \"<<endl;" >>$tmpfilename
echo "    cout <<  \"                 commit: ${14} \"<<endl;" >>$tmpfilename
echo "    cout <<  \"INC_PHIPROF:     ${15} \"<<endl;" >>$tmpfilename
echo "    cout <<  \"                 commit: ${16} \"<<endl;" >>$tmpfilename

        echo "    cout << endl << \"----------- git branch --------- \"<<endl;" >>$tmpfilename
git branch  | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> $tmpfilename


echo "    cout << endl << \"----------- git log (last 10 commits) --------- \"<<endl;" >>$tmpfilename
git log --pretty=oneline | head | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> $tmpfilename


echo "    cout << endl << \"----------- module list --------- \"<<endl;" >>$tmpfilename
module list 2>&1 | gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> $tmpfilename


echo "    cout << endl << \"----------- git status --------- \"<<endl;" >>$tmpfilename
git status | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g'  |gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> $tmpfilename

echo "    cout << endl << \"----------- git diff ---------- \"<<endl;" >>$tmpfilename

echo "    const char diff_data[] = {" >> $tmpfilename
DIFF=$(git diff `git diff --name-only |grep -v generate_version.sh` | xxd -i | sed "s/0x\([0-9a-f]\{2\}\)/'\\\\x\1'/g")
if [[ -n $DIFF ]]; then
   echo -n $DIFF >> $tmpfilename
   echo "    ,0 };" >> $tmpfilename
else
   echo "    0 };" >> $tmpfilename
fi
echo "    cout << diff_data << endl;" >> $tmpfilename

cat >> $tmpfilename <<EOF
  }
  return true;
}
EOF



cat >> $tmpfilename <<EOF

std::string getVersion() {
  std::string  versionInfo;

EOF

echo "  versionInfo+=\"----------- Compilation --------- \n\";" >>$tmpfilename
echo "  versionInfo+=\"date:           $(date)\n\";" >>$tmpfilename
echo "  versionInfo+=\"CMP:            $1 \n\";" >>$tmpfilename
echo "  versionInfo+=\"CXXFLAGS:       $2 \n\";" >>$tmpfilename
echo "  versionInfo+=\"FLAGS:          $3 \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_MPI:        $4 \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_ZOLTAN:     $5 \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_BOOST:      $6 \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_DCCRG:      $7 \n\";" >>$tmpfilename
echo "  versionInfo+=\"                commit: $8 \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_FSGRID:     $9 \n\";" >>$tmpfilename
echo "  versionInfo+=\"                commit: ${10} \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_VLSV:       ${11} \n\";" >>$tmpfilename
echo "  versionInfo+=\"                commit: ${12} \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_HASHINATOR: ${13} \n\";" >>$tmpfilename
echo "  versionInfo+=\"                commit: ${14} \n\";" >>$tmpfilename
echo "  versionInfo+=\"INC_PHIPROF:    ${15} \n\";" >>$tmpfilename
echo "  versionInfo+=\"                commit: ${16} \n\";" >>$tmpfilename


echo "     versionInfo+= \"----------- git branch ---------n\";" >>$tmpfilename
git branch  | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s""\"%s\n","  versionInfo+=",$0"\\n"," ;")}' >> $tmpfilename


echo "   versionInfo+= \"----------- git log (last 10 commits) --------- \";" >>$tmpfilename
git log --pretty=oneline | head | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g' | gawk '{printf("%s\"%s\"%s\n","     versionInfo+= ",$0"\\n"," ;")}' >> $tmpfilename


echo "     versionInfo+=\"----------- module list --------- \";" >>$tmpfilename
module list 2>&1 | gawk '{printf("%s\"%s\"%s\n","   versionInfo+= ",$0"\\n"," ;")}' >> $tmpfilename


echo "     versionInfo+=\"----------- git status --------- \";" >>$tmpfilename
git status | sed 's/\"/\\"/g' | sed 's/\\\"/\\"/g'  |gawk '{printf("%s\"%s\"%s\n","   versionInfo+= ",$0"\\n"," ;")}' >> $tmpfilename


echo "   versionInfo+=\"----------- git diff ---------- \";" >>$tmpfilename

echo "    const char diff_data[] = {" >> $tmpfilename
DIFF=$(git diff `git diff --name-only |grep -v generate_version.sh` | xxd -i | sed "s/0x\([0-9a-f]\{2\}\)/'\\\\x\1'/g")
if [[ -n $DIFF ]]; then
   echo -n $DIFF >> $tmpfilename
   echo "    ,0 };" >> $tmpfilename
else
   echo "    0 };" >> $tmpfilename
fi

echo "std::string buffer;" >> $tmpfilename
echo "buffer+=diff_data;" >> $tmpfilename
echo "versionInfo+=buffer;" >> $tmpfilename

cat >> $tmpfilename <<EOF
  
  return versionInfo;
}
EOF



cat >> $tmpfilename <<EOF

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

mv $tmpfilename version.cpp
