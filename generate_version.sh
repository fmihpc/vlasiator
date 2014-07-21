#!/bin/bash

. /etc/profile

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

cat >> version.cpp <<EOF
  }
  return true;
}
EOF

