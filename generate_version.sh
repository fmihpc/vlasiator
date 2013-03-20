#!/bin/bash


cat > version.cpp <<EOF
#include <iostream>
#include "mpi.h"

using namespace std;

bool printVersion() {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0){ 
EOF




echo "    cout <<endl;" >> version.cpp
echo "    cout << \"svn info\"<<endl;" >>version.cpp
svn info |gawk '{printf("%s\"%s\"%s\n","    cout << ",$0," << endl;")}' >> version.cpp








cat >> version.cpp <<EOF
  }
}
EOF