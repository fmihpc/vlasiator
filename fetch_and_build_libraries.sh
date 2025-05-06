#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
else
   PLATFORM=""
fi
echo "Using platform $PLATFORM"

if [[ $PLATFORM == "-hile_cpu" || $PLATFORM=="-hile_cpu" ]]; then
   # Regular subdirectory instead of on /tmp
   rm -rf library-build
   mkdir -p library-build
   BUILDDIR=library-build
else
   BUILDDIR=`mktemp -d "${TMPDIR:-/tmp}/vlasiator-library-build-XXXXX"`
   ln -s $BUILDDIR library-build
fi

bash ./fetch_libraries.sh $1
bash ./build_fetched_libraries.sh $1

# Clean up build directory
rm -rf $BUILDDIR
