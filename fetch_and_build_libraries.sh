#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$1 != "z" ]]; then
   PLATFORM=-$1
   echo "Using platform $PLATFORM as provided to the script."
else
   if [[ z$VLASIATOR_ARCH != "z" ]]; then
      PLATFORM=-$VLASIATOR_ARCH
      echo "Using platform $PLATFORM as detected from VLASIATOR_ARCH."
   else
      PLATFORM=""
      echo "No explicit $PLATFORM set. If this is not intended, pass an argument to this script or set VLASIATOR_ARCH."
   fi
fi

if [[ $PLATFORM == "-hile_cpu" || $PLATFORM == "-hile_gpu" ]]; then
   # Regular subdirectory instead of on /tmp
   rm -rf library-build
   mkdir -p library-build
   BUILDDIR=library-build
else
   BUILDDIR=`mktemp -d "${TMPDIR:-/tmp}/vlasiator-library-build-XXXXX"`
   ln -s $BUILDDIR library-build
fi

bash ./fetch_libraries.sh ${PLATFORM:1}
if [[ -f ./modules/${PLATFORM:1}.sh ]]
then
   source modules/${PLATFORM:1}.sh
fi
source build_fetched_libraries.sh ${PLATFORM:1}

# Clean up build directory
rm -rf $BUILDDIR
