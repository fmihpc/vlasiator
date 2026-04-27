#!/bin/bash

set -e   # Abort on error

WORKSPACE=`pwd`

if [[ z$VLASIATOR_ARCH != "z" ]]; then
   PLATFORM=-$VLASIATOR_ARCH
else
   if [[ z$1 != "z" ]]; then
      PLATFORM=-$1
   else
      PLATFORM=""
   fi
fi
echo "Using platform $PLATFORM"

if [[ $PLATFORM == "-hile_cpu" || $PLATFORM=="-hile_cpu" ]]; then
   # Regular subdirectory instead of on /tmp
   rm -rf library-build
   mkdir -p library-build
   BUILDDIR=library-build
else
   if [[ $PLATFORM == "-roihu_cpu" ]]; then
      BUILDDIR=`mktemp -d "/scratch/project_2001659/pfaukemp/vlasiator-library-build-XXXXX"`
   else
      BUILDDIR=`mktemp -d "${TMPDIR:-/tmp}/vlasiator-library-build-XXXXX"`
   fi
   ln -s $BUILDDIR library-build
fi

bash ./fetch_libraries.sh ${PLATFORM:1}
if [[ -f ./modules/${PLATFORM:1}.sh ]]
then
   bash ./modules/${PLATFORM:1}.sh
fi
bash ./build_fetched_libraries.sh ${PLATFORM:1}

# Clean up build directory
rm -rf $BUILDDIR
