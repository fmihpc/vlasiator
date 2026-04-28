#!/bin/bash

if [[ z$VLASIATOR_ARCH != "z" ]]
then
   echo -e "\nLoading modules for VLASIATOR_ARCH="$VLASIATOR_ARCH"\n"
   source modules/${VLASIATOR_ARCH}.sh
else
   if [[ z$1 != "z" ]]
   then
      echo -e "\nLoading modules for passed architecture "$1"\n"
      source modules/${1}.sh
   else
      echo -e "\nDefine VLASIATOR_ARCH or pass as argument the target architecture for which to load modules!\n"
   fi
fi

echo -e "Module list now:\n"
module list

