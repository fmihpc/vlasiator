#!/bin/bash

if [[ z$VLASIATOR_ARCH != "z" ]]
then
   if [[ -f modules/${VLASIATOR_ARCH}.sh ]]
   then
      echo -e "\nLoading modules for VLASIATOR_ARCH="$VLASIATOR_ARCH"\n"
      source modules/${VLASIATOR_ARCH}.sh
      echo -e "Module list now:\n"
      module list
   else
      echo -e "\nERROR, modules/"${VLASIATOR_ARCH}".sh not found!!\n"
   fi
else
   if [[ z$1 != "z" ]]
   then
      if [[ -f modules/${1}.sh ]]
      then
         echo -e "\nLoading modules for passed architecture "$1"\n"
         source modules/${1}.sh
         echo -e "Module list now:\n"
         module list
      else
         echo -e "\nERROR, modules/"${1}".sh not found!!\n"
      fi
   else
      echo -e "\nDefine VLASIATOR_ARCH or pass as argument the target architecture for which to load modules!\n"
   fi
fi


