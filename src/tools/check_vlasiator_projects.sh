#!/bin/bash

if [[ ($# == 0) || ($# -ge 3) || \
  ( ( ( ($1 != "") && ($1 != "comp") )     && \
      ( ($1 != "") && ($1 != "cfg")  ) )   || \
    ( ( ($2 != "") && ($2 != "comp") )     && \
      ( ($2 != "") && ($2 != "cfg")  ) ) ) ]]
then
   cat <<EOF

Use this script to compile all projects in a row (using aprun and 12 threads) and/or
check the validity of the default cfg file in each project folder using tools/check_vlasiator_cfg.sh.

Stores result into ascii files in a directory check_projects_compil_logs/ respectively check_projects_cfg_logs/.

Usage: give argument "comp" to compile, "cfg" to check the cfg files, both to do both.

EOF
   exit
fi

COMP=0
CFG=0

if [[ ($1 == "comp") || ($2 == "comp") ]]
then
   echo "Will compile these projects:"
   COMP=1
   rm -rf check_projects_compil_logs
   mkdir check_projects_compil_logs
fi

if [[ ($1 == "cfg") || ($2 == "cfg") ]]
then
   echo "Will check the default cfg file of these projects:"
   CFG=1
   rm -rf check_projects_cfg_logs
   mkdir check_projects_cfg_logs
fi

if [[ ($COMP == 1) || ($CFG == 1) ]]
then
   cd projects/
   LIST=$(find . -maxdepth 1 -type d | grep -v "\./\." | cut -c3- | sort)
   cd ..
   
   echo $LIST


   for i in $LIST
   do
   if [[ $COMP == 1 ]]
   then
      echo "Compiling $i"
      make clean &> /dev/null
      aprun -n 1 -d 12 make -j 12 vlasiator PROJ=${i} &> check_projects_compil_log/${i}
   fi
   
   if [[ $CFG == 1 ]]
   then
      echo "Checking default cfg of $i"
      tools/check_vlasiator_cfg.sh vlasiator_*_${i} projects/${i}/${i}.cfg &> check_projects_cfg_log/${i}
   fi
   done
   
   echo "Finished checks."
fi
