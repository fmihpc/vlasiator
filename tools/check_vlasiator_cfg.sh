#!/bin/bash

mpirun_cmd="aprun -n 1"
vlasiator=$1
cfg=$2


cat $cfg |gawk '{if ( $1 ~ /\[/) {prefix=substr($1,2,length($1)-2);prefix=sprintf("%s.",prefix);} else if(NF>0) printf("%s%s\n",prefix,$0)}' > .cfg_variables

$mpirun_cmd $vlasiator --help|grep "\-\-" | sed 's/--//g'  > .vlasiator_variables

cat .vlasiator_variables | gawk '{print $1}'|sort -u >.vlasiator_variable_names
cat .cfg_variables | gawk '{print $1}'|sort -u >.cfg_variable_names


echo "------------------------------------------------------------------------------------------------------------" 
echo "               Vlasiator options                              |             CFG  options" 
echo "------------------------------------------------------------------------------------------------------------" 

diff --side-by-side --suppress-common-lines .vlasiator_variable_names .cfg_variable_names
echo "------------------------------------------------------------------------------------------------------------" 