#!/bin/bash

mpirun_cmd="aprun -n 1"
#mpirun_cmd="mpirun -np 1"
vlasiator=$1
cfg=$2


if [ $# -ne 2 ]
then
    cat    <<EOF
Prints out differences between parameters in a cfg file and the options that the vlasiator executable understands.

Usage: $0 vlasiator_executable cfg_file

EOF
    exit
fi


if [ ! -x $vlasiator ]
then
    echo "ERROR: Vlasiator executable $vlasiator does not exist or is not executable"
    exit
fi


if [ ! -e $cfg ]
then
    echo "ERROR: cfg file $cfg does not exist"
    exit
fi


# Extract the project name to filter out these options below
project=$( cat $cfg | grep "^project" | cut --delimiter="=" -f 2 | tr -d " " )

# Extract the loaded system boundaries to filter out these options below
boundaries=""
if [[ $( grep "^boundary" $cfg | grep Ionosphere | wc -l ) -eq 1 ]]
then
   boundaries=ionosphere
fi

if [[ $( grep "^boundary" $cfg | grep Maxwellian | wc -l ) -eq 1 ]]
then
   boundaries=$boundaries" maxwellian"
fi

if [[ $( grep "^boundary" $cfg | grep Outflow | wc -l ) -eq 1 ]]
then
   boundaries=$boundaries" outflow"
fi

# Extract the populations to filter out these options below
populations=$( cat $cfg | grep "^ParticlePopulations" | cut --delimiter="=" -f 2 | tr -d " " )

for pop in $populations
do
   for category in ionosphere $project $boundaries properties sparse vspace
   do
      population_prefixes=$population_prefixes" "$pop"_"$category
   done
done



# List of prefixes to allow (excludes all but the active project's project options)
for prefix in $project $boundaries $population_prefixes AMR bailout boundaries fieldsolver gridbuilder io loadBalance Project_common restart sparse variables vlasovsolver
do
   echo "${prefix}\."
done > .allowed_prefixes


#long one-liner. First remove comments, then add prefix to each name and only print if line is not empty
cat $cfg |  grep -v "^[ ]*#" |gawk '{if ( $1 ~ /\[/) {prefix=substr($1,2,length($1)-2);prefix=sprintf("%s.",prefix);} else if(NF>0) printf("%s%s\n",prefix,$0)}' > .cfg_variables

$mpirun_cmd $vlasiator --help | grep "\-\-" | sed 's/--//g'  > .vlasiator_variables


# Replace <population> with loaded populations
cat .vlasiator_variables | grep "\<population\>" | while read opt
do
   option_line=""
   for pop in $populations
   do
      option=${opt/\<population\>/$pop}
      option=${option//\//\\/}
      option=${option//\[/\\[}
      option=${option//\]/\\]}
      option_line+=${option}"\n"
   done
   option_line=${option_line%"\n"}
   opt_protected=${opt//\//\\/}
   opt_protected=${opt_protected//\[/\\[}
   opt_protected=${opt_protected//\]/\\]}
   sed .vlasiator_variables -i'' -e "s/${opt_protected}/${option_line}/g"
done


cat .vlasiator_variables | gawk '{print $1}'|sort -u >.vlasiator_variable_names
cat .vlasiator_variables | gawk '{if( substr($3,1,2)=="(=") { print $1,$3}  else{ print $1}}'|sort -u >.vlasiator_variable_names_default_val
cat .cfg_variables | gawk '{print $1}'|sort -u >.cfg_variable_names


echo "------------------------------------------------------------------------------------------------------------"
echo "Available unused options"
echo "------------------------------------------------------------------------------------------------------------"
comm -23 .vlasiator_variable_names .cfg_variable_names | grep -v "\." > .unused_variables
comm -23 .vlasiator_variable_names .cfg_variable_names | grep -f .allowed_prefixes  >> .unused_variables
grep -f .unused_variables .vlasiator_variable_names_default_val
echo "------------------------------------------------------------------------------------------------------------"

output=$( comm -13 .vlasiator_variable_names .cfg_variable_names )
if [ ${#output} -ne 0 ]
then
   echo "Invalid options"
   echo "------------------------------------------------------------------------------------------------------------"
   comm -13 .vlasiator_variable_names .cfg_variable_names
   echo "------------------------------------------------------------------------------------------------------------"
else
   echo "No invalid options"
fi


rm .cfg_variables .cfg_variable_names .vlasiator_variables .vlasiator_variable_names .allowed_prefixes .unused_variables  .vlasiator_variable_names_default_val
