#!/bin/bash

# use this when running locally/on a frontend
mpirun_cmd=" "
# or pick the launcher you need
#mpirun_cmd="srun"
#mpirun_cmd="mpirun"
#mpirun_cmd="aprun"

vlasiator=$1
cfg=$2

# Return a nonzero exit code if anything invalid was found
retval=0

if [ $# -ne 2 ]
then
    cat    <<EOF
Prints out differences between parameters in a cfg file and the options that the vlasiator executable understands.

Usage: $0 vlasiator_executable cfg_file

EOF
    exit 127
fi


if [ ! -x $vlasiator ]
then
    echo "ERROR: Vlasiator executable $vlasiator does not exist or is not executable"
    exit 127
fi


if [ ! -e $cfg ]
then
    echo "ERROR: cfg file $cfg does not exist"
    exit 127
fi


# Extract the project name to filter out these options below.
project=$( cat $cfg | grep "^project" | cut --delimiter="=" -f 2 | tr -d " " )

# Check valid boundary names
if [[ $( grep "^boundary" $cfg | grep -iv Ionosphere | grep -iv Maxwellian | grep -iv Outflow | grep -iv Copysphere | wc -l ) -gt 0 ]]
then
   echo "Invalid below boundary type(s) listed below, not checking further until these are fixed"
   grep "^boundary" $cfg | grep -iv Ionosphere | grep -iv Maxwellian | grep -iv Outflow | grep -iv copysphere
   retval=1
   exit $retval
fi

# Extract the loaded system boundaries to filter out these options below.
boundaries=""
if [[ $( grep "^boundary" $cfg | grep -i Ionosphere | wc -l ) -eq 1 ]]
then
   boundaries=ionosphere
fi

if [[ $( grep "^boundary" $cfg | grep -i Maxwellian | wc -l ) -eq 1 ]]
then
   boundaries=$boundaries" maxwellian"
fi

if [[ $( grep "^boundary" $cfg | grep -i Outflow | wc -l ) -eq 1 ]]
then
   boundaries=$boundaries" outflow"
fi

if [[ $( grep "^boundary" $cfg | grep -i Copysphere | wc -l ) -eq 1 ]]
then
   boundaries=$boundaries" copysphere"
fi

# Extract the populations to filter out these options below.
populations=$( cat $cfg | grep "^ParticlePopulations" | cut --delimiter="=" -f 2 | tr -d " " )

for pop in $populations
do
   for category in ionosphere $project $boundaries properties sparse vspace
   do
      population_prefixes=$population_prefixes" "$pop"_"$category
   done
done



# List of prefixes to allow (excludes all but the active project's project options)
for prefix in $project $boundaries $population_prefixes AMR bailout boundaries fieldsolver fieldtracing gridbuilder io loadBalance Project_common restart variables vlasovsolver
do
   echo "${prefix}\."
done > .allowed_prefixes


# Long one-liner. First remove comments, then add prefix to each name and only print if line is not empty.
cat $cfg |  grep -v "^[ ]*#" |gawk '{if ( $1 ~ /\[/) {prefix=substr($1,2,length($1)-2);prefix=sprintf("%s.",prefix);} else if(NF>0) printf("%s%s\n",prefix,$0)}' > .cfg_variables

# Extract variables from Vlasiator help output.
$mpirun_cmd $vlasiator --help | grep "^  " | tr "\n" " " | sed 's/--/\n/g' | sed -e 's/ \+/ /g' | grep -v -e '^ \?$' > .vlasiator_variables


# Replace <population> with loaded populations.
# Protect special characters which otherwise (don't) trip grep in the last line of the loop.
cat .vlasiator_variables | grep "\<population\>" | while read opt
do
   option_line=""
   for pop in $populations
   do
      option=${opt/\<population\>/$pop}
      option=${option//\//\\/}
      option=${option//\[/\\[}
      option=${option//\]/\\]}
      option=${option//\*/\\*}
      option_line+=${option}"\n"
   done
   option_line=${option_line%"\n"}
   opt_protected=${opt//\//\\/}
   opt_protected=${opt_protected//\[/\\[}
   opt_protected=${opt_protected//\]/\\]}
   opt_protected=${opt_protected//\*/\\*}
   sed .vlasiator_variables -i'' -e "s/${opt_protected}/${option_line}/g"
done

# Extract option names.
cat .vlasiator_variables | gawk '{print $1}'|sort -u >.vlasiator_variable_names
# Extract default values
cat .vlasiator_variables | gawk '{if( substr($3,1,2)=="(=") { print $1,$3}  else{ print $1}}'|sort -u >.vlasiator_variable_names_default_val
# Extract option names from cfg.
cat .cfg_variables | gawk '{print $1}'|sort -u >.cfg_variable_names


# Process output and diagnostic variables
# Extract them.
cat .cfg_variables | grep "variables.output" | sort -u > .cfg_output_variable_names
cat .cfg_variables | grep "variables.diagnostic" | sort -u > .cfg_diagnostic_variable_names
# Use the intentional : character to extract the end of the help output which contains the list of variables for output and diagnostic.
cat .vlasiator_variables | grep "variables.output" | cut --delimiter=":" -f 2 | sed 's/^ //' | sed 's/ $//' | sed 's/ /\n/g' | sed 's/^/variables.output = /g' | sort -u > .vlasiator_output_variable_names
cat .vlasiator_variables | grep "variables.diagnostic" | cut --delimiter=":" -f 2 | sed 's/^ //' | sed 's/ $//' | sed 's/ /\n/g' | sed 's/^/variables.diagnostic = /g' | sort -u > .vlasiator_diagnostic_variable_names
# Extract output and diagnostic update dates
output_update=`expr match "$( cat .vlasiator_variables | grep "variables.output" )" '.*\([0-9]\{8\}\).*'`
diagnostic_update=`expr match "$( cat .vlasiator_variables | grep "variables.diagnostic" )" '.*\([0-9]\{8\}\).*'`


if [ z$PRINT_ONLY_ERRORS == "z" ]; then
   echo "------------------------------------------------------------------------------------------------------------"
   echo "Available unused options"
   echo "------------------------------------------------------------------------------------------------------------"
   comm -23 .vlasiator_variable_names .cfg_variable_names | grep -v "\." > .unused_variables
   comm -23 .vlasiator_variable_names .cfg_variable_names | grep -f .allowed_prefixes  >> .unused_variables
   grep -f .unused_variables .vlasiator_variable_names_default_val
   echo "------------------------------------------------------------------------------------------------------------"

   echo "------------------------------------------------------------------------------------------------------------"
   echo "Available unused output and diagnostic variables (as of "$output_update" resp. "$diagnostic_update")"
   echo "------------------------------------------------------------------------------------------------------------"
   comm -23 .vlasiator_output_variable_names .cfg_output_variable_names
   comm -23 .vlasiator_diagnostic_variable_names .cfg_diagnostic_variable_names
   echo "------------------------------------------------------------------------------------------------------------"
fi

output=$( comm -13 .vlasiator_variable_names .cfg_variable_names )
if [ ${#output} -ne 0 ]
then
   echo "Invalid options"
   echo "------------------------------------------------------------------------------------------------------------"
   comm -13 .vlasiator_variable_names .cfg_variable_names
   echo "------------------------------------------------------------------------------------------------------------"
   retval=1
else
   echo "------------------------------------------------------------------------------------------------------------"
   echo "No invalid options"
   echo "------------------------------------------------------------------------------------------------------------"
fi

output=$( comm -13 .vlasiator_output_variable_names .cfg_output_variable_names )
diagnostic=$( comm -13 .vlasiator_diagnostic_variable_names .cfg_diagnostic_variable_names )
if [ ${#output} -ne 0 ] || [ ${#diagnostic} -ne 0 ]
then
   echo "------------------------------------------------------------------------------------------------------------"
   echo "Invalid output or diagnostic variables (as of "$output_update" resp. "$diagnostic_update")"
   echo "------------------------------------------------------------------------------------------------------------"
   if [ ${#output} -ne 0 ]
   then
      comm -13 .vlasiator_output_variable_names .cfg_output_variable_names
   fi
   if [ ${#diagnostic} -ne 0 ]
   then
      comm -13 .vlasiator_diagnostic_variable_names .cfg_diagnostic_variable_names
   fi
   echo "------------------------------------------------------------------------------------------------------------"
   retval=2
else
   echo "------------------------------------------------------------------------------------------------------------"
   echo "No invalid output or diagnostic variables (as of "$output_update" resp. "$diagnostic_update")"
   echo "------------------------------------------------------------------------------------------------------------"
fi


rm -f .cfg_variables .cfg_variable_names .vlasiator_variables .vlasiator_variable_names .allowed_prefixes .unused_variables  .vlasiator_variable_names_default_val .cfg_output_variable_names .cfg_diagnostic_variable_names .vlasiator_diagnostic_variable_names .vlasiator_output_variable_names

exit $retval
