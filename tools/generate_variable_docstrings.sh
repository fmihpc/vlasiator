#!/bin/bash

# get a list of documented output variables
#
# run the vlasiator binary to get the current help text
# replace all newlines with spaces to concat multiline text for each option
# reintroduce linebreaks at the -- that start new options
# grab the variables.output option
# break at the colon
# remove duplicate spaces that come from multicolumn format
# remove remaining spaces with newlines to get one documented DRO per line
# sort alphabetically
# run through uniq to get rid of redundent empty newlines
# dump to a file
./vlasiator --help |\
   tr '\n' ' ' |\
   sed -e "s/--/\n--/g" |\
   grep variables.output |\
   cut -d : -f 2- |\
   tr -s " " |\
   tr ' ' '\n' |\
   sort |\
   uniq > /tmp/listed_output.txt

# get a list of existing output variables
#
# grep for the common pattern in datareduction/datareducer.cpp
# take the last option in quotes
# sort alphabetically
# run through uniq to get rid of redundent empty newlines
# dump to a file
grep "if(P::systemWriteAllDROs || lowercase ==" ./datareduction/datareducer.cpp |\
   rev | cut -d '"' -f 2 | rev |\
   sort |\
   uniq > /tmp/grep_output.txt

echo "variables.output string for parameters.cpp"
echo "Available ("$(date +%Y%m%d)"): "$(cat /tmp/grep_output.txt | tr '\n' ' ')


echo " "


# get a list of documented diagnostic variables
#
# run the vlasiator binary to get the current help text
# replace all newlines with spaces to concat multiline text for each option
# reintroduce linebreaks at the -- that start new options
# grab the variables.diagnostic option
# break at the colon
# remove duplicate spaces that come from multicolumn format
# remove remaining spaces with newlines to get one documented DRO per line
# sort alphabetically
# run through uniq to get rid of redundent empty newlines
# dump to a file
./vlasiator --help |\
   tr '\n' ' ' |\
   sed -e "s/--/\n--/g" |\
   grep variables.diagnostic |\
   cut -d : -f 2- |\
   tr -s " " |\
   tr ' ' '\n' |\
   sort |\
   uniq > /tmp/listed_diagnostic.txt

# get a list of existing diagnostic variables
#
# grep for the common pattern in datareduction/datareducer.cpp
# take the last option in quotes
# sort alphabetically
# run through uniq to get rid of redundent empty newlines
# dump to a file
grep "if(P::diagnosticWriteAllDROs || lowercase ==" ./datareduction/datareducer.cpp |\
   rev | cut -d '"' -f 2 | rev |\
   sort |\
   uniq > /tmp/grep_diagnostic.txt

echo "variables.diagnostic string for parameters.cpp"
echo "Available ("$(date +%Y%m%d)"): "$(cat /tmp/grep_diagnostic.txt | tr '\n' ' ')
