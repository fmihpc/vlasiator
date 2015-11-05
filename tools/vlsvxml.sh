#!/bin/bash

if [ $# -ne 1 ] 
then
    cat    <<EOF 

Prints out XML footer of a vlsv file

Usage: $0 file.vlsv

EOF
    exit
fi

if [ -e $1 ]
then
    tail -c $(( 100 * 1024))  $1 |strings |gawk 'BEGIN {xmlstarted=0} {if($1 == "<VLSV>") xmlstarted=1; if(xmlstarted) print $0;}' 
else
    echo "File $1 does not exist!"
fi
