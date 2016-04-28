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
    tail -c $(( 10 * 1024))  $1 |strings |gawk 'BEGIN {xmlstarted=0} {if( index($1, "<VLSV>")>0) xmlstarted=1; if(xmlstarted) print $0;}' 
else
    echo "File $1 does not exist!"
fi
