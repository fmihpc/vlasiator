#!/bin/bash

if [ $# -ne 1 ] 
then
    cat    <<EOF 

Prints out parameters in a vlsv file

Usage: $0 file.vlsv

EOF
    exit
fi



if [ -e $1 ]
then
echo "Tag          Name           Value      Vectorsize" 
echo "------------------------------------------------"
tail -1000  $1 |
grep arraysize |
sed 's/[<\"=>]/ /g' |
gawk  '{ 
printf("%-14s",$1); 

for (i=0;i<NF;i++)
   if ($(i)=="name" )
      printf("%-16s",$(i+1));  

hadValue=0;
for (i=0;i<NF;i++)
  if ($(i)=="value" ) {
     hadValue=1;
     printf("%-15s",$(i+1)); 
  }

if(hadValue==0)
     printf("%-15s"," "); 

for (i=0;i<NF;i++)
   if ($(i)=="vectorsize" ) 
     printf("%s",$(i+1));   

printf "\n";
   
}' 
else

echo "File $1 does not exist!"

fi