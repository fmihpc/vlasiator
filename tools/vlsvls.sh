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
echo "Tag                     Name                          Mesh                          Vectorsize" 
echo "----------------------------------------------------------------------------------------------"
tail -c $(( 100 * 1024))  $1 |strings |gawk 'BEGIN {xmlstarted=0} {if($1 == "<VLSV>") xmlstarted=1; if(xmlstarted) print $0;}' |
sed 's|[<\"=>/]| |g' |
gawk  '{ 

printf("%-24s",$1); 

hasname=0
for (i=0;i<NF;i++)
   if ($(i)=="name" ) {
      printf("%-30s",$(i+1));  
      hasname=1;
   }
  
if ( hasname == 0)
     printf("%-30s"," ");  

hasmesh=0
for (i=0;i<NF;i++)
  if ($(i)=="mesh" ) {
     hasmesh=1;
     printf("%-30s",$(i+1)); 
  }

if ( hasmesh == 0)
    printf("%-30s"," ");  



for (i=0;i<NF;i++)
   if ($(i)=="vectorsize" ) 
     printf("%s",$(i+1));   

printf "\n";
}' 

else

echo "File $1 does not exist!"

fi
