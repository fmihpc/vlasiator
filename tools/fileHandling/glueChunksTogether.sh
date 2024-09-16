#!/bin/bash

if [ ! $# -eq 6 ]
then
    cat <<EOF
chunkDataToTape.sh file destination start end
    Script for writing a large file in chunks to destination using dd.
    
    source         Path to the chunks to read.
    destination    Path to which the glued chunks have to be written.
    file           File name
    start          Index of first chunk to transfer
    end            Index of last chunk to transfer
    chunkNumber    Maximum index
EOF
    exit
fi

source=$1
destination=$2
file=$3
start=$4
end=$5
chunkNumber=$6

chunkSize=$((5 * 1024*1024*1024))

echo $(date) "Transferring ${file}, with a total of $chunkNumber chunks of $chunkSize bytes. Now handling chunks $start to $end."

if [ $start -gt $end ]
then
   echo "Please swap start and end indices"
   exit
fi

if [ $start -gt $chunkNumber ]
then
   echo "Start index too large."
   exit
fi

if [ $end -gt $chunkNumber ]
then
   echo "End index too large."
   exit
fi

for chunk in $( seq $start $end )
do
   echo $(date) "Handling chunk $chunk/$chunkNumber"
   dd iflag=fullblock bs=${chunkSize} seek=$chunk count=1 if=${source}/${file}.${chunkNumber}_${chunk} of=${destination}/${file} 2>> dd_chunk.err
done

echo $(date) "Transferring ${file}, with a total of $chunkNumber chunks of $chunkSize bytes. Done handling chunks $start to $end."
