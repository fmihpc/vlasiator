#!/bin/bash
# Sebastian von Alfthan sebastian.von.alfthan@fmi.fi

server=gsiftp://gridftp-fr1.hww.de:2812
path=$1
inputfile=$2
echo $inputfile


echo "Downloading inputfile" 
globus-url-copy  -rst  gsiftp://gridftp-fr1.hww.de:2812/${path}/$inputfile ./$inputfile

while read line; do

#inputfile produced with ls -la
    file=$(echo $line | gawk '{print $9}')  
    size=$(echo $line | gawk '{print $5}')
  
#chunksize
    chunkSize=1000000000
    retval=0
    totalChunks=$(( 1+size/chunkSize )) 
    if [ -e $file ] 
	then
	localSize=$( ls -la  $file | gawk '{print $5}' )
	#Start from next possible chunkposition, some data may be lost from incomplete chunk
	i=$(( localSize / chunkSize ))
	if [ $localSize -eq  $size ]
	then
	    retval=1
	    echo "File $file is already transferred"
	fi
	    
    else
	#nothing has been transferred, start from beginning
	i=0
    fi

    retryIndex=0
    while [ $retval -eq 0 ]
    do
	offset=$(echo $i $chunkSize|gawk '{print $1*$2}')
	
	echo "Downloading $f chunk $((i+1))/$totalChunks at $offset for $file" 
	globus-url-copy  -rst -len $chunkSize  -off $offset  ${server}/${path}/$file ./

	localSize=$( ls -la  $file | gawk '{print $5}' )
	if [ $localSize -eq $size ] 
	then
	    echo "Done"
	    retval=1
	    retryIndex=0
	fi
	
    

	if [ $localSize -lt $(( chunkSize + offset )) ]
	then
	    if [ $retval -eq 0 ]
	    then
		echo "Last transfer failed, retry"
	        retryIndex=$(( retryIndex+1 ))
	    fi
	else
	    i=$(( i+1 ))
	    retryIndex=0
	fi

	if [ $retryIndex -eq 2 ]
	then
	    echo "Too many retries, abort"
	    echo "Failed on reading to offset $offset"
	    retVal=2
	fi
	    
    done
done < $inputfile



