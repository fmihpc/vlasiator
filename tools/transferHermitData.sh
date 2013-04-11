#!/bin/bash
# Sebastian von Alfthan sebastian.von.alfthan@fmi.fi

server=gsiftp://gridftp-fr1.hww.de:2812
path=$1
inputfile=$2

if [ ! $# -eq 2 ]
then
cat <<EOF
transferHermitData path transfer_file

    Script for transferring data using gridFTP. Please start up the proxy using grid_proxy_init first

    path          is a path on hermit (e.g. /univ_1/ws1/ws/iprsalft-paper1-runs-0/2D/ecliptic/AAE)"
    transfer_file is a file in the path on hermit created using ls -la *myfiles* > tranfer.txt"       
EOF

exit

fi

export GLOBUS_TCP_SOURCE_RANGE=20000,20500
export GLOBUS_TCP_PORT_RANGE=20000,20500



echo "Downloading inputfile $inputfile from $path" 
globus-url-copy  -rst  gsiftp://gridftp-fr1.hww.de:2812/${path}/$inputfile ./$inputfile
cat $inputfile |gawk '{s=s+$5} END {print "Downloading ", NR, " files with ",s/(1024*1024*1024),"GB of data"}'


while read line; do
    #inputfile produced with ls -la, get name and size
    file=$(echo $line | gawk '{print $9}')  
    size=$(echo $line | gawk '{print $5}')    

    #chunksize
    chunkSize=1000000000
    totalChunks=$(( 1+size/chunkSize )) 
    retval=0
    #compute where to start download
    if [ -e $file ] 
	then
	#file exists
	localSize=$( ls -la  $file | gawk '{print $5}' )
	#Start from next possible chunkposition, some data may be lost from incomplete chunk
	i=$(( localSize / chunkSize ))
	if [ $localSize -eq  $size ]
	then
	    #file complete
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
	#offset into file where we start to download data
	offset=$(echo $i $chunkSize|gawk '{print $1*$2}')

	echo "Downloading $f chunk $((i+1))/$totalChunks at $offset for $file ... " 
        localStartSize=$offset
	startTime=$( date +"%s" )
	globus-url-copy  -rst -len $chunkSize  -off $offset  ${server}/${path}/$file ./
	localEndSize=$( ls -la  $file | gawk '{print $5}' )
	endTime=$( date +"%s" )
	echo $startTime $endTime $localStartSize $localEndSize | 
	gawk '{
             dataMb=($4-$3)/(1024*1024);
             times=($2-$1); 
             print "  ... chunk downloaded at", dataMb," MB in ",times " s : ", dataMb/times, "MB/s"
            }'



	localSize=$( ls -la  $file | gawk '{print $5}' )
	if [ $localSize -eq $size ] 
	then
            #the whole file has been downloaded, excellent!
	    echo "Done"
	    retval=1
	    retryIndex=0
	else
	    #file not complete
	    if [ $localSize -lt $(( chunkSize + offset )) ]
	    then
     	        #we failed to download the whole chunk
	        retryIndex=$(( retryIndex+1 ))
		echo "Chunk transfer failed, retry number $retryIndex "
		if [ $retryIndex -gt 2 ]
		then
		    echo "Too many retries, abort"
		    echo "Failed on reading to offset $offset"
		    retval=2
		fi
	    else
		#chunk downloaded, lets get the next one
		i=$(( i+1 ))
		retryIndex=0
	    fi 
	fi
    done

   
    
done < $inputfile



