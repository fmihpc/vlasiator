#!/bin/bash
# Sebastian von Alfthan sebastian.von.alfthan@csc.fi
# Yann Kempf yann.kempf@fmi.fi
# Urs Ganse urs.ganse@utu.fi

# FUNCTIONS ##################################

function transferFileList {
    server=$1
    path=$2
    inputfile=$3
    localTapePath=$4

    export GLOBUS_TCP_SOURCE_RANGE=20000,20500
    export GLOBUS_TCP_PORT_RANGE=20000,20500


    while read line; do
    #inputfile produced with ls -la, get name and size. sed one-liner to remove color-codes
    file=$(echo $line| sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" | gawk '{print $9}')  
    size=$(echo $line | gawk '{print $5}')    

    #chunksize
    chunkSize=1000000000
    totalChunks=$(( 1+size/chunkSize )) 
    retval=0

    #file exists on archive folder, if it is incomplete on archive server then that is not taken into account in any way
    if [ -e ${localTapePath}/${file} ] 
    then
        tapeSize=$( ls -la  ${localTapePath}/${file} | gawk '{print $5}' )
        if [ $tapeSize -eq  $size ]
        then
          #file complete
          retval=1
          echo "$(date) ${file}: File is already transferred and on archive"
        fi
    fi
    
    
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
        echo "$(date) ${file}: File is already transferred"
        fi
    else
        #nothing has been transferred, start from beginning
        i=0
    fi



    retryIndex=0
    while [ $retval -eq 0 ]
    do
        sleep 1 #short sleep to make it easier to cancel..
        #offset into file where we start to download data
        offset=$(echo $i $chunkSize|gawk '{print $1*$2}')
        
        echo "$(date) ${file}: Starting download of chunk $((i+1))/$totalChunks at $offset " 
            localStartSize=$offset
        startTime=$( date +"%s.%N" )
        globus-url-copy  -rst -len $chunkSize  -off $offset  ${server}/${path}/${file} ./
        rc=$?
        if [[ $rc != 0 ]] ; then
        echo "Failed: globus-url-copy  -rst -len $chunkSize  -off $offset  ${server}/${path}/$file ./"
        exit $rc
        fi
        localEndSize=$( ls -la  $file | gawk '{print $5}' )
        endTime=$( date +"%s.%N" )
        echo $startTime $endTime $localStartSize $localEndSize $file $((i+1)) "$(date)" | 
        gawk '{
             dataMb=($4-$3)/(1024*1024);
             times=($2-$1); 
             print $7,$5,": chunk ",$6," downloaded at", dataMb," MB in ",times " s : ", dataMb/times, "MB/s"
            }'



        localSize=$( ls -la  $file | gawk '{print $5}' )
        if [ $localSize -eq $size ] 
        then
        #the whole file has been downloaded, excellent!
        echo "$(date) ${file}: Done"
        if [ -e ${localTapePath}/${file} ] 
        then
            echo "$(date) ${file}: WARNING file with the same name already exists on ${localTapePath} - file not moved from staging at $(pwd)"
        else
            mv ${file} ${localTapePath}/
            echo "$(date) ${file}: Moved from staging at $( pwd ) to ${localTapePath}"
        fi

        retval=1
        retryIndex=0
        else
        #file not complete
        if [ $localSize -lt $(( chunkSize + offset )) ]
        then
                    #we failed to download the whole chunk
            retryIndex=$(( retryIndex+1 ))
            echo "$(date) ${file}: Chunk transfer failed, retry number $retryIndex "
            if [ $retryIndex -gt 2 ]
            then
                echo "$(date) ${file}: Too many retries, abort. Failed on reading to offset $offset"
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
}


function transferFileListDdSsh {
    user=$1
    server=$2
    path=$3
    inputfile=$4
    localTapePath=$5


    while read line
    do
        #inputfile produced with ls -la, get name and size. sed one-liner to remove color-codes
        file=$(echo $line| sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" | gawk '{print $9}')  
        size=$(echo $line | gawk '{print $5}')    
        
        #chunksize
        chunkSize=$(( 1024 * 1024 * 1024 ))
        totalChunks=$(( 1+size/chunkSize )) 
        retval=0

        #file exists on archive folder, if it is incomplete on archive server then that is not taken into account in any way
        if [ -e ${localTapePath}/${file} ] 
        then
            tapeSize=$( ls -la  ${localTapePath}/${file} | gawk '{print $5}' )
            if [ $tapeSize -eq  $size ]
            then
            #file complete
                retval=1
                echo "$(date) ${file}: File is already transferred and on archive"
            fi
        fi
        
        
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
                echo "$(date) ${file}: File is already transferred"
            fi
        else
        #nothing has been transferred, start from beginning
        i=0
        fi
        
        retryIndex=0
        while [ $retval -eq 0 ]
        do
            sleep 1 #short sleep to make it easier to cancel..
            #Size of current transfer. chunksize for all, except last transfer
            transferSize=$(echo $i $chunkSize $size|gawk '{if(($1+1)*$2 > $3) print $3-$1*$2; else print $2;}')
            echo transferSize $transferSize
            echo "$(date) ${file}: Starting download of chunk $((i+1))/$totalChunks " 
            startTime=$( date +"%s.%N" )
       if [ $server == "localhost" ]
       then
            dd iflag=fullblock bs=${chunkSize} skip=$i count=1 if=${path}/${file} 2>> dd_read.err | tee >(md5sum > .transfer_${file}_chunk_${i}_checksum.txt) |\
               dd iflag=fullblock bs=${chunkSize} seek=$i count=1 of=${file} 2>> dd_write.err
            sourceChecksum=`cat .transfer_${file}_chunk_${i}_checksum.txt`
       else
            ssh -o Compression=no ${user}@${server} "dd iflag=fullblock bs=${chunkSize} skip=$i count=1 if=${path}/${file} | tee >(md5sum > ${path}/.transfer_${file}_chunk_${i}_checksum.txt)" 2>> dd_read.err |\
               dd iflag=fullblock bs=${chunkSize} seek=$i count=1 of=${file} 2>> dd_write.err
            sourceChecksum=`ssh ${user}@${server} "cat ${path}/.transfer_${file}_chunk_${i}_checksum.txt"`
       fi


       endTime=$( date +"%s.%N" )

            echo $startTime $endTime $chunkSize $file $((i+1)) "$(date)" |
            gawk '{
                dataMb=($3)/(1024*1024);
                times=($2-$1);
                print $6,$4,": chunk ",$5," downloaded at", dataMb," MB in ",times " s : ", dataMb/times, "MB/s"
            }'
            
            #Test if file is complete
            # Check checksums
            targetChecksum=`dd iflag=fullblock bs=${chunkSize} skip=$i count=1 if=${file} 2>>/dev/null | md5sum`
            if [[ $sourceChecksum == $targetChecksum ]]; then
               echo "Chunk $i transferred successfully"
               i=$(( i+1 ))
               retryIndex=0
            else
               retryIndex=$(( retryIndex+1 ))
               echo "$(date) ${file}: Chunk $i checksum inconsistent, retry number $retryIndex "
               if [ $retryIndex -gt 10 ]
               then
                   echo "$(date) ${file}: Too many retries, abort. Failed on reading to offset $offset"
                   retval=2
               fi
            fi

            # Initially it breaks if there is no file.
            touch $file
            localSize=$( ls -la  $file | gawk '{print $5}' )
            if [ $localSize -eq $size ] 
            then
                #the whole file has been downloaded, excellent!
                echo "$(date) ${file}: Done"
                if [ -e ${localTapePath}/${file} ] 
                then
                    echo "$(date) ${file}: WARNING file with the same name already exists on ${localTapePath} - file not moved from staging at $(pwd)"
                else
                    mv ${file} ${localTapePath}/
                    echo "$(date) ${file}: Moved from staging at $( pwd ) to ${localTapePath}"
                fi
            retval=1
            retryIndex=0
            fi
        done
    done < $inputfile
}



function transferFileListRsync {
    user=$1
    server=$2
    path=$3
    inputfile=$4
    localTapePath=$5

        while read line; do
        #inputfile produced with ls -la, get name and size. sed one-liner to remove color-codes
        file=$(echo $line| sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" | gawk '{print $9}')  
        size=$(echo $line | gawk '{print $5}')    
        retval=0
        
        #file exists on archive folder, if it is incomplete on archive server then that is not taken into account in any way
        if [ -e ${localTapePath}/${file} ] 
        then
            tapeSize=$( ls -la  ${localTapePath}/${file} | gawk '{print $5}' )
            if [ $tapeSize -eq  $size ]
            then
                #file complete
                retval=1
                echo "$(date) ${file}: File is already transferred and on archive"
            fi
        fi
        
        retryIndex=0
        while [ $retval -eq 0 ]
        do
            sleep 1 #short sleep to make it easier to cancel..
            echo "$(date) ${file}: Starting download ($retryIndex retries)" 
            #create empty file
            if [ ! -e $file ] 
            then
                touch $file
            fi
            
            startTime=$( date +"%s.%N" )
            localStartSize=$( ls -la  $file | gawk '{print $5}' ) 
            if [ $localStartSize -ne $size ] 
            then
                #start download if file is not complete
                rsync --inplace --partial   ${user}@${server}:${path}/${file} ./
                rc=$?
                if [[ $rc != 0 ]] ; then
                    echo "Failed:    rsync --inplace --partial   ${user}@${server}:${path}/${file} ./"
                fi
                localEndSize=$( ls -la  $file | gawk '{print $5}' )
                endTime=$( date +"%s.%N" )
                echo $startTime $endTime $localStartSize $localEndSize $file "$(date)" | 
                gawk '{
                     dataMb=($4-$3)/(1024*1024);
                     times=($2-$1); 
                     print $6,$5,": downloaded at", dataMb," MB in ",times " s : ", dataMb/times, "MB/s"
                    }'
            else
                echo "$(date) ${file}: File is already transferred and in staging area"
            fi
            
            localSize=$( ls -la  $file | gawk '{print $5}' )
            if [ $localSize -eq $size ] 
            then
                #the whole file has been downloaded, excellent!
                echo "$(date) ${file}: Done"
                if [ -e ${localTapePath}/${file} ] 
                then
                    echo "$(date) ${file}: WARNING file with the same name already exists on ${localTapePath} - file not moved from staging at $(pwd)"
                else
                    mv ${file} ${localTapePath}/ 2>> errlog
                    echo "$(date) ${file}: Moved from staging at $( pwd ) to ${localTapePath}"
                fi
                retval=1
                retryIndex=0
            else
                echo  "$(date) ${file}: File is not complete; $localEndSize / $size"
                retryIndex=$(( retryIndex+1 ))
                if [ $retryIndex -gt 50 ]
                then
                    echo "$(date) ${file}: Too many retries, abort. (50 max)"
                    retval=2
                fi
            fi
        done
    done < $inputfile
}


# MAIN PROGRAM ##################################
export GLOBUS_TCP_SOURCE_RANGE=20000,20500
export GLOBUS_TCP_PORT_RANGE=20000,20500

if [ ! $# -eq 5 ]
then
    cat <<EOF
transferPraceData userserver path transfer_file local_storage_path
    Script for transferring data using gridFTP or rsync (depends on machine). 
    Please run grid_proxy_init first when using the gridFTP backend.
   
    user             Username, option not used for gridftp or local-dd transfers (put arbitrary name)
    server           One of: Hermit (gridftp), Hazelhen-r (rsync), Hazelhen-ds (dd|ssh), Abel (gridftp), Sisu-g (gridftp) Sisu-r (rsync) Sisu-ds (dd|ssh) localhost-dd (local-dd) localhost-rp (rsync with special port)
    path             is a path on remote machine (e.g. /univ_1/ws1/ws/iprsalft-paper1-runs-0/2D/ecliptic/AAE)"
    transfer_file    is a file in the path on the remote machine created using ls -la *myfiles_to_transfer* > transfer_list.txt"       
    local_storage_path  is the folder where the files are ultimately copied after transfer, e.g., a tape drive. During transfer they go to the current folder. "." is also allowed.
EOF
    exit
fi
user=$1
machine=$2
path=$3
inputfile=$4
localTapePath=$5
parallelTransfers=5

#read in command line variables
if [ $machine == "Abel" ]
then
    server=gsiftp://gridftp1.prace.uio.no:2811 
    method=gridftp
elif [ $machine == "Hermit" ]
then
    server=gsiftp://gridftp-fr1.hww.de:2812
    method=gridftp
elif [ $machine == "Hazelhen-r" ]
then
    server=hazelhen.hww.de
    method=rsync
elif [ $machine == "Hazelhen-ds" ]
then
    server=hazelhen.hww.de
    method=ddssh
elif [ $machine == "Sisu-g" ]
then
    server=gsiftp://gridftp.csc.fi:2811
    method=gridftp
elif [ $machine == "Sisu-r" ]
then
    server=sisu.csc.fi
    method=rsync
elif [ $machine == "Sisu-ds" ]
then
    server=sisu.csc.fi
    method=ddssh
elif [ $machine == "localhost-dd" ]
then
    server=localhost
    method=ddssh
elif [ $machine == "localhost-rp" ]
then
    server=localhost
    method=rsync
    export RSYNC_RSH="ssh -p 1235"
else
    echo "Allowed server values are Hermit, Hazelhen-r, Hazelhen-ds, Abel, Sisu-g, Sisu-r, Sisu-ds, localhost-dd, localhost-rp"
    exit 1
fi



if [ ! -d $localTapePath ]
then
    echo "Tape path $localTapePath does not exist!"
    exit
fi


#clean up old inputfile

if [ -e $inputfile ]
then
    rm $inputfile 
fi


echo "Downloading inputfile $inputfile from $path" 
if [ $method == "gridftp" ]
then
    globus-url-copy  -rst  ${server}/${path}/$inputfile ./$inputfile
    rc=$?
    if [[ $rc != 0 ]] ; then
        echo "Failed: globus-url-copy -rst  ${server}/${path}/$inputfile ./$inputfile"
        echo "Could not download list of files"
        exit $rc
    fi
elif [ $method == "rsync" ]
then
    rsync -P --inplace  ${user}@${server}:${path}/$inputfile ./$inputfile
    rc=$?
    if [[ $rc != 0 ]] ; then
        echo "Failed:     rsync -P --inplace  ${user}@${server}:${path}/$inputfile ./$inputfile"
        echo "Could not download list of files"
        exit $rc
    fi
elif [ $method == "ddssh" ]
then
    if [ $server == "localhost" ]
    then
   cp ${path}/$inputfile ./$inputfile
   rc=$?
    else
   rsync -P --inplace  ${user}@${server}:${path}/$inputfile ./$inputfile
   rc=$?
    fi
    if [[ $rc != 0 ]] ; then
        echo "Failed: Could not download list of files"
        exit $rc
    fi
else
    echo "Failed: Unknown method $method"
    echo "Could not download list of files"
    exit $rc
fi


echo "Transferring files in $inputfile at $path" >> transferLog.txt
echo "Files staged at $( pwd) and archived at ${localTapePath}" >> transferLog.txt
cat $inputfile |gawk '{s=s+$5} END {print  NR, " files with ",s/(1024*1024*1024),"GB of data"}' >> transferLog.txt


#how many files (max) per transfer
filesPerTransfer=$( wc -l $inputfile |gawk -v p=$parallelTransfers '{print int($1/p)+1}' )
echo "$parallelTransfers parallel transfers with up to $filesPerTransfer per transfer" >> transferLog.txt

#split into transferfiles
rm -f .para_$inputfile_*
split -l $filesPerTransfer -d $inputfile .para_$inputfile_


i=0
for paraInput in .para_$inputfile_*
do
    if [ $method == "gridftp" ]
    then
        transferFileList  $server $path $paraInput $localTapePath >> transferLog.txt &
    elif [ $method == "rsync" ]
    then
        transferFileListRsync  $user $server $path $paraInput $localTapePath >> transferLog.txt &
        echo "Started background transfer-job $!"
    elif [ $method == "ddssh" ]
    then
        transferFileListDdSsh  $user $server $path $paraInput $localTapePath >> transferLog.txt &
        echo "Started background transfer-job $!"
    fi
    transferPids[$i]=$!
    i=$((i+1))
done 
