#!/bin/bash

if [ $# -ne 1 ]; then
    echo "ERROR: USAGE: $0   your_Email-address"
    echo "  This little script may help you to take care of your work space"
    echo "  It will extend any workspace which has less than one week left and send a message to the email address"
    exit 2
fi

EMAIL_ADR=$1
MACHINE="Hornet"
CURDIR=$(pwd)
HOUR=3600
DAY=86400
WEEK=604800
MAILER=mail



for WS_NAME in $(ws_list -s)
do
    L=`ws_list | grep ^${WS_NAME}`
    DURATION=`echo $L | awk '{printf ( "%d %s %d %s", $(NF - 3), $(NF - 2),  $(NF - 1), $NF) }'`
    END_STR=`date --date="${DURATION}"  +%c`
    SECONDS_LEFT=$(($(date --date="$DURATION" +%s) - $(date +%s)))


    if (($SECONDS_LEFT < $WEEK ))
    then
	msg=$(ws_extend $WS_NAME 31 2>&1 )
	${MAILER}  ${EMAIL_ADR}  -s "Attempted to extend $WS_NAME owned by $USER "  <<EOF_2
   $USER workspace ${WS_NAME} on host ${HOST}
   $msg
EOF_2
	printf "             $WS_NAME is expiring in $DURATION. Tried to extend work space\n"
	printf "$msg \n"
    fi
    
done

cd $CURDIR