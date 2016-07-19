#!/bin/bash

if [ $# -ne 1 ]; then
    echo "ERROR: USAGE: $0   your_Email-address"
    echo "  This little script may help you to take care of your work space"
    echo "  It will send a reminder to the email address, please add to crontab with crontab -e"
    echo "  */5 * * * * vlasiator/tools/ws_reminder.sh vlasiator-runs@fmihpc.flowdock.com"
    exit 2
fi

EMAIL_ADR=$1
MACHINE="Hornet"
CURDIR=$(pwd)
HOUR=3600
DAY=86400
WEEK=604800
MAILER=mail


cd $HOME
if [ ! -d .ws_remind ]
then
    mkdir .ws_remind
    
fi
cd .ws_remind

if [ -e last_check ]
then
    SECONDS_SINCE_LAST_CHECK=$(($(date +%s) - $(cat last_check)))
else
    SECONDS_SINCE_LAST_CHECK=$(date +%s)
fi

if (( $SECONDS_SINCE_LAST_CHECK < 3 ))
then
    #avoid repeated tests when, e.g., bashrc is called multiple times
    exit
fi
#store when this check was done
date +%s > last_check


printf "Workspace status:\n"
for WS_NAME in $(ws_list -s)
do
    L=`ws_list | grep ^${WS_NAME}`
    DURATION=`echo $L | awk '{printf ( "%d %s %d %s", $(NF - 3), $(NF - 2),  $(NF - 1), $NF) }'`
    END_STR=`date --date="${DURATION}"  +%c`
    SECONDS_LEFT=$(($(date --date="$DURATION" +%s) - $(date +%s)))

    if [ -e ${WS_NAME}_last_email ]
    then
	SECONDS_SINCE_WARNING=$(($(date +%s) - $(cat ${WS_NAME}_last_email)))
    else
        #more than 30 days of seconds
	SECONDS_SINCE_WARNING=3000000
    fi

    if (( $SECONDS_LEFT < $WEEK ))
    then
	printf "    WARNING: "
    else
	printf "             "
    fi

    printf "$WS_NAME expires in $DURATION "
    if ((  ( $SECONDS_LEFT > $WEEK && $SECONDS_SINCE_WARNING > $WEEK ) || 
	   ( $SECONDS_LEFT < $WEEK && $SECONDS_LEFT > $DAY  && $SECONDS_SINCE_WARNING > $DAY ) ||
	   ( $SECONDS_LEFT < $DAY  && $SECONDS_SINCE_WARNING > $HOUR )))
    then
	echo $(date +%s ) >  ${WS_NAME}_last_email

	${MAILER}  ${EMAIL_ADR}  -s "WS DELETE WARNING: $WS_NAME in $DURATION (owned by $USER) "  <<EOF_2
   $USER workspace ${WS_NAME} on host ${HOST}
   will be deleted on: ${END_STR}
EOF_2
	printf "(reminder sent to $EMAIL_ADR)"
    fi

    printf "\n"
done


cd $CURDIR

