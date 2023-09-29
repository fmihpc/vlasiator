#!/bin/bash

if [[ $# < 1 ]]
then
    echo "Usage: ./submit.sh JOB, where JOB is the job script."
    exit 1
fi

ERR=0

BIN=$(sed -n 's/BIN\s*=\s*//p' < $1)
CFG=$(sed -n 's/CFG\s*=\s*//p' < $1)

if ! [[ -f $BIN ]]
then
   echo "Binary $BIN not found!"
   exit 1
fi

if ! [[ -f $BIN ]]
then
   echo "Config $CFG not found!"
   exit 1
fi

./check_vlasiator_cfg.sh $BIN $CFG > check.txt

if grep -q "Invalid options" check.txt 
then
    echo "Invalid options in config file, check check.txt"
    echo ""
    ERR=1
fi

SW=$(sed -n s/^file_..\\s*=\\s*//p < $CFG)

if ! test -f $SW
then
    echo "Solar wind file $SW not found!"
    echo ""
    ERR=1
fi

MF=$(sed -n s/^atmosphericModelFile\\s*=\\s*//p < $CFG)
if [[ -z $MF ]]
then
    MF=NRLMSIS.dat
fi

if ! test -f $MF
then
    echo "Atmospheric model file $MF not found!"
    echo ""
    ERR=1
fi

IO=$(sed -n s/^system_write_path\\s*=\\s*//p < $CFG)
if [[ -z $IO ]]
then
    IO=.
fi

if ! test -d $IO
then
    echo "System write path $IO not found!"
    echo ""
    ERR=1
else
    echo $IO
    lfs getstripe -d $IO
fi

RW=$(sed -n s/^restart_write_path\\s*=\\s*//p < $CFG)
if [[ -z $RW ]]
then
    RW=.
fi

if ! test -d $RW
then
    echo "Restart write path $RW not found!"
    echo ""
    ERR=1
elif [[ $RW != $IO ]]
then
    echo $RW
    lfs getstripe -d $RW
fi

PROJECT=$(sed -n 's/#SBATCH\s\+--account\s*=\s*//p' < $1)
if which csc-workspaces &> /dev/null
then
    csc-workspaces quota | head -n 1
    csc-workspaces quota | grep $PROJECT
elif which lumi-workspaces &> /dev/null
then
    lumi-workspaces | head -n 5 | tail -n 2
    lumi-workspaces | grep /$PROJECT
    echo ""
    lumi-workspaces | grep -A4 -m1 -e 'Status of'
    lumi-workspaces | grep 
    lumi-workspaces | grep $PROJECT | tail -n 1
fi

if [[ $ERR -gt 0 ]]
then
    echo "Errors found, not submitting job."
    exit 1
else
    echo "Submitting $1"
    sbatch $1
fi
