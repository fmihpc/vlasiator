#!/bin/sh

LAST_RESTART=$(ls -t ../restart_write_acc/restart*.vlsv | head -n 1)
test -e $LAST_RESTART && ln -s $LAST_RESTART ./restart.vlsv
