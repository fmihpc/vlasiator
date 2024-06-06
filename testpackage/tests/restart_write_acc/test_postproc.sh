#!/bin/sh

LAST_RESTART=$(ls -t restart*.vlsv | head -n 1)
test -e $LAST_RESTART && ln -s $LAST_RESTART ./restart.vlsv
