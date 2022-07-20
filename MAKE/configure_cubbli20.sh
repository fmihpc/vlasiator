#!/bin/bash
# This script is assumed to be loacted in MAKE folder

MAKE_FOLDER=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Create external_libs folder

NAME_OF_EXTERNAL_LIBS_FOLDER="external_libs_test"

mkdir -p $MAKE_FOLDER/../$NAME_OF_EXTERNAL_LIBS_FOLDER

EXTERNAL_LIBS_FOLDER=$( cd $(dirname ${MAKE_FOLDER}); cd $MAKE_FOLDER/../$NAME_OF_EXTERNAL_LIBS_FOLDER; pwd -P )

# Change directory in a subshell
(
    cd $EXTERNAL_LIBS_FOLDER
    # Zoltan
    
    # VLSV 
    # DCCRG
    # FSGRID
    # PHIPROF
    # Eigen


    
)