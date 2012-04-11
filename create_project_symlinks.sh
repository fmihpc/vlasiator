#!/bin/bash -e
# Checks whether project symlinks point to given location,
# and updates them if needed.
# The name of the project should be given as an argument
# Removes the old project.o if cpp link wasn't correct

EXPECTED_H=projects/"$1"/"$1".h
EXPECTED_C=projects/"$1"/"$1".cpp

CURRENT_H=$(ls -l project.h | awk '{print $NF}')
CURRENT_C=$(ls -l project.cpp | awk '{print $NF}')

if [ "$EXPECTED_H" != "$CURRENT_H" ]; then
    ln -svf "$EXPECTED_H" project.h
fi

if [ "$EXPECTED_C" != "$CURRENT_C" ]; then
    ln -svf "$EXPECTED_C" project.cpp
    rm -f project.o
fi

