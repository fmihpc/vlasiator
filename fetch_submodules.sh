#!/bin/bash

set -e   # Abort on error

# Header-only / source dependencies that used to be git submodules.
# Cloned into ./submodules/, matching the paths expected by the Makefile.

mkdir -p submodules
cd submodules

rm -rf fsgrid
git clone https://github.com/fmihpc/fsgrid.git

rm -rf dccrg
git clone -b vlasiator-version https://github.com/fmihpc/dccrg.git

rm -rf eigen
git clone -b master https://gitlab.com/libeigen/eigen.git

rm -rf vectorclass
git clone https://github.com/vectorclass/version2 vectorclass

rm -rf vectorclass-addon
git clone https://github.com/vectorclass/add-on vectorclass-addon

rm -rf hashinator
git clone https://github.com/fmihpc/hashinator.git
