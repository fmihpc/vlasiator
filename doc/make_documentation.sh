#!/bin/bash

HOME_DIR=$(realpath ../../)

BUILD_DIR=$(realpath ./_build_doxygen/)


PUBLIC_HTML_DIR=${BUILD_DIR}/public_html
DOXYGEN=doxygen #/home/users/alfthan/doxygen-1.8.1.2/bin/doxygen

VLASIATOR_DIR=${HOME_DIR}/vlasiator/
VLASIATOR_COMMITLOG=${BUILD_DIR}/vlasiator/commitlog.html
VLASIATOR_DOXYGENDOCS=${BUILD_DIR}/vlasiator/doc

mkdir -p $VLASIATOR_DOXYGENDOCS

cd $VLASIATOR_DIR


#create doxygen docs
rm -rf $VLASIATOR_DOXYGENDOCS

git checkout master
git log > $VLASIATOR_COMMITLOG # unformatted git log
sed s+^OUTPUT_DIRECTORY.*+"OUTPUT_DIRECTORY = ${VLASIATOR_DOXYGENDOCS}\n"+ <Doxyfile.template >Doxyfile

$DOXYGEN


git checkout dev

VLASIATOR_COMMITLOG=${BUILD_DIR}/vlasiator-dev/commitlog.html
VLASIATOR_DOXYGENDOCS=${BUILD_DIR}/vlasiator-dev/doc

mkdir -p $VLASIATOR_DOXYGENDOCS
git log > $VLASIATOR_COMMITLOG # unformatted git log
sed s+^OUTPUT_DIRECTORY.*+"OUTPUT_DIRECTORY = ${VLASIATOR_DOXYGENDOCS}\n"+ <Doxyfile.template >Doxyfile

$DOXYGEN



